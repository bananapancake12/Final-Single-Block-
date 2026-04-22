! *************************************************************************************************************** !
!
! Last modified by C 12/05/2016
!   - Set up for SHS
!
! TODO
!   - Linear interpolation for interp?...
!
! *************************************************************************************************************** !
!
! Contains: nonlinear  - Called from newribs : Several called functions within file
!                        - Solves the nonlinear advective term
!                        - Records everything to file
!                        - Implementation of immersed boundaries
!           interp_u/v - Called from nonlinear
!                        - Interperlates the grids onto each other to calculate the nonlinear term
!                          - (linear interpolarion)
!           der_x/z    - Called from nonlinear
!                        - Calcualtes the x/z derrivative
!           der_yu/v_h - Called from nonlinear
!                        - Calcualtes the y derrivative
!           dtc_calc   - Called form nonlinear
!                        - Calculates the convective(?) timestep
!           imm_bounds - Called from nonlinear
!                        - Immersed boundaries
!           record_out - Called from nonlinear
!                        - Records to file
!
! Also includes four_to_phys_.. for IFFT
!   phys_to_four_.. for FFT
!   modes_to_planes_.. for modes to planes
!   planes_to_modes_.. for planes to modes
!
!
! *************************************************************************************************************** !

! ****** Modified C 31/08/2015 ****** !
! Advective terms now calculated in conservation form
!	Might not actually be quicker when using immersed boundaries

! module FOU3D_mod
!   use declaration
!   ! use rec_out
!   ! use transpose
!   ! use error_mod
!   ! use spectra_mod
!   implicit none
! contains


subroutine nonlinear(Nu1,Nu2,Nu3,u1,u2,u3,du1,du2,du3,p,div,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!   NONLINEAR TERMS  !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid
  
  integer flagst,flagwr,flagslinst,flagqwr,j,i,k,column
  real(8) C1

  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  complex(8) :: div( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  complex(8) :: Nu1(jlim(1,ugrid)+1     : jlim(2,ugrid)-1,    columns_num(myid) ) 
  complex(8) :: Nu2(jlim(1,vgrid)+1     : jlim(2,vgrid)-1,    columns_num(myid) )
  complex(8) :: Nu3(jlim(1,ugrid)+1     : jlim(2,ugrid)-1,    columns_num(myid) )
  complex(8) :: du1( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: du2( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) ) 
  complex(8) :: du3( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) ) 

  
  if (iter-iter0>=nstat .and. kRK==1) then
    flagst = 1
    iter0  = iter
  else
    flagst = 0
  end if

  if (iter>=iwrite .and. kRK==1) then
    flagwr = 1
    iwrite = iwrite+nwrite
  else
    flagwr = 0
  end if

  if (t>=nextqt) then
    flagqwr = 1
    nextqt = nextqt+10.0d0
  else
    flagqwr = 0
  end if
  
  u1PL_itp = 0d0
  u2PL_itp = 0d0
  u3PL_itp = 0d0

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Interpolating"
  end if

  !C! Interpolate the grid velocities to the other grid points
  call interp_u(u1_itp,u1,myid)
  call interp_v(u2_itp,u2,myid)
  call interp_u(u3_itp,u3,myid) 



  u1PL  = 0d0
  u2PL  = 0d0
  u3PL  = 0d0
  Nu1PL = 0d0
  Nu2PL = 0d0
  Nu3PL = 0d0

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Modes to planes"
  end if

  !C! Shift 6 velocity fields into planes
  !!!!!!!!!!  modes to planes: !!!!!!!!!!
  call modes_to_planes_UVP ( u1PL,    u1,    ugrid,nyu,nyu_LB,myid,status,ierr)
  call modes_to_planes_UVP ( u2PL,    u2,    vgrid,nyv,nyv_LB,myid,status,ierr)
  call modes_to_planes_UVP ( u3PL,    u3,    ugrid,nyu,nyu_LB,myid,status,ierr)
  call modes_to_planes_UVP ( u1PL_itp,u1_itp,vgrid,nyv,nyv_LB,myid,status,ierr)
  call modes_to_planes_UVP ( u2PL_itp,u2_itp,ugrid,nyu,nyu_LB,myid,status,ierr)
  call modes_to_planes_UVP ( u3PL_itp,u3_itp,vgrid,nyv,nyv_LB,myid,status,ierr)


  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Spectra"
  end if

  !!!!!!!!!!!!!   spectra:  !!!!!!!!!!!!!
  if (flagst==1) then
    call spectra(u1,u2,u2_itp,u3,p,myid)
  end if
  

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Record Out"
  end if

  !!!!!!!!!!!!! record out: !!!!!!!!!!!!!
  if (flagwr==1) then
    call error(div,myid,ierr)
    ppPL = 0d0
    call modes_to_planes_UVP(ppPL,p,3,nyp,nyp_LB,myid,status,ierr)
    !call modes_to_planes_UVP(ppPL,div,3,myid,status,ierr) !Output divergence for checking
    call record_out(u1,myid)
  end if

! if (myid == 4) then
  ! ! u1PL
  ! u1PL(:,:,28) = 0
  ! u1PL(:,129,28) = 0
  ! u1PL(15,15,28) = 1
  ! u1PL(39,3,28) = 2
  ! u1PL(67,187,28) = 3
  ! u1PL(122,172,28) = 5
  ! u1PL(1,1,28) = 7
  ! u1PL(128,63,28) = 11


  ! ! u2PL
  ! u2PL(:,:,28) = 0
  ! u2PL(:,129,28) = 0
  ! u2PL(15,15,28) = 1
  ! u2PL(39,3,28) = 2
  ! u2PL(67,187,28) = 3
  ! u2PL(122,172,28) = 5
  ! u2PL(1,1,28) = 7
  ! u2PL(128,63,28) = 11


  ! ! u3PL
  ! u3PL(:,:,28) = 0
  ! u3PL(:,129,28) = 0
  ! u3PL(15,15,28) = 1
  ! u3PL(39,3,28) = 2
  ! u3PL(67,187,28) = 3
  ! u3PL(122,172,28) = 5
  ! u3PL(1,1,28) = 7
  ! u3PL(128,63,28) = 11

  ! ! u1PL_itp
  ! u1PL_itp(:,:,28) = 0
  ! u1PL_itp(:,129,28) = 0
  ! u1PL_itp(15,15,28) = 1
  ! u1PL_itp(39,3,28) = 2
  ! u1PL_itp(67,187,28) = 3
  ! u1PL_itp(122,172,28) = 5
  ! u1PL_itp(1,1,28) = 7
  ! u1PL_itp(128,63,28) = 11

  ! ! u2PL_itp
  ! u2PL_itp(:,:,28) = 0
  ! u2PL_itp(:,129,28) = 0
  ! u2PL_itp(15,15,28) = 1
  ! u2PL_itp(39,3,28) = 2
  ! u2PL_itp(67,187,28) = 3
  ! u2PL_itp(122,172,28) = 5
  ! u2PL_itp(1,1,28) = 7
  ! u2PL_itp(128,63,28) = 11

  ! ! u3PL_itp
  ! u3PL_itp(:,:,28) = 0
  ! u3PL_itp(:,129,28) = 0
  ! u3PL_itp(15,15,28) = 1
  ! u3PL_itp(39,3,28) = 2
  ! u3PL_itp(67,187,28) = 3
  ! u3PL_itp(122,172,28) = 5
  ! u3PL_itp(1,1,28) = 7
  ! u3PL_itp(128,63,28) = 11
! end if

  ! u1PL(:,kLkup(-64),:) = 0
  ! u2PL(:,kLkup(-64),:) = 0
  ! u3PL(:,kLkup(-64),:) = 0
  ! u1PL_itp(:,kLkup(-64),:) = 0
  ! u2PL_itp(:,kLkup(-64),:) = 0
  ! u3PL_itp(:,kLkup(-64),:) = 0

  ! u1PL(iLkup(64):iLkup(64)+1,:,:) = 0
  ! u2PL(iLkup(64):iLkup(64)+1,:,:) = 0
  ! u3PL(iLkup(64):iLkup(64)+1,:,:) = 0
  ! u1PL_itp(iLkup(64):iLkup(64)+1,:,:) = 0
  ! u2PL_itp(iLkup(64):iLkup(64)+1,:,:) = 0
  ! u3PL_itp(iLkup(64):iLkup(64)+1,:,:) = 0

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Ops in Planes"
  end if

  !!!!!!!!! four to ops: !!!!!!!!!
  call ops_in_planes2(myid,flagst) !C! ops in planes to compute velocity products and x/z derriatives

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1, "=====> Planes to modes UVP"
  end if


  !C! Shift y derrivative products to modes
  call planes_to_modes_UVP(uv_f,uv_fPL,vgrid,nyv,nyv_LB,myid,status,ierr)
  call planes_to_modes_UVP(vv_c,vv_cPL,ugrid,nyu,nyu_LB,myid,status,ierr)
  call planes_to_modes_UVP(wv_f,wv_fPL,vgrid,nyv,nyv_LB,myid,status,ierr)
  
  
  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Planes to modes NUVP"
  end if


  !C! Shift x/z derrivatives to modes
  call planes_to_modes_NUVP(Nu1,Nu1PL,ugrid,nyu,nyu_LB,myid,status,ierr)
  call planes_to_modes_NUVP(Nu2,Nu2PL,vgrid,nyv,nyv_LB,myid,status,ierr)
  call planes_to_modes_NUVP(Nu3,Nu3PL,ugrid,nyu,nyu_LB,myid,status,ierr)

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Derivatives"
  end if

  !C! Calculate y derrivatives
    call der_yu_h(Nu1_dy,uv_f,myid)
    call der_yv_h(Nu2_dy,vv_c,myid)
    call der_yu_h(Nu3_dy,wv_f,myid)
  
  !C! Calculate final advective term
    do column = 1,columns_num(myid)
      do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
        Nu1(j,column) = Nu1(j,column)+Nu1_dy(j,column)
        Nu3(j,column) = Nu3(j,column)+Nu3_dy(j,column)
      end do
      do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
        Nu2(j,column) = Nu2(j,column)+Nu2_dy(j,column)
      enddo
    enddo
  ! enddo   

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> CFL and Stats"
  end if


  !!!!!!!!!!!  CFL and stats: !!!!!!!!!!!
  if (kRK==1) then
  ! Calculating timestep
    call dtc_calc(u1PL,u2PL,u3PL,myid)
    call MPI_ALLREDUCE(dt,dtc,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)

    ! write(6,*) " Finished MPI reduce =====> dt and dti"

    dt  = min(2.5d0*dtv,CFL*dtc)
    ! write(6,*) "dt"

    dti = 1d0/dt
    ! write(6,*) "dti"

    t   = t+dt
    ! write(6,*) "t+dt"

    ! Calculating and writing stats and spectra
    if (flagst==1) then

    ! if(myid==0) then
    !   write(6,*) "t=", MPI_Wtime() - t1," CFL and Stats =====> calc omega"
    ! end if
    
    !C! Calculate Omega_x
    u2PLN=0d0
    call modes_to_planes_phys (u2PLN,u2,vgrid,nyv,nyv_LB,myid,status,ierr)    
    
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call der_z_N(u2PLN(:,:,j),wx(:,:,j),k1F_z) !C! u2PL in Fourier space
      ! call der_z_N(u2PLN(:,:,j),wx(:,:,j),k1F_z) !C! u2PL in Fourier space
      if(j == 10 ) then
        !write(6,*) "wx", wx(:,10,j)
      end if 


    end do
  
    call der_yv_h_wx(du3dy_columns,u3,myid)

    call modes_to_planes_phys(du3dy_planes,du3dy_columns,vgrid,nyv,nyv_LB,myid,status,ierr)
    
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      do i = 1,Nspec_x+2
        do k = 1,Nspec_z
          wx(i,k,j) = -wx(i,k,j)+du3dy_planes(i,k,j) !omega_x at faces!
        end do
      end do
      call four_to_phys_N(wx(1,1,j))
    end do

      call stats(myid,status,ierr) 
      ! (u1,u3,myid,status,ierr)
      
    end if
    
    if (flagwr==1) then
      call write_spect(myid,status,ierr) 
      call write_stats(myid,status,ierr) 
      ! call write_sl_stats(myid,status,ierr) 
    end if    
    
    ! if (flagslinst==1) then
    !   call inst_sl_stats(u1,u3,myid,status,ierr)
    ! endif
        
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Update RHS with new advective term. du now full RHS
  C1 = -gRK(kRK)
  do column = 1,columns_num(myid)
    do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      du1(j,column) = du1(j,column)+C1*Nu1(j,column)
      du3(j,column) = du3(j,column)+C1*Nu3(j,column)
    enddo
    do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
      du2(j,column) = du2(j,column)+C1*Nu2(j,column)
    enddo
  enddo
  ! enddo

end subroutine

subroutine interp_u(u_itp,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  interp u/w grid to v grid  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none
  
  integer j,column,myid
  complex(8)  u( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )
  complex(8)  u_itp( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) )
  
  do column = 1,columns_num(myid)
    ! We interpolate everything. vgrid has got one less point than ugrid
    do j = jlim(1,vgrid),jlim(2,vgrid)  
      u_itp(j,column) = ((yv(j)-yu(j))*u(j+1,column)+(yu(j+1)-yv(j))*u(j,column))/(yu(j+1)-yu(j))
    end do
  end do

  
end subroutine

subroutine interp_v(u_itp,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  interp v grid to u/w grid  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none
  
  integer j,column,myid
  complex(8)   u  ( jlim(1,vgrid) : jlim(2,vgrid), columns_num(myid) )
  complex(8)   u_itp( jlim(1,ugrid):jlim(2,ugrid), columns_num(myid) )

  do column = 1,columns_num(myid)
    do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      u_itp(j,column) = ((yu(j)-yv(j-1))*u(j,column)+(yv(j)-yu(j))*u(j-1,column))/(yv(j)-yv(j-1))
    end do
    u_itp(jlim(1,ugrid),column) = &
&      gridweighting_interp(1)*(u(jlim(1,vgrid),column)-u_itp(jlim(1,ugrid)+1,column)) &
&      + u_itp(jlim(1,ugrid)+1,column)
    u_itp(jlim(2,ugrid),column) = &
&      gridweighting_interp(2)*(u(jlim(2,vgrid),column)-u_itp(jlim(2,ugrid)-1,column)) &
&      + u_itp(jlim(2,ugrid)-1,column)
  end do

                  
end subroutine

subroutine der_x(u,dudx,kx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  use declaration
  implicit none

  integer i,k
  complex(8) u   (0:Ngal_x/2, Ngal_z)
  complex(8) dudx(0:Ngal_x/2, Ngal_z)
  complex(8) kx(0:Nspec_x/2)

  do k = 1,Nspec_z/2
    do i = 0,Nspec_x/2
      dudx(i,k) = kx(i)*u(i,k)
    end do
    do i = Nspec_x/2+1,Ngal_x/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudx(i,k) = 0d0                               !!!!!!!!!!!!  must be explicitly specified
    end do
  end do

  !Zeros includes mode Nz/2+1 mode (zero for advection derrivatives)
  do k = Nspec_z/2+1,Ngal_z-Nspec_z/2+1!!!!!!!!!!!!  Zeros for the antialiasing region (z-dir)
    do i = 0,Ngal_x/2                        !!!!!!!!!!!!
      dudx(i,k) = 0d0                               !!!!!!!!!!!!  must be explicitly specified
    end do
  end do
  do k = Ngal_z-Nspec_z/2+2,Ngal_z
    do i = 0,Nspec_x/2
      dudx(i,k) = kx(i )*u(i,k)
    end do
    do i = Nspec_x/2+1,Ngal_x/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudx(i,k) = 0d0                               !!!!!!!!!!!!  must be explicitly specified
    end do
  end do

end subroutine


subroutine der_x_N(u,dudx,kx) !For N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer i,k,dk2,k2
  complex(8) u   (0:Nspec_x/2,Nspec_z)
  complex(8) dudx(0:Nspec_x/2,Nspec_z)
  complex(8) kx(0:Nspec_x/2)

  do k = 1,Nspec_z
    do i = 0,Nspec_x/2
      dudx(i,k) = kx(i)*u(i,k)
    end do
  end do
  
end subroutine


subroutine der_z(u,dudz,kz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  use declaration
  implicit none

  integer i,k,dk2,k2
  complex(8) u   (0:Ngal_x/2,Ngal_z)
  complex(8) dudz(0:Ngal_x/2,Ngal_z)
  complex(8) kz(1:Nspec_z)

  dk2=Nspec_z-Ngal_z

  do k = 1,Nspec_z/2
    do i = 0,Nspec_x/2
      dudz(i,k) = kz(k)*u(i,k)
    end do
    do i = Nspec_x/2+1,Ngal_x/2       !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudz(i,k) = 0d0
    end do
  end do
  !Zeros includes mode Nz/2+1 mode (zero for advection derrivatives)
  do k = Nspec_z/2+1,Ngal_z-Nspec_z/2+1!!!!!!!!!!!!  Zeros for the antialiasing region (z-dir)
    do i = 0,Ngal_x/2 
      dudz(i,k) = 0d0
    end do
  end do
  do k = Ngal_z-Nspec_z/2+2,Ngal_z
    k2 = k+dk2
    do i = 0,Nspec_x/2
      dudz(i,k) = kz(k2)*u(i,k)
    end do
    do i = Nspec_x/2+1,Ngal_x/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudz(i,k) = 0d0
    end do
  end do

end subroutine

subroutine der_z_N(u,dudz,kz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  use declaration
  implicit none

  integer i,k,dk2,k2
  complex(8) u   (0:Nspec_x/2,Nspec_z)
  complex(8) dudz(0:Nspec_x/2,Nspec_z) 
  complex(8) kz(1:Nspec_z)

  ! doing the operation in complex space 
    do k = 1,Nspec_z
      do i = 0,Nspec_x/2
        dudz(i,k) = kz(k)*u(i,k)
      end do
    end do


end subroutine


subroutine der_yu_h(dudy,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      f-->c     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,myid
  complex(8)  u   (jlim(1,vgrid)  :jlim(2,vgrid)  ,columns_num(myid))
  complex(8) dudy (jlim(1,ugrid)+1:jlim(2,ugrid)-1,columns_num(myid))

  do column = 1,columns_num(myid)
    do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      dudy(j,column) = (u(j,column)-u(j-1,column))*dthdyu(j)*ddthetavi
    end do
  end do

end subroutine

subroutine der_yv_h(dudy,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      c-->f     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,myid
  complex(8)  u   (jlim(1,ugrid)  :jlim(2,ugrid)  ,columns_num(myid))
  complex(8) dudy (jlim(1,vgrid)+1:jlim(2,vgrid)-1,columns_num(myid))

  do column = 1,columns_num(myid)
    do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
      dudy(j,column) = (u(j+1,column)-u(j,column))*dthdyv(j)*ddthetavi
    end do
  end do

end subroutine

subroutine der_yv_h_wx(dudy,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      c-->f     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,myid
  complex(8)  u   (jlim(1,ugrid):jlim(2,ugrid),columns_num(myid))
  complex(8) dudy (jlim(1,vgrid):jlim(2,vgrid),columns_num(myid))

  do column = 1,columns_num(myid)
    do j = jlim(1,vgrid),jlim(2,vgrid)
      dudy(j,column) = (u(j+1,column)-u(j,column))*dthdyv(j)*ddthetavi
    end do
  end do

end subroutine

subroutine dtc_calc(u1,u2,u3,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!       CFL      !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,myid
  real(8) u1mloc,u2mloc,u3mloc
  real(8) u1(igal,kgal,jgal(2,1):jgal(2,2))
  real(8) u2(igal,kgal,jgal(1,1):jgal(1,2))
  real(8) u3(igal,kgal,jgal(2,1):jgal(2,2))

  u1mloc = maxval(abs(u1))
  u2mloc = 1d-6

  do j = jgal(1,1),jgal(1,2)
    do k = 1,kgal
      do i = 1,igal-2
        u2mloc = max(u2mloc,abs(u2(i,k,j))*dthdyv(j))
      end do
    end do
  end do

  u3mloc = maxval(abs(u3))

  ! write(6,*) "u2max", u2mloc
  ! write(6,*) "dthetavi", dthetavi
  
  dt = min(1d0/(alp*(Nspec_x/2)*u1mloc),1d0/(bet*(Nspec_z/2)*u3mloc),1d0/(u2mloc*dthetavi))

  ! write(6,*) "u1max", u1mloc
  ! write(6,*) "u2max", u2mloc
  ! write(6,*) "u3max", u3mloc

  ! write(6,*) "cfl u1", 1d0/(alp*(Nspec_x/2)*u1mloc), myid
  ! write(6,*) "=====> finished DTC calc"

end subroutine

subroutine four_to_phys_u(u1,u2,u3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Fourier Space to Physical Space u and its derivatives
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none

  real(8) u1 (Ngal_x+2,Ngal_z)
  real(8) u2 (Ngal_x+2,Ngal_z)
  real(8) u3 (Ngal_x+2,Ngal_z)
  
  u1(:,Ngal_z/2+1)=0d0
  u2(:,Ngal_z/2+1)=0d0
  u3(:,Ngal_z/2+1)=0d0
  
  call cft(u1,Ngal_x+2,2,(Nspec_x+2)/2,1,buffCal_z)
  call rft(u1,Ngal_x+2,Ngal_z,1,buffRal_x)
  call cft(u2,Ngal_x+2,2,(Nspec_x+2)/2,1,buffCal_z)
  call rft(u2,Ngal_x+2,Ngal_z,1,buffRal_x)
  call cft(u3,Ngal_x+2,2,(Nspec_x+2)/2,1,buffCal_z)
  call rft(u3,Ngal_x+2,Ngal_z,1,buffRal_x)
  
end subroutine

subroutine four_to_phys_du(du)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Fourier Space to Physical Space u and its derivatives
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none
  real(8) du  (Ngal_x+2,Ngal_z)

  du(:,Ngal_z/2+1)=0d0
  
  call cft(du   ,Ngal_x+2,2,(Nspec_x+2)/2,1,buffCal_z)
  call rft(du   ,Ngal_x+2,Ngal_z,1,buffRal_x)

end subroutine

subroutine four_to_phys_N(du)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Fourier Space to Physical Space u and its derivatives
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none

  real(8) du  (Nspec_x+2,Nspec_z)
  call cft(du   ,Nspec_x+2,2,(Nspec_x+2)/2,1,buffC_z)
  call rft(du   ,Nspec_x+2,Nspec_z,1,buffR_x)

end subroutine


subroutine modes_to_planes_UVP (xPL,x,grid,nygrid,nygrid_LB,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,grid
  integer column
  integer inode,yourid
  integer msizeR,msizeS
  complex(8)   x(jlim(1,grid):jlim(2,grid),columns_num(myid))
  real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1)
  integer :: nygrid, nygrid_LB
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  ! Loop for itself
  ! Transpose the cube that it already owns

  yourid = myid

  jminR = max(planelim(grid,1,  myid),jlim(1,grid)+1)
  jmaxR = min(planelim(grid,2,  myid),jlim(2,grid)-1)

  if (jminR==nygrid_LB+1 .and. jmaxR>=jminR) then
    jminR = jminR-1
  end if
  if (jmaxR==nygrid   .and. jmaxR>=jminR) then
    jmaxR = jmaxR+1
  end if
  do j = jminR,jmaxR
    do column = 1,columns_num(yourid)
      i = columns_i(column,yourid)
      k = columns_k(column,yourid) - dk(column,yourid)
      xPL(2*i+1,k,j) = dreal(x(j,column))
      xPL(2*i+2,k,j) = dimag(x(j,column))
    end do
  end do


  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then

      jminS = max(planelim(grid,1,yourid),jlim(1,grid)+1)
      jmaxS = min(planelim(grid,2,yourid),jlim(2,grid)-1)
      jminR = max(planelim(grid,1,  myid),jlim(1,grid)+1)
      jmaxR = min(planelim(grid,2,  myid),jlim(2,grid)-1)
      if (jminS==nygrid_LB+1  ) then
        jminS = jminS-1
      end if
      if (jmaxS==nygrid) then
        jmaxS = jmaxS+1
      end if
      if (jminR==nygrid_LB+1  ) then
        jminR = jminR-1
      end if
      if (jmaxR==nygrid) then
        jmaxR = jmaxR+1
      end if
      allocate(buffS(jminS:jmaxS,columns_num(  myid)))
      allocate(buffR(jminR:jmaxR,columns_num(yourid)))
      msizeS = 2*(columns_num(  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
      msizeR = 2*(columns_num(yourid)*(jmaxR-jminR+1))
      msizeS = max(msizeS,0)
      msizeR = max(msizeR,0)

      do j = jminS,jmaxS
        do column = 1,columns_num(myid)
          buffS(j,column) = x(j,column)

        end do
      end do

      call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
&                         MPI_COMM_WORLD,status,ierr)
      do j = jminR,jmaxR
        do column = 1,columns_num(yourid)
          i = columns_i(column,yourid)
          k = columns_k(column,yourid) - dk(column,yourid)
          xPL(2*i+1,k,j) = dreal(buffR(j,column))
          xPL(2*i+2,k,j) = dimag(buffR(j,column))
        end do
      end do

      deallocate(buffR,buffS)

      ! end do
    end if
  end do

end subroutine

subroutine modes_to_planes_phys (xPL,x,grid,nygrid,nygrid_LB,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,grid !myiband
  integer column

  integer inode,yourid
  integer msizeR,msizeS
  ! complex(8), intent(in) :: x(jlim(1,grid):,:)
  complex(8)   x( jlim(1,grid) : jlim(2,grid), columns_num(myid) )
  real(8)      xPL(Nspec_x+2,Nspec_z,jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable :: buffS(:,:),buffR(:,:)
  integer, intent(in) :: nygrid, nygrid_LB


  yourid = myid
  jminR = max(planelim(grid,1,  myid),jlim(1,grid)+1)
  jmaxR = min(planelim(grid,2,  myid),jlim(2,grid)-1)
  if (jminR==nygrid_LB+1  ) then
        jminR = jminR-1
      end if
      if (jmaxR==nygrid) then
        jmaxR = jmaxR+1
      end if
  do j = jminR,jmaxR
    do column = 1,columns_num(yourid)
      i = columns_i(column,yourid)
      k = columns_k(column,yourid) - dk_phys(column,yourid)


      ! write(6,*) "k", k, "columns_k",columns_k(column,yourid), "dk_phys", dk_phys(column,yourid), yourid


      xPL(2*i+1,k,j) = dreal(x(j,column))
      xPL(2*i+2,k,j) = dimag(x(j,column))
    end do
  end do
  ! end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      jminS = max(planelim(grid,1,yourid),jlim(1,grid)+1)
      jmaxS = min(planelim(grid,2,yourid),jlim(2,grid)-1)
      jminR = max(planelim(grid,1,  myid),jlim(1,grid)+1)
      jmaxR = min(planelim(grid,2,  myid),jlim(2,grid)-1)
      if (jminS==nygrid_LB+1  ) then
        jminS = jminS-1
      end if
      if (jmaxS==nygrid) then
        jmaxS = jmaxS+1
      end if
      if (jminR==nygrid_LB+1  ) then
        jminR = jminR-1
      end if
      if (jmaxR==nygrid) then
        jmaxR = jmaxR+1
      end if
      allocate(buffS(jminS:jmaxS,columns_num(  myid)))
      allocate(buffR(jminR:jmaxR,columns_num(yourid)))
      msizeS = 2*(columns_num(  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
      msizeR = 2*(columns_num(yourid)*(jmaxR-jminR+1))
      msizeS = max(msizeS,0)
      msizeR = max(msizeR,0)

      do j = jminS,jmaxS
        do column = 1,columns_num(myid)
          buffS(j,column) = x(j,column)

        end do
      end do

      call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
&                         MPI_COMM_WORLD,status,ierr)
      do j = jminR,jmaxR
        do column = 1,columns_num(yourid)
          i = columns_i(column,yourid)
          k = columns_k(column,yourid) - dk_phys(column,yourid)
          xPL(2*i+1,k,j) = dreal(buffR(j,column))
          xPL(2*i+2,k,j) = dimag(buffR(j,column))
        end do
      end do

      deallocate(buffR,buffS)

      ! end do
    end if
  end do

end subroutine

! subroutine modes_to_planes_phys_lims (xPL,x,nystart,nyend,grid,nygrid,nygrid_LB,myid,myiband,status,ierr)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!! MODES TO PLANES USED IN FFT TRID LU!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   use declaration
!   implicit none

!   include 'mpif.h'             ! MPI variables
!   integer status(MPI_STATUS_SIZE),ierr,myid

!   integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband,grid,myiband
!   integer column,nystart,nyend

!   integer inode,yourid
!   integer msizeR,msizeS
!   type(cfield) x  
!   real(8)      xPL(Nspec_x+2,Nspec_z,limPL_FFT(grid,1,myid):limPL_FFT(grid,2,myid))
!   complex(8), allocatable :: buffS(:,:),buffR(:,:)

!   plband = bandPL_FFT(myid)
!   yourid = myid
!   ! do iband = sband,eband
!     ! jband = iband
! !     jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)+1),nystart)
! !     jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)-1),nyend)
! jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid)),nystart)
! jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid)),nyend)
! !     if (jminR==Ny(grid,0)+1  ) then
! !           jminR = jminR-1
! !         end if
! !         if (jmaxR==Ny(grid,nband)) then
! !           jmaxR = jmaxR+1
! !         end if
!     do j = jminR,jmaxR
!       do column = 1,columns_num(yourid)
!         i = columns_i(column,yourid)
!         k = columns_k(column,yourid) - dk_phys(column,yourid)
!         xPL(2*i+1,k,j) = dreal(x%f(j,column))
!         xPL(2*i+2,k,j) = dimag(x%f(j,column))
!       end do
!     end do
!   ! end do

!   do inode = 1,pnodes-1
!     yourid = ieor(myid,inode)
!     if (yourid<np) then
!       ! do iband = sband,eband
!         !jband=crossband(iband,yourid)
!         ! jband = iband
! !         jminS = max(max(limPL_FFT(grid,1,yourid),jlim(1,grid,iband)+1),nystart)
! !         jmaxS = min(min(limPL_FFT(grid,2,yourid),jlim(2,grid,iband)-1),nyend)
! !         jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)+1),nystart)
! !         jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)-1),nyend)
! jminS = max(max(limPL_FFT(grid,1,yourid),jlim(1,grid)),nystart)
! jmaxS = min(min(limPL_FFT(grid,2,yourid),jlim(2,grid)),nyend)
! jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid)),nystart)
! jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid)),nyend)
! !         if (jminS==Ny(grid,0)+1  ) then
! !           jminS = jminS-1
! !         end if
! !         if (jmaxS==Ny(grid,nband)) then
! !           jmaxS = jmaxS+1
! !         end if
! !         if (jminR==Ny(grid,0)+1  ) then
! !           jminR = jminR-1
! !         end if
! !         if (jmaxR==Ny(grid,nband)) then
! !           jmaxR = jmaxR+1
! !         end if
!         allocate(buffS(jminS:jmaxS,columns_num(  myid)))
!         allocate(buffR(jminR:jmaxR,columns_num(yourid)))
!         msizeS = 2*(columns_num(  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
!         msizeR = 2*(columns_num(yourid)*(jmaxR-jminR+1))
!         msizeS = max(msizeS,0)
!         msizeR = max(msizeR,0)

!         do j = jminS,jmaxS
!           do column = 1,columns_num(myid)
!             buffS(j,column) = x%f(j,column)

!           end do
!         end do

!         call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*nband+11*nband, &
! &                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*nband+7*nband, &
! &                         MPI_COMM_WORLD,status,ierr)
!         do j = jminR,jmaxR
!           do column = 1,columns_num(yourid)
!             i = columns_i(column,yourid)
!             k = columns_k(column,yourid) - dk_phys(column, yourid)
!             xPL(2*i+1,k,j) = dreal(buffR(j,column))
!             xPL(2*i+2,k,j) = dimag(buffR(j,column))
!           end do
!         end do

!         deallocate(buffR,buffS)

!       ! end do
!     end if
!   end do

! end subroutine

! subroutine modes_to_planes_phys_lims_2 (xPL,x,nystart,nyend,grid,nygrid,nygrid_LB,myid,status,ierr)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! used in slip stats so not really needed


!   use declaration
!   implicit none

!   include 'mpif.h'             ! MPI variables
!   integer status(MPI_STATUS_SIZE),ierr,myid

!   integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,grid
!   integer column,nystart,nyend
!   integer inode,yourid
!   integer msizeR,msizeS
!   type(cfield) x
!   !complex(8), intent(in) :: x(jlim(1,grid):,:)
!   real(8)      xPL(Nspec_x+2,Nspec_z,jgal(grid,1)-1:jgal(grid,2)+1)
!   complex(8), allocatable :: buffS(:,:),buffR(:,:)
!   integer, intent(in) :: nygrid,nygrid_LB 
  
!   yourid = myid

!     jminR = max(max(limPL_incw(grid,1,  myid),jlim(1,grid)+1),nystart)
!     jmaxR = min(min(limPL_incw(grid,2,  myid),jlim(2,grid)-1),nyend)
!     if (jminR==nygrid_LB+1  ) then
!           jminR = jminR-1
!         end if
!         if (jmaxR==nygrid) then
!           jmaxR = jmaxR+1
!         end if
        
!     do j = jminR,jmaxR
!       do column = 1,columns_num(yourid)
!         i = columns_i(column,yourid)
!         k = columns_k(column,yourid) - dk_phys(column,yourid)
!         xPL(2*i+1,k,j) = dreal(x%f(j,column))
!         xPL(2*i+2,k,j) = dimag(x%f(j,column))
!       end do
!     end do

!   do inode = 1,pnodes-1
!     yourid = ieor(myid,inode)
!     if (yourid<np) then
!         jminS = max(max(limPL_incw(grid,1,yourid),jlim(1,grid)+1),nystart)
!         jmaxS = min(min(limPL_incw(grid,2,yourid),jlim(2,grid)-1),nyend)
!         jminR = max(max(limPL_incw(grid,1,  myid),jlim(1,grid)+1),nystart)
!         jmaxR = min(min(limPL_incw(grid,2,  myid),jlim(2,grid)-1),nyend)
!         if (jminS==nygrid_LB+1  ) then
!           jminS = jminS-1
!         end if
!         if (jmaxS==nygrid) then
!           jmaxS = jmaxS+1
!         end if
!         if (jminR==nygrid_LB+1  ) then
!           jminR = jminR-1
!         end if
!         if (jmaxR==nygrid) then
!           jmaxR = jmaxR+1
!         end if
                
!         allocate(buffS(jminS:jmaxS,columns_num(  myid)))
!         allocate(buffR(jminR:jmaxR,columns_num(yourid)))
!         msizeS = 2*(columns_num(  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
!         msizeR = 2*(columns_num(yourid)*(jmaxR-jminR+1))
!         msizeS = max(msizeS,0)
!         msizeR = max(msizeR,0)

!         do j = jminS,jmaxS
!           do column = 1,columns_num(myid)
!             buffS(j,column) = x%f(j,column)
!           end do
!         end do

!         call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &
! &                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
! &                         MPI_COMM_WORLD,status,ierr)

!         do j = jminR,jmaxR
!           do column = 1,columns_num(yourid)
!             i = columns_i(column,yourid)
!             k = columns_k(column,yourid) - dk_phys(column,yourid)
!             xPL(2*i+1,k,j) = dreal(buffR(j,column))
!             xPL(2*i+2,k,j) = dimag(buffR(j,column))
!           end do
!         end do

!         deallocate(buffR,buffS)

!       ! end do
!     end if
!   end do

! end subroutine

subroutine ops_in_planes(myid,flagst)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! OPS IN PLANES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  include 'mpif.h'
  
  integer :: myid,flagst
  integer i,k,j,l,temp,jidx
  integer ia, ip, ka, kp, la, lp
  real(8) ddyi
  real(8), allocatable:: du1dx(:,:),du1dz(:,:),du2dx(:,:),du2dz(:,:),du3dx(:,:),du3dz(:,:), buff(:,:)
  integer :: rank_true, nprocs_true, ierr_local
  real(8) :: t0

  allocate(du1dx(igal,kgal),du1dz(igal,kgal))
  allocate(du2dx(igal,kgal),du2dz(igal,kgal))
  allocate(du3dx(igal,kgal),du3dz(igal,kgal))
  allocate(buff(igal,kgal))

  t0 = MPI_Wtime()
  
  du1dx = 0d0
  du1dz = 0d0
  du2dx = 0d0
  du2dz = 0d0
  du3dx = 0d0
  du3dz = 0d0

  uu_cPL = 0d0
  uw_cPL = 0d0
  vv_cPL = 0d0
  wu_cPL = 0d0
  ww_cPL = 0d0




  do j = limPL_excw(ugrid,1,myid),limPL_excw(ugrid,2,myid)
    ! nonlinear interaction into 0th mode (bad, much easier to do if multiply by 4 and do one quadrant but need to think about k = 0 terms carefully)
    do i = -(Nspec_x/2),Nspec_x/2
      do k = -(Nspec_z/2),Nspec_z/2-1
        la = i
        lp = -i
        ia = iLkup(la)
        ip = iLkup(lp)
        ka = kLkup(iNeg(la)*k)
        kp = kLkup(iNeg(lp)*(-k))
        vv_cPL(1,1,j) = vv_cPL(1,1,j) + u2PL_itp(ia,ka,j)*u2PL_itp(ip,kp,j) - &
              & iNeg(i)*iNeg(-i)*u2PL_itp(ia+1,ka,j)*u2PL_itp(ip+1,kp,j)
        vv_cPL(2,1,j) = vv_cPL(2,1,j) + u2PL_itp(ia,ka,j)*u2PL_itp(ip+1,kp,j)*iNeg(-i) + &
              & iNeg(i)*u2PL_itp(ia+1,ka,j)*u2PL_itp(ip,kp,j)
        ! really the imaginary part should be 0 so no need to compute but anyway just in case
      end do
    end do

    ! if(j==10) then
    !   write(6,*) " Finished 0th mode nonlin =====> Linear advection", myid
    !   write(6,*) "vv_cPL", vv_cPL(1,1,j), vv_cPL(2,1,j), j 
    ! end if 

    
    ! if(j==limPL_excw(ugrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1," Finished 0th mode nonlin =====> Linear advection", myid
    end if 


    ! if(j==150) then
    !   write(6,*) "vv_cPL", vv_cPL(1,1,j), vv_cPL(2,1,j), j 
    ! end if 
    
    
    ! linear advection
    ! wrong for 0th mode, 0,0 interaction should not be counted twice, but doesn't matter since differentiate = 0
    ! vv_cPL correct since 0th mode computed earlier
    buff(:,:) = 2*(u1PL(1,1,j)*u1PL(:,:,j))
    uu_cPL(:,:,j) = uu_cPL(:,:,j) + buff(:,:)
    
    ! do k = 1,1
    !   do i = 1, igal
    !     write(6,*) "i", i, "uu_cPL", uu_cPL(i,k,1)
    !   end do 
    ! end do 


    buff(:,:) = 2*(u3PL(1,1,j)*u3PL(:,:,j))
    ww_cPL(:,:,j) = ww_cPL(:,:,j) + buff(:,:)

    buff(:,:) = 2*(u2PL_itp(1,1,j)*u2PL_itp(:,:,j))
    buff(1:2,1) = 0
    vv_cPL(:,:,j) = vv_cPL(:,:,j) + buff(:,:)
    
    buff(:,:) = u1PL(1,1,j)*u3PL(:,:,j) + u3PL(1,1,j)*u1PL(:,:,j)  
    wu_cPL(:,:,j) = wu_cPL(:,:,j) + buff(:,:)
    uw_cPL(:,:,j) = uw_cPL(:,:,j) + buff(:,:)



    ! if (j == 28) then
    !   write(ext4,'(i5.5)') int(1000d0*(t-700)+kRK)
    !   fnameima = 'output/20_'//ext4//'.dat'
    !   write(*,*) fnameima
    !   open(11,file=fnameima,form='unformatted')
    !   write(11) uu_cPL(:,:,j)
    !   write(11) vv_cPL(:,:,j)
    !   write(11) ww_cPL(:,:,j)
    !   write(11) uw_cPL(:,:,j)
    !   write(*,*) u1PL(1,1,j)
    !   write(11) wu_cPL(:,:,j)
    ! end if  

    ! if (j == 28) then
    !   do i= 1,10
    !     write(6,*) "uu_cPL", uu_cPL(i,1,22)
    !   end do 
    ! end if 


    ! nonlinear advection: go through a list
    ! fields = {'uu', 'uv', 'uw', 'vu', 'vv', 'vw', 'wu', 'wv', 'ww'}; in the order of 1 to 9, where the first is the passive

    ! write(6,*) "=====> Non Linear advection", myid

    if (j > (nyv+1)/2) then
      jidx = j-1
    else 
      jidx = j
    end if

    ! if(j==10) then
    !   write(6,*) "=====> nonlinear interaction into 0th mode - U", myid
    ! end if 

    !if(j==limPL_excw(ugrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1,"=====> nonlinear interaction into 0th mode - U", myid
    end if 

    

    call nonlinInter(nonlin(jidx, 1), uu_cPL(1,1,j), u1PL(1,1,j), u1PL(1,1,j))
    call nonlinInter(nonlin(jidx, 3), uw_cPL(1,1,j), u3PL(1,1,j), u1PL(1,1,j))
    call nonlinInter(nonlin(jidx, 5), vv_cPL(1,1,j), u2PL_itp(1,1,j), u2PL_itp(1,1,j))
    call nonlinInter(nonlin(jidx, 7), wu_cPL(1,1,j), u1PL(1,1,j), u3PL(1,1,j))
    call nonlinInter(nonlin(jidx, 9), ww_cPL(1,1,j), u3PL(1,1,j), u3PL(1,1,j))

    ! do k = 1,1
    !   do i = 1, igal
    !     write(6,*) "i", i, "uu_cPL", uu_cPL(i,k,1)
    !   end do 
    ! end do 

    ! if(j==10) then
    !   write(6,*) "=====> Finished Nonlin Inter U", myid
    ! end if 

    ! if(j==limPL_excw(ugrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1,"=====> Finished Nonlin Inter U", myid
    end if 

    ! if (j == 28) then
    !   write(ext4,'(i5.5)') int(1000d0*(t-700)+kRK)
    !   fnameima = 'output/20_'//ext4//'.dat'
    !   write(*,*) fnameima
    !   open(11,file=fnameima,form='unformatted')
    !   write(11) uu_cPL(:,:,j)
    !   write(11) vv_cPL(:,:,j)
    !   write(11) ww_cPL(:,:,j)
    !   write(11) uw_cPL(:,:,j)
    !   ! write(*,*) u1PL(1,1,j)
    !   ! write(11) wu_cPL(:,:,j)
    ! end if  
    ! write(6,*) "k1F_x", k1F_x, "k1F_z", k1F_z

    ! write(6,*) "uw_cPL size   =", size(uw_cPL,1), size(uw_cPL,2), size(uw_cPL,3)
    ! write(6,*) "igal", igal, "kgal", kgal

    ! differentiate in x and z
    call der_x(uu_cPL(1,1,j),du1dx,k1F_x)
    call der_z(uw_cPL(1,1,j),du1dz,k1F_z)
    call der_x(wu_cPL(1,1,j),du3dx,k1F_x)
    call der_z(ww_cPL(1,1,j),du3dz,k1F_z)

    ! if(j==limPL_excw(ugrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1, "=====> Finished Derivatives -U", myid
    end if 
  
    ! if (j == 28) then
    !   write(ext4,'(i5.5)') int(1000d0*(t-700)+kRK)
    !   fnameima = 'output/30_'//ext4//'.dat'
    !   write(*,*) fnameima
    !   open(11,file=fnameima,form='unformatted')
    !   write(11) du1dx
    !   write(11) du1dz
    !   write(11) du3dx
    !   write(11) du3dz
    !   ! write(11) wu_cPL(:,:,j)
    ! end if  
    ! write(6,*) "j=", j, myid

    do k = 1,Ngal_z
      do i = 1,Ngal_x
        Nu1PL(i,k,j) = du1dx(i,k)+du1dz(i,k)
        Nu3PL(i,k,j) = du3dx(i,k)+du3dz(i,k)
      end do
    end do
    ! if (j == 28) then
    !   write(ext4,'(i5.5)') int(1000d0*(t-700)+kRK)
    !   fnameima = 'output/10_'//ext4//'.dat'
    !   write(*,*) fnameima
    !   open(11,file=fnameima,form='unformatted')
    !   write(11) Nu1PL(:,:,j)
    !   write(11) Nu3PL(:,:,j)
    ! end if
    call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j))

    ! if(j==limPL_excw(ugrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1, "=====> Finished U", myid
    end if 

  end do
    
  uv_fPL = 0d0
  vu_fPL = 0d0
  vw_fPL = 0d0
  wv_fPL = 0d0

  ! if(myid==0) then
  !   write(6,*) "=====> nonlinear interaction into 0th mode - V", myid
  ! end if 

  ! if(j==limPL_excw(ugrid,2,myid)) then
  if(j==10) then
    write(6,*) "t=", MPI_Wtime() - t1, "=====> Begining V", myid
  end if 

  do j = limPL_excw(vgrid,1,myid),limPL_excw(vgrid,2,myid)
    ! nonlinear interaction into 0th mode
    do i = -(Nspec_x/2),Nspec_x/2
      do k = -(Nspec_z/2-1),Nspec_z/2-1
        la = i
        lp = -i
        ia = iLkup(la)
        ip = iLkup(lp)
        ka = kLkup(iNeg(la)*k)
        kp = kLkup(iNeg(lp)*(-k))

        uv_fPL(1,1,j) = uv_fPL(1,1,j) + u1PL_itp(ia,ka,j)*u2PL(ip,kp,j) - &
              & iNeg(i)*iNeg(-i)*u1PL_itp(ia+1,ka,j)*u2PL(ip+1,kp,j)
        uv_fPL(2,1,j) = uv_fPL(2,1,j) + u1PL_itp(ia,ka,j)*u2PL(ip+1,kp,j)*iNeg(-i) + &
              & iNeg(i)*u1PL_itp(ia+1,ka,j)*u2PL(ip,kp,j)

        ! write(6,*) 'myid=', myid, ' j=', j, ' uv_fPL(1,1,j)=', uv_fPL(1,1,j)
        ! write(6,*) 'myid=', myid, ' j=', j, ' uv_fPL(2,1,j)=', uv_fPL(2,1,j)

        wv_fPL(1,1,j) = wv_fPL(1,1,j) + u3PL_itp(ia,ka,j)*u2PL(ip,kp,j) - &
              & iNeg(i)*iNeg(-i)*u3PL_itp(ia+1,ka,j)*u2PL(ip+1,kp,j)
        wv_fPL(2,1,j) = wv_fPL(2,1,j) + u3PL_itp(ia,ka,j)*u2PL(ip+1,kp,j)*iNeg(-i) + &
              & iNeg(i)*u3PL_itp(ia+1,ka,j)*u2PL(ip,kp,j)
      end do 
    end do

    ! if(j==limPL_excw(vgrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1," Finished 0th mode nonlin =====> Linear advection", myid
    end if 
    
    ! linear advection
    buff(:,:) = u1PL_itp(1,1,j)*u2PL(:,:,j) + u2PL(1,1,j)*u1PL_itp(:,:,j)
    vu_fPL(:,:,j) = vu_fPL(:,:,j) + buff(:,:)
    buff(1:2,1) = 0;
    uv_fPL(:,:,j) = uv_fPL(:,:,j) + buff(:,:)


    buff(:,:) = u2PL(1,1,j)*u3PL_itp(:,:,j) + u3PL_itp(1,1,j)*u2PL(:,:,j)
    vw_fPL(:,:,j) = vw_fPL(:,:,j) + buff(:,:)
    buff(1:2,1) = 0;
    wv_fPL(:,:,j) = wv_fPL(:,:,j) + buff(:,:)

    ! if(j==limPL_excw(vgrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1,"=====> Finished Linear Advection - V", myid
    end if 

    ! nonlinear advection: go through a list
    ! fields = {'uu', 'uv', 'uw', 'vu', 'vv', 'vw', 'wu', 'wv', 'ww'}; in the order of 1 to 9, where the first is the passive
    call nonlinInter(nonlin(j, 2), uv_fPL(1,1,j), u2PL(1,1,j), u1PL_itp(1,1,j))
    call nonlinInter(nonlin(j, 4), vu_fPL(1,1,j), u1PL_itp(1,1,j), u2PL(1,1,j))
    call nonlinInter(nonlin(j, 6), vw_fPL(1,1,j), u3PL_itp(1,1,j), u2PL(1,1,j))
    call nonlinInter(nonlin(j, 8), wv_fPL(1,1,j), u2PL(1,1,j), u3PL_itp(1,1,j))

    ! if(j==10) then
    !   write(6,*) "=====> Finished Nonlin Inter v", myid
    ! end if 

    ! if(j==limPL_excw(vgrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1, "=====> Finished Nonlin Inter v", myid
    end if 


    ! if (j == 130) then
    !   write(11) uv_fPL(:,:,j)
    !   write(11) wv_fPL(:,:,j)
    !   ! write(11) vu_fPL(:,:,j)
    !   ! write(11) vw_fPL(:,:,j)
    !   close(11)
    ! end if    


    call der_x(vu_fPL(1,1,j),du2dx,k1F_x)
    call der_z(vw_fPL(1,1,j),du2dz,k1F_z)

    ! if(j==10) then
    !   write(6,*) "=====> Finished Derivatives", myid
    ! end if 

    ! if(j==limPL_excw(vgrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1, "=====> Finished Derivatives -V", myid
    end if 

    do k = 1,Ngal_z
      do i = 1,Ngal_x
        Nu2PL(i,k,j) = du2dx(i,k)+du2dz(i,k)   
      end do
    end do
    
    call four_to_phys_u(u1PL_itp(1,1,j),u2PL(1,1,j),u3PL_itp(1,1,j))

    ! if(j==limPL_excw(vgrid,2,myid)) then
    if(j==10) then
      write(6,*) "t=", MPI_Wtime() - t1, "Finished V =====> Finished Ops in planes", myid
    end if 
  
  end do
  
  ! do myid= 0,np-1
  !   write(6,*) "ops in planes done, myid =", myid
  ! end do 

  ! call MPI_Comm_rank(MPI_COMM_WORLD, rank_true, ierr_local)
  ! call MPI_Comm_size(MPI_COMM_WORLD, nprocs_true, ierr_local)

  ! write(6,*) "MYID CHECK:", "passed myid=", myid, &
  !           "true rank=", rank_true, &
  !           "nprocs=", nprocs_true
  

  deallocate(du1dx,du1dz,du2dx,du2dz,du3dx,du3dz,buff)

end subroutine

subroutine nonlinInter(jlist, x_cPL, uaPL, upPL)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! nonlinear Interactions !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none

  type(nonlinList) jlist
  integer ia, ip, it, ka, kp, kt, la, lp, p, pn
  integer l, iia, kka, iit, kkt, IkkcNeg
  real(8) x_cPL(igal,kgal), uaPL(igal,kgal), upPL(igal,kgal), x, xx

  do l = 1,size(jlist%list,1)
    la = jlist%list(l,1)
    lp = jlist%list(l,3)
    ia = iLkup(la)
    ip = iLkup(lp)
    it = iLkup(la+lp)
    do p = 0,1 
      if (jlist%list(l,2)+jlist%list(l,4) ==0 .and. p ==1) cycle
      pn = p*2-1
      ka = kLkup(pn*iNeg(la)*jlist%list(l,2))
      kp = kLkup(pn*iNeg(lp)*jlist%list(l,4))
      kt = kLkup(pn*(jlist%list(l,2)+jlist%list(l,4)))

      ! write(6,*) "jlist%list(l,2)", jlist%list(l,2), "jlist%list(l,4)",  jlist%list(l,4), kLkup(0)

      
      x_cPL(it,kt) = x_cPL(it,kt) + uaPL(ia,ka)*upPL(ip,kp) - iNeg(la)*iNeg(lp)*uaPL(ia+1,ka)*upPL(ip+1,kp) 
      x_cPL(it+1,kt) = x_cPL(it+1,kt) + iNeg(lp)*uaPL(ia,ka)*upPL(ip+1,kp) + iNeg(la)*uaPL(ia+1,ka)*upPL(ip,kp)
      ! x = uaPL(ia,ka)*upPL(ip,kp) - iNeg(la)*iNeg(lp)*uaPL(ia+1,ka)*upPL(ip+1,kp)
      ! xx = iNeg(lp)*uaPL(ia,ka)*upPL(ip+1,kp) + iNeg(la)*uaPL(ia+1,ka)*upPL(ip,kp)
      ! if (kt == 15 .and. it == 15 .and. abs(x)>1e-12 ) then
      !   write(*,*) x, x_cPL(it,kt), la, lp, jlist%list(l,2), jlist%list(l,4)
      ! end if
    end do
  end do

  ! write(6,*) "kLkup", kLkup(:)

end subroutine

subroutine ops_in_planes2(myid,flagst)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! OPS IN PLANES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! the original ops_in_planes where everything is converted into physical space etc etc

  use declaration
  implicit none
  
  integer i,k,j,myid,flagst,temp
  real(8) ddyi
  real(8), allocatable:: du1dx(:,:),du1dz(:,:),du2dx(:,:),du2dz(:,:),du3dx(:,:),du3dz(:,:)

  allocate(du1dx(igal,kgal),du1dz(igal,kgal))
  allocate(du2dx(igal,kgal),du2dz(igal,kgal))
  allocate(du3dx(igal,kgal),du3dz(igal,kgal))
  
  du1dx = 0d0
  du1dz = 0d0
  du2dx = 0d0
  du2dz = 0d0
  du3dx = 0d0
  du3dz = 0d0


  do j = limPL_excw(ugrid,1,myid),limPL_excw(ugrid,2,myid)

    call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j))

    do k = 1,Ngal_z
      do i = 1,Ngal_x
        uu_cPL(i,k,j) = u1PL    (i,k,j)*u1PL    (i,k,j)
        uw_cPL(i,k,j) = u1PL    (i,k,j)*u3PL    (i,k,j)
        vv_cPL(i,k,j) = u2PL_itp(i,k,j)*u2PL_itp(i,k,j)
        ww_cPL(i,k,j) = u3PL    (i,k,j)*u3PL    (i,k,j)
      end do
    end do
    
    call phys_to_four_du(uu_cPL(1,1,j))
    call phys_to_four_du(uw_cPL(1,1,j))    
    call phys_to_four_du(vv_cPL(1,1,j))    
    call phys_to_four_du(ww_cPL(1,1,j))  

    ! if (myid ==0 .and. j == 5) then  
    !   write(6,*) "uu_cPL(i,k,j)", uu_cPL(:,10,j)
    ! end if 

    if (myid ==0 .and. j == 5) then  
      write(6,*) "uw_cPL(i,k,j)", uw_cPL(:,10,j)
    end if 

    ! if (myid ==0 .and. j == 5) then  
    !   write(6,*) "ww_cPL(i,k,j)", ww_cPL(:,10,j)
    ! end if 

  
    call der_x(uu_cPL(1,1,j),du1dx,k1F_x)
    call der_z(uw_cPL(1,1,j),du1dz,k1F_z)
    call der_x(uw_cPL(1,1,j),du3dx,k1F_x)
    call der_z(ww_cPL(1,1,j),du3dz,k1F_z)

    ! write(6,*) "ugrid j-range:", limPL_excw(ugrid,1,0), limPL_excw(ugrid,2,0)

    ! if (myid ==0 .and. j == 5) then  
    !   write(6,*) "du1dx(i,k,j)", du1dx(:,10)
    ! end if 
  
    do k = 1,Ngal_z
      do i = 1,Ngal_x
        Nu1PL(i,k,j) = du1dx(i,k)+du1dz(i,k)
        Nu3PL(i,k,j) = du3dx(i,k)+du3dz(i,k)
      end do
    end do

    ! if (j == 28) then
    !   write(ext4,'(i5.5)') int(1000d0*(t-700)+kRK)
    !   fnameima = 'output/10_'//ext4//'.dat'
    !   write(*,*) fnameima
    !   open(11,file=fnameima,form='unformatted')
    !   write(11) du1dx
    !   write(11) du1dz
    !   write(11) du3dx
    !   write(11) du3dz
    !   ! write(11) wu_cPL(:,:,j)
    ! end if  

    ! if (j == 28) then
    !   write(ext4,'(i5.5)') int(1000d0*(t-700)+kRK)
    !   fnameima = 'output/20_'//ext4//'.dat'
    !   write(*,*) fnameima
    !   open(11,file=fnameima,form='unformatted')
    !   write(11) Nu1PL(:,:,j)
    !   write(11) Nu3PL(:,:,j)
    ! end if 

  end do
    
  do j = limPL_excw(vgrid,1,myid),limPL_excw(vgrid,2,myid)
    
    call four_to_phys_u(u1PL_itp(1,1,j),u2PL(1,1,j),u3PL_itp(1,1,j))
    
    do k = 1,Ngal_z
      do i = 1,Ngal_x
        uv_fPL(i,k,j) = u1PL_itp(i,k,j)*u2PL(i,k,j)
        wv_fPL(i,k,j) = u3PL_itp(i,k,j)*u2PL(i,k,j)
      end do
    end do
    
    call phys_to_four_du(uv_fPL(1,1,j))
    call phys_to_four_du(wv_fPL(1,1,j))

    ! uv_fPL(:,129,:) = 0
    ! uv_fPL(129:130,:,:) = 0
    ! wv_fPL(:,129,:) = 0
    ! wv_fPL(129:130,:,:) = 0

    ! if (j == 124) then
    !   write(11) uv_fPL(:,:,j)
    !   write(11) wv_fPL(:,:,j)
    !   close(11)
    ! end if    
    
    call der_x(uv_fPL(1,1,j),du2dx,k1F_x)
    call der_z(wv_fPL(1,1,j),du2dz,k1F_z)

    ! if (myid ==0 .and. j == 5) then  
    !   write(6,*) "du2dz(i,k,j)", du2dz(:,10)
    ! end if 


    do k = 1,Ngal_z
      do i = 1,Ngal_x
        Nu2PL(i,k,j) = du2dx(i,k)+du2dz(i,k)   
      end do
    end do
    
  end do
  
  deallocate(du1dx,du1dz,du2dx,du2dz,du3dx,du3dz)

end subroutine

! subroutine ops_in_planes2(myid,flagst)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!! OPS IN PLANES !!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! the original ops_in_planes where everything is converted into physical space etc etc

!   use declaration
!   implicit none
  
!   integer i,k,j,myid,flagst,temp
!   real(8) ddyi
!   real(8), allocatable:: du1dx(:,:),du1dz(:,:),du2dx(:,:),du2dz(:,:),du3dx(:,:),du3dz(:,:),buff(:,:), u1_buff(:,:), u3_buff(:,:)
!   integer la,lp, ia, ip, ka, kp
!   real(8), allocatable :: u1_tmp(:,:), u2_tmp(:,:), u3_tmp(:,:)
!   real(8) :: u1_00_re, u1_00_im, u2_00_re, u2_00_im, u3_00_re, u3_00_im
  

!   allocate(du1dx(igal,kgal),du1dz(igal,kgal))
!   allocate(du2dx(igal,kgal),du2dz(igal,kgal))
!   allocate(du3dx(igal,kgal),du3dz(igal,kgal))
!   allocate(buff(igal,kgal))
!   allocate(u1_buff(igal,kgal))
!   allocate(u3_buff(igal,kgal))
!   allocate(u1_tmp(igal, kgal))
!   allocate(u2_tmp(igal, kgal))
!   allocate(u3_tmp(igal, kgal))

!   u1_tmp = 0d0
!   u2_tmp = 0d0
!   u3_tmp = 0d0
  
!   du1dx = 0d0
!   du1dz = 0d0
!   du2dx = 0d0
!   du2dz = 0d0
!   du3dx = 0d0
!   du3dz = 0d0

!   uu_cPL = 0d0
!   uw_cPL = 0d0
!   vv_cPL = 0d0
!   wu_cPL = 0d0
!   ww_cPL = 0d0
!   uu_cPL_f = 0d0
!   uw_cPL_f = 0d0
!   vv_cPL_f = 0d0
!   ww_cPL_f = 0d0
!   uv_fPL_f = 0d0
!   wv_fPL_f = 0d0

!   !--------------------------------------------------------------------!
!   !!!!!!!!!!!!!!! comment out to go back to normal DNS !!!!!!!!!!!!!!!!
  
!   ! zero mode treatment seperate from nonlinear stuff! 

!   do j = limPL_excw(ugrid,1,myid),limPL_excw(ugrid,2,myid)
!     ! nonlinear interaction into 0th mode (bad, much easier to do if multiply by 4 and do one quadrant but need to think about k = 0 terms carefully)
!     do i = -(Nspec_x/2),Nspec_x/2
!       do k = -(Nspec_z/2),Nspec_z/2-1
!         la = i
!         lp = -i
!         ia = iLkup(la)
!         ip = iLkup(lp)
!         ka = kLkup(iNeg(la)*k)
!         kp = kLkup(iNeg(lp)*(-k))
!         vv_cPL_f(1,1,j) = vv_cPL_f(1,1,j) + u2PL_itp(ia,ka,j)*u2PL_itp(ip,kp,j) - &
!               & iNeg(i)*iNeg(-i)*u2PL_itp(ia+1,ka,j)*u2PL_itp(ip+1,kp,j)
!         vv_cPL_f(2,1,j) = vv_cPL_f(2,1,j) + u2PL_itp(ia,ka,j)*u2PL_itp(ip+1,kp,j)*iNeg(-i) + &
!               & iNeg(i)*u2PL_itp(ia+1,ka,j)*u2PL_itp(ip,kp,j)
!         ! really the imaginary part should be 0 so no need to compute but anyway just in case
!       end do
!     end do
    
!     ! linear advection
!     ! wrong for 0th mode, 0,0 interaction should not be counted twice, but doesn't matter since differentiate = 0
!     ! vv_cPL correct since 0th mode computed earlier
!     ! buff(:,:) = 2*(u1PL(1,1,j)*u1PL(:,:,j))
!     ! uu_cPL_f(:,:,j) = uu_cPL_f(:,:,j) + buff(:,:)

!     buff(:,:) = 2d0 * (u1PL(1,1,j) * u1PL(:,:,j))
!     buff(1:2,1) = 0d0                          ! remove UU so we dont double count- now correct for 0 mode 
!     uu_cPL_f(:,:,j) = uu_cPL_f(:,:,j) + buff(:,:)
!     uu_cPL_f(1,1,j) = uu_cPL_f(1,1,j) + u1PL(1,1,j)*u1PL(1,1,j)  

!     buff(:,:) = 2*(u3PL(1,1,j)*u3PL(:,:,j))
!     ww_cPL_f(:,:,j) = ww_cPL_f(:,:,j) + buff(:,:)

!     buff(:,:) = 2*(u2PL_itp(1,1,j)*u2PL_itp(:,:,j))
!     buff(1:2,1) = 0
!     vv_cPL_f(:,:,j) = vv_cPL_f(:,:,j) + buff(:,:)
    
!     buff(:,:) = u1PL(1,1,j)*u3PL(:,:,j) + u3PL(1,1,j)*u1PL(:,:,j)  
!     ! wu_cPL_f(:,:,j) = wu_cPL_f(:,:,j) + buff(:,:)
!     uw_cPL_f(:,:,j) = uw_cPL_f(:,:,j) + buff(:,:)


!     ! if (myid == 0 .and. j == 5) then
!     !   uw_cPL_f(:,:,j) = 0d0
!     !   buff(:,:)       = 0d0
!     !   u1_buff(:,:)    = 0d0
!     !   u3_buff(:,:)    = 0d0

!     !   ! --- Means (kx=0,kz=0) in packed storage ---
!     !   u1_buff(1,1) = 1d0
!     !   u3_buff(1,1) = 1d0

!     !   ! --- kx=0, kz=3  (k=4) and its negative partner (k = Ngal_z-2) ---
!     !   u1_buff(1,4)         = 2d0
!     !   u1_buff(1,Ngal_z-2)  = 2d0
!     !   u3_buff(1,4)         = 5d0
!     !   u3_buff(1,Ngal_z-2)  = 5d0

!     !   ! --- kx=0, kz=5  (k=6) and its negative partner (k = Ngal_z-4) ---
!     !   u1_buff(1,6)         = 22d0
!     !   u1_buff(1,Ngal_z-4)  = 22d0
!     !   u3_buff(1,6)         = 31d0
!     !   u3_buff(1,Ngal_z-4)  = 31d0

!     !   ! --- kx=1, kz=2  (i=3 is Re(kx=1), k=3 is kz=2) and negative kz partner (k = Ngal_z-1) ---
!     !   u1_buff(3,3)         = 3d0
!     !   u1_buff(3,Ngal_z-1)  = 3d0
!     !   u3_buff(3,3)         = 11d0
!     !   u3_buff(3,Ngal_z-1)  = 11d0

!     !   ! --- linear addback: U*w + W*u with U=u1_buff(1,1), W=u3_buff(1,1) ---
!     !   buff(:,:) = u1_buff(1,1) * u3_buff(:,:) + u3_buff(1,1) * u1_buff(:,:)

!     !   uw_cPL_f(:,:,j) = uw_cPL_f(:,:,j) + buff(:,:)

!     !   ! sanity: print the *linear* Fourier targets
!     !   write(6,*) "DBG uw_f kz=3 kx=0:", uw_cPL_f(1:2,4,j)
!     !   write(6,*) "DBG uw_f kz=5 kx=0:", uw_cPL_f(1:2,6,j)
!     !   write(6,*) "DBG uw_f kz=2 kx=1:", uw_cPL_f(3:4,3,j)
!     ! end if


!   end do 


!   do j = limPL_excw(vgrid,1,myid),limPL_excw(vgrid,2,myid)
!     ! nonlinear interaction into 0th mode
!     do i = -(Nspec_x/2),Nspec_x/2
!       do k = -(Nspec_z/2-1),Nspec_z/2-1
!         la = i
!         lp = -i
!         ia = iLkup(la)
!         ip = iLkup(lp)
!         ka = kLkup(iNeg(la)*k)
!         kp = kLkup(iNeg(lp)*(-k))

!         uv_fPL_f(1,1,j) = uv_fPL_f(1,1,j) + u1PL_itp(ia,ka,j)*u2PL(ip,kp,j) - &
!               & iNeg(i)*iNeg(-i)*u1PL_itp(ia+1,ka,j)*u2PL(ip+1,kp,j)
!         uv_fPL_f(2,1,j) = uv_fPL_f(2,1,j) + u1PL_itp(ia,ka,j)*u2PL(ip+1,kp,j)*iNeg(-i) + &
!               & iNeg(i)*u1PL_itp(ia+1,ka,j)*u2PL(ip,kp,j)
!         wv_fPL_f(1,1,j) = wv_fPL_f(1,1,j) + u3PL_itp(ia,ka,j)*u2PL(ip,kp,j) - &
!               & iNeg(i)*iNeg(-i)*u3PL_itp(ia+1,ka,j)*u2PL(ip+1,kp,j)
!         wv_fPL_f(2,1,j) = wv_fPL_f(2,1,j) + u3PL_itp(ia,ka,j)*u2PL(ip+1,kp,j)*iNeg(-i) + &
!               & iNeg(i)*u3PL_itp(ia+1,ka,j)*u2PL(ip,kp,j)
!       end do 
!     end do

!     ! linear advection
!     buff(:,:) = u1PL_itp(1,1,j)*u2PL(:,:,j) + u2PL(1,1,j)*u1PL_itp(:,:,j)
!     ! vu_fPL_f(:,:,j) = vu_fPL_f(:,:,j) + buff(:,:)
!     buff(1:2,1) = 0;
!     uv_fPL_f(:,:,j) = uv_fPL_f(:,:,j) + buff(:,:)


!     buff(:,:) = u2PL(1,1,j)*u3PL_itp(:,:,j) + u3PL_itp(1,1,j)*u2PL(:,:,j)
!     ! vw_fPL_f(:,:,j) = vw_fPL_f(:,:,j) + buff(:,:)
!     buff(1:2,1) = 0;
!     wv_fPL_f(:,:,j) = wv_fPL_f(:,:,j) + buff(:,:)

    
!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "u2PL(i,k,j)", u2PL(:,1,j)
!     ! end if 
  
!   end do 

!   !------------------------------------------------------------------------!

  
!   do j = limPL_excw(ugrid,1,myid),limPL_excw(ugrid,2,myid)

!     ! !--------------------------------------------------------------------!
!     ! !!!!!!!!!!!!!! comment out to go back to normal DNS !!!!!!!!!!!!!!!!
!     ! ! apply masks before the transform to physical space
!     ! do k = 1,Ngal_z
!     !   do i = 1,Ngal_x+2
!     !     u1PL(i,k,j) = u1PL(i,k,j)*mask_U(i,k,j)
!     !     u3PL(i,k,j) = u3PL(i,k,j)*mask_W(i,k,j)
!     !     u2PL_itp(i,k,j) = u2PL_itp(i,k,j)*mask_V_itp(i,k,j)
!     !   end do
!     ! end do
!     ! !--------------------------------------------------------------------!

!     ! u1PL(1:2,1,j) = 0d0
!     ! u2PL_itp(1:2,1,j) = 0d0
!     ! u3PL(1:2,1,j) = 0d0

!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "u3PL(i,k,j)", u3PL(:,10,j)
!     ! end if 
     

!     ! call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j))

!     ! if (myid==0 .and. j==5) then
!     !   do k = 1,Ngal_z
!     !     do i = 1,Ngal_x
!     !       u1PL(i,k,j) = 2d0*cos(2d0*pi*3d0*(k-1)/Ngal_z) &
!     !         + 3d0*cos(2d0*pi*1d0*(i-1)/Ngal_x)*cos(2d0*pi*2d0*(k-1)/Ngal_z) &
!     !         + 22d0*cos(2d0*pi*5d0*(k-1)/Ngal_z)

!     !       u3PL(i,k,j) = 5d0*cos(2d0*pi*3d0*(k-1)/Ngal_z) &
!     !         + 11d0*cos(2d0*pi*1d0*(i-1)/Ngal_x)*cos(2d0*pi*2d0*(k-1)/Ngal_z) &
!     !         + 31d0*cos(2d0*pi*5d0*(k-1)/Ngal_z)
!     !     end do
!     !   end do
!     ! endif

!     ! if (myid==0 .and. j==5) then
!     !   do k = 1, Ngal_z
!     !     do i = 1, Ngal_x
!     !       u1PL(i,k,j) = 4d0  * cos(2d0*pi*3d0*(k-1)/Ngal_z) &
!     !         + 44d0 * cos(2d0*pi*5d0*(k-1)/Ngal_z)

!     !       u3PL(i,k,j) = 10d0 * cos(2d0*pi*3d0*(k-1)/Ngal_z) &
!     !         + 62d0 * cos(2d0*pi*5d0*(k-1)/Ngal_z)
!     !     end do
!     !   end do
!     ! end if

!     ! do k = 1,Ngal_z
!     !   do i = 1,Ngal_x
!     !     uu_cPL(i,k,j) = u1PL    (i,k,j)*u1PL    (i,k,j)
!     !     uw_cPL(i,k,j) = u1PL    (i,k,j)*u3PL    (i,k,j)
!     !     vv_cPL(i,k,j) = u2PL_itp(i,k,j)*u2PL_itp(i,k,j)
!     !     ww_cPL(i,k,j) = u3PL    (i,k,j)*u3PL    (i,k,j)
        
!     !   end do
!     ! end do 


!     ! Copy Fourier-plane data into temps (including padding cols 1..Ngal_x+2)
!     u1_tmp(:,:) = u1PL(1:Ngal_x,1:Ngal_z,j)        ! will be overwritten by phys data anyway
!     u2_tmp(:,:) = u2PL_itp(1:Ngal_x,1:Ngal_z,j)
!     u3_tmp(:,:) = u3PL(1:Ngal_x,1:Ngal_z,j)

!     u1_00_re = u1PL(1,1,j); u1_00_im = u1PL(2,1,j)
!     u2_00_re = u2PL_itp(1,1,j); u2_00_im = u2PL_itp(2,1,j)
!     u3_00_re = u3PL(1,1,j); u3_00_im = u3PL(2,1,j)

!     ! Zero mean mode ONLY in the working arrays (will restore right after)
!     u1PL(1:2,1,j)      = 0d0
!     u2PL_itp(1:2,1,j)  = 0d0
!     u3PL(1:2,1,j)      = 0d0

!     call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j))

!     ! Restore the original Fourier mean slots so later solver steps see unchanged spectral fields
!     u1PL(1,1,j) = u1_00_re;  u1PL(2,1,j) = u1_00_im
!     u2PL_itp(1,1,j) = u2_00_re;  u2PL_itp(2,1,j) = u2_00_im
!     u3PL(1,1,j) = u3_00_re;  u3PL(2,1,j) = u3_00_im

!     ! Copy the physical fields out to temp arrays (so uu/uw/vv/ww don’t depend on u1PL being reused later)
!     do k = 1, Ngal_z
!       do i = 1, Ngal_x
!         u1_tmp(i,k) = u1PL(i,k,j)
!         u2_tmp(i,k) = u2PL_itp(i,k,j)
!         u3_tmp(i,k) = u3PL(i,k,j)
!       end do
!     end do

!     ! Now build the products from the temp physical fields
!     do k = 1, Ngal_z
!       do i = 1, Ngal_x
!         uu_cPL(i,k,j) = u1_tmp(i,k) * u1_tmp(i,k)
!         uw_cPL(i,k,j) = u1_tmp(i,k) * u3_tmp(i,k)
!         vv_cPL(i,k,j) = u2_tmp(i,k) * u2_tmp(i,k)
!         ww_cPL(i,k,j) = u3_tmp(i,k) * u3_tmp(i,k)
!       end do
!     end do
    
!     call phys_to_four_du(uu_cPL(1,1,j))
!     call phys_to_four_du(uw_cPL(1,1,j))    
!     call phys_to_four_du(vv_cPL(1,1,j))    
!     call phys_to_four_du(ww_cPL(1,1,j))  

!     !--------------------------------------------------------------------!
!     !!!!!!!!!!!!!! comment out to go back to normal DNS !!!!!!!!!!!!!!!!
!     ! add zero mode +linear stuff back in before calculating derivatives (_f is 0 mode stuff and Linear advec)

!     ! wu_cPl(:,:,j) =  uw_cPL(:,:,j)

!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "uw_cPL(i,k,j)", uw_cPL(:,10,j)
!     ! end if 

!     do k = 1,Ngal_z
!       do i = 1,Ngal_x +2
!         uu_cPL(i,k,j) = uu_cPL(i,k,j) + uu_cPL_f(i,k,j)
!         uw_cPL(i,k,j) = uw_cPL(i,k,j) + uw_cPL_f(i,k,j)
!         vv_cPL(i,k,j) = vv_cPL(i,k,j) + vv_cPL_f(i,k,j)
!         ww_cPL(i,k,j) = ww_cPL(i,k,j) + ww_cPL_f(i,k,j)
!         ! wu_cPl(i,k,j) = wu_cPl(i,k,j) + wu_cPl_f(i,k,j)
!       end do
!     end do
!     !--------------------------------------------------------------------!

!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "uu_cPL(i,k,j)", uu_cPL(:,10,j)
!     ! end if 

!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "ww_cPL(i,k,j)", ww_cPL(:,10,j)
!     ! end if 
    
!     if (myid ==0 .and. j == 5) then  
!       ! write(6,*) "uw_cPL_f(i,k,j)", uw_cPL_f(:,4,j)
!       write(6,*) "uw_cPL(i,k,j)", uw_cPL(:,10,j)
!       write(6,*) 
!     end if 


!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "uw_cPL(i,k,j)", uw_cPL(:,9,j)
!     !   write(6,*) "uw_cPL(i,k,j)", uw_cPL(:,Ngal_z-7,j)
!     ! end if 

!     call der_x(uu_cPL(1,1,j),du1dx,k1F_x)
!     call der_z(uw_cPL(1,1,j),du1dz,k1F_z)
!     call der_x(uw_cPL(1,1,j),du3dx,k1F_x)
!     call der_z(ww_cPL(1,1,j),du3dz,k1F_z)

!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "du3dx(i,k,j)", du3dx(:,10)
!     ! end if 

      
!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "du3dz(i,k,j)", du3dz(:,10)
!     ! end if 
  
!     do k = 1,Ngal_z
!       do i = 1,Ngal_x
!         Nu1PL(i,k,j) = du1dx(i,k)+du1dz(i,k)
!         Nu3PL(i,k,j) = du3dx(i,k)+du3dz(i,k)
!       end do
!     end do
!   end do
    
!   do j = limPL_excw(vgrid,1,myid),limPL_excw(vgrid,2,myid)

!     ! ! apply mask to v grid things 
!     ! do k = 1,Ngal_z
!     !   do i = 1,Ngal_x+2
!     !     u1PL(i,k,j) = u1PL_itp(i,k,j)*mask_U_itp(i,k,j)
!     !     u3PL(i,k,j) = u3PL_itp(i,k,j)*mask_W_itp(i,k,j)
!     !     u2PL_itp(i,k,j) = u2PL(i,k,j)*mask_V(i,k,j)
!     !   end do
!     ! end do

!     u1PL_itp(1:2,1,j) =  0
!     u2PL(1:2,1,j) =  0
!     u3PL_itp(1:2,1,j) =  0
    
!     call four_to_phys_u(u1PL_itp(1,1,j),u2PL(1,1,j),u3PL_itp(1,1,j))
    
!     do k = 1,Ngal_z
!       do i = 1,Ngal_x
!         uv_fPL(i,k,j) = u1PL_itp(i,k,j)*u2PL(i,k,j)
!         wv_fPL(i,k,j) = u3PL_itp(i,k,j)*u2PL(i,k,j)
!       end do
!     end do
    
!     call phys_to_four_du(uv_fPL(1,1,j))
!     call phys_to_four_du(wv_fPL(1,1,j))

!     !--------------------------------------------------------------------!
!     !!!!!!!!!!!!!! comment out to go back to normal DNS !!!!!!!!!!!!!!!!

!     do k = 1,Ngal_z
!       do i = 1,Ngal_x
!         uv_fPL(i,k,j) = uv_fPL(i,k,j) + uv_fPL_f(i,k,j)
!         ! vu_fPL(i,k,j) = vu_fPL(i,k,j) + vu_fPL_f(i,k,j)
!         wv_fPL(i,k,j) = wv_fPL(i,k,j) + wv_fPL_f(i,k,j)
!         ! vw_fPL(i,k,j) = vw_fPL(i,k,j) + vw_fPL_f(i,k,j)
!       end do
!     end do
!     !--------------------------------------------------------------------!

    
!     call der_x(uv_fPL(1,1,j),du2dx,k1F_x)
!     call der_z(wv_fPL(1,1,j),du2dz,k1F_z)

!     ! if (myid ==0 .and. j == 5) then  
!     !   write(6,*) "du2dz(i,k,j)", du2dz(:,10)
!     ! end if 

!     do k = 1,Ngal_z
!       do i = 1,Ngal_x
!         Nu2PL(i,k,j) = du2dx(i,k)+du2dz(i,k)   
!       end do
!     end do
    
!   end do
  
!   deallocate(du1dx,du1dz,du2dx,du2dz,du3dx,du3dz)
!   deallocate(buff)

! end subroutine


subroutine phys_to_four_du(duPL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Physical Space to Fourier Space u 
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none

  real(8) duPL(Ngal_x+2,Ngal_z)

  call rft(duPL,Ngal_x+2,Ngal_z,-1,buffRal_x)
  call cft(duPL,Ngal_x+2,2,(Nspec_x+2)/2,-1,buffCal_z)
  
  duPL(:,Ngal_z/2+1)=0d0 !oddball advective term = 0

end subroutine

subroutine phys_to_four_N(duPL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Physical Space to Fourier Space u 
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none


  real(8) duPL(Nspec_x+2,Nspec_z)

  call rft(duPL,Nspec_x+2,Nspec_z,-1,buffR_x)
  call cft(duPL,Nspec_x+2,2,(Nspec_x+2)/2,-1,buffC_z)
  
  !duPL(:,N(2,iband)/2+1)=0d0

end subroutine
subroutine planes_to_modes_UVP (x,xPL,grid,nygrid,nygrid_LB,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES  NEW !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  real :: t_start, t_end

  integer :: i,k,j,jminS,jmaxS,jminR,jmaxR,grid
  integer :: column, nprocs
  integer inode,yourid
  integer msizeR,msizeS
  complex(8)   x( jlim(1,grid) : jlim(2,grid), columns_num(myid) )
  real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable:: buffS(:,:),buffR(:,:)
  integer, intent(in) :: nygrid,nygrid_LB

  ! write(6,*) "starting self transpose"

  ! Loop for itself
  ! Transpose the cube that it already owns


  ! call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  ! write(6,*) "DEBUG myid=", myid, "MPI nprocs=", nprocs, "np var=", np
  ! call flush(6)



  jminR = max(limPL_excw(grid,1,myid),jlim(1,grid)+1)  ! Select the planes to transpose 
  jmaxR = min(limPL_excw(grid,2,myid),jlim(2,grid)-1)

  ! write(*,*) "rank", myid, "x bounds:", lbound(x,1), ubound(x,1), "cols:", lbound(x,2), ubound(x,2)
  

  
  if (jminR==nygrid_LB+1 .and. jmaxR>=jminR) then   ! Special cases: walls
    jminR = jminR-1
  end if
  if (jmaxR==nygrid   .and. jmaxR>=jminR) then
    jmaxR = jmaxR+1
  end if
  ! write(*,*) "rank", myid, "jminR/jmaxR:", jminR, jmaxR, "grid", grid

  do j = jminR,jmaxR
    do column = 1,columns_num(myid)
      i = columns_i(column,myid)
      k = columns_k(column,myid) - dk(column,myid)
      x(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
    end do
  end do

  !end do

  ! write(6,*) "finished self transpose"

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
        jminS = max(limPL_excw(grid,1,  myid),jlim(1,grid)+1)  ! Select the planes to be SENT.
        jmaxS = min(limPL_excw(grid,2,  myid),jlim(2,grid)-1)  ! max and min because maybe this proc needs less planes that the other proc has
        jminR = max(limPL_excw(grid,1,yourid),jlim(1,grid)+1)  ! Select the planes to be RECEIVED
        jmaxR = min(limPL_excw(grid,2,yourid),jlim(2,grid)-1)

        ! Adding the walls =

        if (jminS==nygrid_LB+1  ) then
          jminS=jminS-1
        end if
        if (jmaxS==nygrid) then
          jmaxS=jmaxS+1
        end if
        if (jminR==nygrid_LB+1  ) then
          jminR=jminR-1
        end if
        if (jmaxR==nygrid) then
          jmaxR=jmaxR+1
        end if
        allocate(buffS(jminS:jmaxS,columns_num(yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(  myid)))

        ! if (myid ==0 .and. yourid == 4 .and. jband == 2) then
        !   write(6,*) "buffs", size(buffS,1), size(buffS,2)
        ! end if 

        msizeS = 2*(columns_num(yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(yourid)
            i = columns_i(column,yourid)
            k = columns_k(column,yourid) - dk(column,yourid)
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do

        ! call cpu_time(t_start)
        
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
&                         MPI_COMM_WORLD,status,ierr)

        ! call cpu_time(t_end)

        ! if (myid ==0 .and. yourid == 4 .and. jband == 2) then
        !   write(6,*) "cpu time", t_start, t_end
        ! end if

        do j=jminR,jmaxR
          do column = 1,columns_num(myid)
            x(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      ! end do
    end if
  end do

end subroutine


! subroutine planes_to_modes_phys_lims (x,xPL,nystart,nyend,grid,myid,status,ierr)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES. USED IN FFT TRID LU !!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! i dont think this is called at all...

! ! Prepare the vectors for the Fourier transform
! ! The procs broadcast the data they have of a plane, and receive the data of a pencil,
! !  they also transpose the data for the Fourier transform XZY -> YXZ

!   use declaration
!   implicit none

!   include 'mpif.h'             ! MPI variables
!   integer status(MPI_STATUS_SIZE),ierr,myid

!   integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid
!   integer column,nystart,nyend
!   integer jband
!   integer inode,yourid
!   integer msizeR,msizeS
!   type(cfield) x
!   real(8)      xPL(Nspec_x+2,Nspec_z,limPL_FFT(grid,1,myid):limPL_FFT(grid,2,myid))
!   complex(8), allocatable:: buffS(:,:),buffR(:,:)

!   ! Loop for itself
!   ! Transpose the cube that it already owns
!   plband = bandPL_FFT(myid) ! Return the band (phys) the proc works at
!   ! do iband = sband,eband
!   jminR = max(max(limPL_FFT(grid,1,myid),jlim(1,grid)),nystart)  ! Select the planes to transpose 
!   jmaxR = min(min(limPL_FFT(grid,2,myid),jlim(2,grid)),nyend  )
  
!   do j = jminR,jmaxR
!     do column = 1,columns_num(myid)
!       i = columns_i(column,myid)
!       k = columns_k(column,myid) - dk_phys(column,myid)
!       x%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
!     end do
!   end do
!   ! end do

!   do inode = 1,pnodes-1
!     yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
!     if (yourid<np) then
!       ! do iband = sband,eband
!       !jband = crossband(iband,yourid)
!       ! jband = iband
!       jminS = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid)),nystart)  ! Select the planes to be SENT.
!       jmaxS = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid)),nyend  )  ! max and min because maybe this proc needs less planes that the other proc has
!       jminR = max(max(limPL_FFT(grid,1,yourid),jlim(1,grid)),nystart)  ! Select the planes to be RECEIVED
!       jmaxR = min(min(limPL_FFT(grid,2,yourid),jlim(2,grid)),nyend  )

!       allocate(buffS(jminS:jmaxS,columns_num(yourid)))
!       allocate(buffR(jminR:jmaxR,columns_num(  myid)))
!       msizeS = 2*(columns_num(yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
!       msizeR = 2*(columns_num(  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
!       msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
!       msizeR = max(msizeR,0)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       do j=jminS,jmaxS
!         do column = 1,columns_num(yourid)
!           i = columns_i(column,yourid)
!           k = columns_k(column,yourid) - dk_phys(column,yourid)
!           buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
!         end do
!       end do
!       call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*nband+11*nband, &   ! SEND_RECV so it can send and receive at the same time
! &                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*nband+7*nband, &
! &                         MPI_COMM_WORLD,status,ierr)
!       do j=jminR,jmaxR
!         do column = 1,columns_num(myid)
!           x%f(j,column) = buffR(j,column)                         ! Store the data received
!         end do
!       end do
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       deallocate(buffR,buffS)
!       ! end do
!     end if
!   end do
! end subroutine


subroutine planes_to_modes_NUVP(x,xPL,grid,nygrid,nygrid_LB,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ


  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,grid
  integer column
  integer inode,yourid
  integer msizeR,msizeS
  complex(8)  x(jlim(1,grid)+1:jlim(2,grid)-1,columns_num(myid))
  real(8)     xPL(igal,kgal,jgal(grid,1):jgal(grid,2))
  complex(8), allocatable:: buffS(:,:),buffR(:,:)
  integer     nygrid, nygrid_LB

  ! Loop for itself
  ! Transpose the cube that it already owns


  jminR = max(limPL_excw(grid,1,myid),jlim(1,grid)+1)  ! Select the planes to transpose 
  jmaxR = min(limPL_excw(grid,2,myid),jlim(2,grid)-1)

  ! write(6,*) "jminR", jminR, "jmaxR", jmaxR, "myid", myid

  do j = jminR,jmaxR
    do column = 1,columns_num(myid)
      i = columns_i(column,myid)
      k = columns_k(column,myid) - dk(column,myid)
      ! write(6,*) "j=", j, "mpi", myid
      ! write(6,*) 'myid=', myid, ' j=', j, ' l3=', lbound(XPL,3), ' u3=', ubound(XPL,3)
      !call flush(6)

      x(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      ! write(6,*) "j=", j, "mpi", myid
    end do
  end do


  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      jminS = max(limPL_excw(grid,1,  myid),jlim(1,grid)+1)  ! Select the planes to be SENT.
      jmaxS = min(limPL_excw(grid,2,  myid),jlim(2,grid)-1)  ! max and min because maybe this proc needs less planes that the other proc has
      jminR = max(limPL_excw(grid,1,yourid),jlim(1,grid)+1)  ! Select the planes to be RECEIVED
      jmaxR = min(limPL_excw(grid,2,yourid),jlim(2,grid)-1)
      allocate(buffS(jminS:jmaxS,columns_num(yourid)))
      allocate(buffR(jminR:jmaxR,columns_num(  myid)))
      msizeS = 2*(columns_num(yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
      msizeR = 2*(columns_num(  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
      msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
      msizeR = max(msizeR,0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=jminS,jmaxS
        do column = 1,columns_num(yourid)
          i = columns_i(column,yourid)
          k = columns_k(column,yourid) - dk(column,yourid)
          buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
        end do
      end do
      call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
&                         MPI_COMM_WORLD,status,ierr)


!       call MPI_SENDRECV( buffS, msizeS, MPI_DOUBLE_COMPLEX, yourid, 77*yourid + 53*myid, &
! &                        buffR, msizeR, MPI_DOUBLE_COMPLEX, yourid, 53*yourid + 77*myid, &
! &                        MPI_COMM_WORLD, status, ierr )


      do j=jminR,jmaxR
        do column = 1,columns_num(myid)
          x(j,column) = buffR(j,column)                         ! Store the data received
        end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      deallocate(buffR,buffS)
      !end do
    end if
  end do

end subroutine

subroutine record_out(u1,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   RECORD OUT   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  ! use littleharsh_mod
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  complex(8)   u1( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )
  integer nx,nz, i 
  integer j,jmax,iproc
  real(8), allocatable:: buffSR(:,:)
  integer, allocatable:: dummint(:)
  real(8) Uslip  

  ! Rebuilding N 
  if (.not. allocated(N)) allocate(N(4,0:4))


  N= 0 
  N(1,1:3) = Nspec_x
  N(1,4) = -2
  N(2,1:3) = Nspec_z
  N(3,3) = nyv
  N(4,3) = nyu

  if (myid == 0) then
      write(6,*) "N:"
      do i = 1,4
      write(6,*) N(i,0:4)
      end do
  end if 

  if (myid/=0) then
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,123*myid,MPI_COMM_WORLD,ierr)
      end do
      do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,124*myid,MPI_COMM_WORLD,ierr)
      end do
      do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,125*myid,MPI_COMM_WORLD,ierr)
      end do
      do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
      call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,126*myid,MPI_COMM_WORLD,ierr)
      end do

      deallocate(buffSR)
  else
      write(ext4,'(i5.5)') int(100d0*(t))!int(t)!
      allocate(dummint(88))
      dummint = 0
      !!!!!!!!!!!!!    u1    !!!!!!!!!!!!!
      fnameima = 'output/u1_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
      open(10,file=fnameima,form='unformatted')
      write(10) t,Re,alp,bet,mpgx,nband,iter,dummint 
      write(10) N
      write(10) yu,ddthetavi,dthdyu
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,1,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
      do iproc = 1,np-1
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
          call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,123*iproc,MPI_COMM_WORLD,status,ierr)
          write(10) j,1,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
      end do
      close(10)
      !!!!!!!!!!!!!    u2    !!!!!!!!!!!!!
      fnameima='output/u2_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
      open(10,file=fnameima,form='unformatted')
      write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
      write(10) N
      write(10) yv,ddthetavi,dthdyv
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,2,nx,nz,yv(j),buffSR
      end do
      deallocate(buffSR)
      do iproc = 1,np-1
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
          call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,124*iproc,MPI_COMM_WORLD,status,ierr)
          write(10) j,2,nx,nz,yv(j),buffSR
      end do
      deallocate(buffSR)
      end do
      close(10)
      !!!!!!!!!!!!!    u3    !!!!!!!!!!!!!
      fnameima = 'output/u3_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
      open(10,file=fnameima,form='unformatted')
      write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
      write(10) N
      write(10) yu,ddthetavi,dthdyu
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,3,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
      do iproc = 1,np-1
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
          call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,status,ierr)
          write(10) j,3,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
      end do
      close(10)
      !!!!!!!!!!!!!    p     !!!!!!!!!!!!!
      fnameima = 'output/p_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
      open(10,file=fnameima,form='unformatted')
      write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
      write(10) N
      write(10) yu,ddthetavi,dthdyu
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
      call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
      write(10) j,4,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
      do iproc = 1,np-1
      nx = Nspec_x+2
      nz = Nspec_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(pgrid,1,iproc),limPL_incw(pgrid,2,iproc)
          call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,126*iproc,MPI_COMM_WORLD,status,ierr)
          write(10) j,4,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
      end do
      close(10)
      deallocate(dummint)
  end if

  if (myid==0) then

      Uslip = ((u1(1,1))*(-1d0-yu(0))+(u1(0,1))*(yu(1)+1d0))/(yu(1)-yu(0))

      ! call flowrateIm(Qx,u1(nyu_LB,1))
      call flowrateIm(Qx,u1(:,1))
      ! call maxvel(u1(nyu_LB,1))
      call maxvel(u1(:,1))
      write(*,*) ''
      write(*,*) 'iter',iter
      write(*,*) 't   ',t
      write(*,*) 'dtv ',dtv
      write(*,*) 'dtc ',dtc
      write(*,*) 'dt  ',dt
      write(*,*) 'err ',err
      write(*,*) 'Qx  ',Qx
      if (flag_ctpress==0) then
      write(*,*) 'QxT ',QxT
      write(*,*) 'mpgx',mpgx
      write(*,*) 'dpgx',dgx
      else
      write(*,*) 'mpgx',mpgx
      end if
      write(*,*) 'Umax',Umax
      write(*,*) 'Uslp',Uslip
  !    write(*,*) 'utau',utau
  ! Save to history file
      if (flag_ctpress==0) then
      write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,QxT,mpgx,dgx,Umax,Uslip
      else
      write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,    mpgx,    Umax,Uslip
      end if
      flush(30)
  end if

end subroutine



subroutine write_Qcrit(myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   RECORD OUT   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer nx,nz
  integer j,jmax,iproc
  real(8), allocatable:: buffSR(:,:)
  integer, allocatable:: dummint(:)


  if (myid/=0) then
    nx = Ngal_x+2
    nz = Ngal_z
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call MPI_SEND(Qcrit(:,:,j),nx*nz,MPI_REAL8,0,127*myid,MPI_COMM_WORLD,ierr)
    end do
  else
    write(ext4,'(i5.5)') int(10000d0*t)!int(t)!
    allocate(dummint(88))
    dummint = 0
    !!!!!!!!!!!!!    Qcrit    !!!!!!!!!!!!!
    fnameima='output/Qcrit_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    open(10,file=fnameima,form='unformatted')
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    write(10) N
    write(10) yv,ddthetavi,dthdyv
    nx = Ngal_x+2
    nz = Ngal_z
    !allocate(buffSR(nx,nz))
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      !call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,2,nx,nz,yv(j),Qcrit(:,:,j)
    end do
    !deallocate(buffSR)
    do iproc = 1,np-1
      nx = Ngal_x+2
      nz = Ngal_z
      allocate(buffSR(nx,nz))
      do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,127*iproc,MPI_COMM_WORLD,status,ierr)
        write(10) j,2,nx,nz,yv(j),buffSR
      end do
      deallocate(buffSR)
    end do
    close(10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(dummint)
  end if

end subroutine

subroutine u_to_buff(buffSR,u,nx,nz,igal,kgal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    U to BUFF   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rearrange z-modes before writing and after reading

    implicit none

    integer nx,nz,igal,kgal
    real(8) u(igal,kgal)
    real(8) buffSR(nx,nz)
    integer i,k,dkk

    do k = 1,nz/2
    do i = 1,nx
        buffSR(i,k) = u(i,k    )
    end do
    end do
    do i = 1,nx
    !buffSR(i,nz/2+1) = 0d0
    end do
    do k = nz/2+1,nz
    dkk = kgal-nz
    do i = 1,nx
        buffSR(i,k) = u(i,k+dkk)
    end do
    end do

end subroutine

subroutine buff_to_u(u,buffSR,nx,nz,igal,kgal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    BUFF to U   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rearrange z-modes before writing and after reading

  implicit none
  integer nx,nz,igal,kgal
  real(8) u(igal,kgal)
  real(8) buffSR(nx,nz)
  integer i,k,dkk

  do k = 1,min(nz/2,kgal/2)
    do i = 1,min(nx,igal)
      u(i,k    ) = buffSR(i,k)
    end do
  end do
  dkk = kgal-nz
  do k = nz-min(nz/2,kgal/2)+1,nz
    do i = 1,min(nx,igal)
      u(i,k+dkk) = buffSR(i,k)
    end do
  end do

end subroutine
