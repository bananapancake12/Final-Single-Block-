! *************************************************************************************************************** !
!
! Last modified by C 12/05/2016
!   - Set up for SHS
!
!  TODO
!    - Pressure correction of ghost points
!    - flowrate_corr LUsol0 boundary conditions if not no-slip?
!
! *************************************************************************************************************** !
!
! Contains:
!       littleharsh - Main program : Calls start, getini, RHS0_u1, RHS0_u2, RHS0_u3, nonlinear, solveU, solveV, meanflow_ctU, SolveP, v_corr, meanflow_ctP, finalize
!                     - Does the main loop
!        
!
!  Everything else: Moved to littleharsh_mod!! 

!            littleharsh      - 
!            divergence       -
!            laplacian_U      -
!            laplacian_V      -
!            error            -
!            finalize         -
!            flowrateIm       -
!            flowrateRe       -
!            flowrate_corr    -
!            maxvel           -
!            meanflow_ctP     -
!            meanflow_ctU     -
!            solveP           -
!            RHS0_u1          -
!            RHS0_u2          -
!            RHS0_u3          -
!            solveU           -
!            solveV           -
!            v_corr           -
!
! *************************************************************************************************************** !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! LITTLE HARSH NEW  !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program littleharsh

  use declaration
  ! use init_mod
  ! use littleharsh_mod
  ! use tridLU_3D
  implicit none

  include 'mpif.h'                         
  integer status(MPI_STATUS_SIZE),ierr,myid 

  complex(8), allocatable :: u1(:,:), u2(:,:), u3(:,:), p(:,:), psi(:,:)
  complex(8), allocatable :: div(:,:)
  complex(8), allocatable :: Nu1(:,:), Nu2(:,:), Nu3(:,:)
  complex(8), allocatable :: du1(:,:), du2(:,:), du3(:,:)

  real(8), allocatable :: a_ugrid(:,:)
  real(8), allocatable :: a_vgrid(:,:)
  
  real(8):: z,sigmaz
  integer column,i,k,j

  real(8):: flagslstat

  


  call MPI_INIT(ierr)
  if (ierr.ne.MPI_SUCCESS) pause 'ierr'
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)      ! np   = Number of processes in the group of comm
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)    ! myid = rank of the calling process

  if (myid==0) then
    write(*,*) ''
    write(*,*) '---------------LITTLE HARSH---------------'
    write(*,*) "         'That's a little harsh'          "
    write(*,*) "                   meow                   "
    write(*,*) '------------------------------------------'
    write(*,*) ''
  end if

  t1 = MPI_Wtime()

  ! Initialise
  call start(myid,status,ierr)

  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Finished Start"
  end if

! itersl=iter0
! nstatsl = nwrite

  allocate( u1( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) ) )
  allocate( u2( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) ) )
  allocate( u3( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) ) )

  allocate( du1( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) ) )
  allocate( du2( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) ) )
  allocate( du3( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) ) )

  allocate( Nu1( jlim(1,ugrid) + 1  : jlim(2,ugrid) - 1,  columns_num(myid) ) )
  allocate( Nu2( jlim(1,vgrid) + 1  : jlim(2,vgrid) - 1,  columns_num(myid) ) )
  allocate( Nu3( jlim(1,ugrid) + 1  : jlim(2,ugrid) - 1,  columns_num(myid) ) )

  allocate( p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) ) )
  allocate( psi( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) ) )
  allocate( div( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) ) )

  ! allocate( a  (2))

  allocate( a_ugrid(3,jlim(1,ugrid):jlim(2,ugrid)) )
  allocate( a_vgrid(3,jlim(1,vgrid):jlim(2,vgrid)) )

  ! write(6,*) "u1", jlim(1,ugrid)    , jlim(2,ugrid)

  u1  = 0d0
  u2  = 0d0
  u3  = 0d0
  du1 = 0d0
  du2 = 0d0
  du3 = 0d0
  Nu1 = 0d0
  Nu2 = 0d0
  Nu3 = 0d0
  p   = 0d0
  psi = 0d0
  div = 0d0
  a_ugrid = 0d0
  a_vgrid = 0d0



  if(myid==0) then
    write(6,*) "t=", MPI_Wtime() - t1,"=====> Calling getini "
  end if


  ! Get the initial conditions
  call getini(u1,u2,u3,p,div,myid,status,ierr)
  write(6,*) 'finished getini', myid

!   if(myid==0) then
!     write(6,*) "u1", u1(jlim(1,ugrid),:)
!   end if 

nextqt = floor(t*10d0)/10d0+0.1d0
!write(*,*) 'nextqt'
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid==0) then
    write(*,*) ''
    write(*,*) 'Initialised'
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! MAIN LOOP 
  do while (t<maxt) ! This is the original condition
  ! do while (t<maxt .AND. iter <4)
    ! Runge-Kutta substeps
    do kRK = 1,3
    

      ! Build     linear terms of right-hand-side of Navier-Stokes equation
      call RHS0_u1(du1,u1,Nu1,p,myid)
      call RHS0_u2(du2,u2,Nu2,p,myid)
      call RHS0_u3(du3,u3,Nu3,p,myid)

        ! if(myid==0) then
        !     write(6,*) "du1", du1(0,:)
        ! end if 

      if(myid==0) then
        write(6,*) "t=", MPI_Wtime() - t1,"finished RHS =====> Building Nonlinear"
      end if 

      ! Build non-linear terms of right-hand-side of Navier-Stokes equation
      call nonlinear(Nu1,Nu2,Nu3,u1,u2,u3,du1,du2,du3,p,div,myid,status,ierr)

    !   if(myid==0) then
    !     write(6,*) "du1 after nonlin", du1(0,:)
    !   end if 
      
      if(myid==0) then
        write(6,*) "t=", MPI_Wtime() - t1,"Finished Nonlinear =====> Solving"
      end if 

      ! Resolve the matricial system (FFT BANDS IMPLEMENTED)
      ! call solveU(u1,du1,1,myid)
      ! call solveV(u2,du2,2,myid)
      ! call solveU(u3,du3,3,myid)

      ! Resolve the matricial system (NO FFT BANDS)
      call solveU_W(u1,du1,u3,du3,a_ugrid,myid)
      call solveV(u2,du2,a_vgrid,myid)

    !   if(myid==0) then
    !     write(6,*) "du1 after solve", du1(0,:)
    !   end if 

      if(myid==0) then
        write(6,*) "t=", MPI_Wtime() - t1,"Finished Solving Velocities =====> Solving Pressure "
      end if 
      
     ! Compute the pressure gradient if constant mass flow condition is set
     if (flag_ctpress==0) then
       if (myid == 0) then
         call meanflow_ctU(u1,myid) ! Computes the convenient mpgx
       endif
     end if
      
      ! Solves pressure matrix
      call solveP(p,psi,u1,u2,u3,div,myid)

      ! Correct velocity for incompressibility
      call v_corr(u1,u2,u3,psi,div,myid,status,ierr)
    
    end do
    
    ! Evaluate mass flow if constant pressure gradient is set
    if (flag_ctpress==1) then
      if (myid == 0) then
        call meanflow_ctP(u1,myid)
      end if
    end if

    iter = iter+1

  end do 

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Deallocate memory and shut communications
  call finalize(u1,u2,u3,p,div,myid,status,ierr)

end program

subroutine divergence(div,u1,u2,u3,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!   DIVERGENCE   !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Divergence in Fourier space. (Pencils)
  ! Get the result of the divergence at the centers (pgrid)
  ! Div(u) =  du/dx +  dw/dz + dv/dy
  !        = i*kx*u + i*kz*w + FinDiff(v) 
  !        = i*kx*u + i*kz*w + (v(i) - v(i-1))/(DeltaTheta)*(dtheta/dyu)

  use declaration
  implicit none

  integer i,k,j,myid,column
  complex(8)  div( jlim(1,pgrid) : jlim(2,pgrid), columns_num(myid) )
  complex(8)   u1( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )
  complex(8)   u2( jlim(1,vgrid) : jlim(2,vgrid), columns_num(myid) )
  complex(8)   u3( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )

  complex(8) kx,kzF

  do column = 1,columns_num(myid)
      i   = columns_i(column,myid)
      k   = columns_k(column,myid)
      kx  = k1F_x(i)  ! im*kx
      kzF = k1F_z(k)  ! im*kz
      !!!!!!!!      First band, finite diffs      !!!!!!!!
      do j = jlim(1,pgrid),jlim(2,pgrid)
      div(j,column)= kx                 *  u1(j  ,column) &
  &                   +kzF                *  u3(j  ,column) &
  &                   +ddthetavi*dthdyu(j)*( u2(j  ,column) &
  &                                         -u2(j-1,column))
      end do
  end do

  end subroutine


  subroutine laplacian_U(Lu,u,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!    LAPLACIAN NEW CHECKED   !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Second-order Central Difference Scheme laplacian using 
  !   a 3-point stencil (Ferziger pp. 50).
  !   Notes:
  !     dy2i = 2/dy2; dy2(1,j) = [y(j  )-y(j-1)]*[y(j+1)-y(j-1)]  Upwind
  !                   dy2(2,j) = [y(j+1)-y(j  )]*[y(j  )-y(j-1)]  Centered
  !                   dy2(3,j) = [y(j+1)-y(j  )]*[y(j+1)-y(j-1)]  Backward

  use declaration
  implicit none

  integer i,k,j,column,myid
  complex(8) Lu(jlim(1,ugrid):jlim(2,ugrid),columns_num(myid))
  complex(8)  u(jlim(1,ugrid):jlim(2,ugrid),columns_num(myid))
  real(8) k2x,k2z


  do column = 1,columns_num(myid)
      i = columns_i(column,myid)
      k = columns_k(column,myid)
      k2x = k2F_x(i)
      k2z = k2F_z(k)
      do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      Lu(j,column) =         dyu2i(1,j) *u(j-1,column) &
  &                  +(k2x+k2z+dyu2i(2,j))*u(j  ,column) &
  &                  +         dyu2i(3,j) *u(j+1,column)
      end do

      ! if (column == 1 .and. myid == 0) then
      !   write(6,*) 'Lu complex, j=1..10:'
      !   do j = 1, 10
      !     write(6,*) real(Lu(j,1)), aimag(Lu(j,1))
      !   end do
      ! end if


  end do




  end subroutine

  subroutine laplacian_V(Lu,u,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!    LAPLACIAN  NEW CHECKED  !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Second-order Central Difference Scheme laplacian using 
  !   a 3-point stencil (Ferziger pp. 50).
  !   Notes:
  !     dy2i = 2/dy2; dy2(1,j) = [y(j  )-y(j-1)]*[y(j+1)-y(j-1)]  Upwind
  !                   dy2(2,j) = [y(j+1)-y(j  )]*[y(j  )-y(j-1)]  Centered
  !                   dy2(3,j) = [y(j+1)-y(j  )]*[y(j+1)-y(j-1)]  Backward

  use declaration
  implicit none

  integer i,k,j,column,myid
  complex(8) Lu(jlim(1,vgrid):jlim(2,vgrid),columns_num(myid))
  complex(8)  u(jlim(1,vgrid):jlim(2,vgrid),columns_num(myid))
  real(8) k2x,k2z

  do column = 1,columns_num(myid)
      i = columns_i(column,myid)
      k = columns_k(column,myid)
      k2x = k2F_x(i)
      k2z = k2F_z(k)
      do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
      Lu(j,column) =         dyv2i(1,j) *u(j-1,column) &
  &                  +(k2x+k2z+dyv2i(2,j))*u(j  ,column) &
  &                  +         dyv2i(3,j) *u(j+1,column)
      end do

      ! if (column == 1 .and. myid == 0) then
      !   write(6,*) 'Lu complex, j=1..10:'
      !   do j = 1, 10
      !     write(6,*) real(Lu(j,1)), aimag(Lu(j,1))
      !   end do
      ! end if

  end do

  end subroutine

  subroutine error(A,myid,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!     ERROR      !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    include 'mpif.h'             ! MPI variables
    integer ierr

    integer j,column,myid
    complex(8), intent(in) :: A( jlim(1,pgrid): jlim(2,pgrid), columns_num(myid) )
    real(8) erri,errband

    erri = 0d0
    errband = 0d0
    do column = 1,columns_num(myid)
        do j = jlim(1,pgrid),jlim(2,pgrid)
        err     = abs(A(j,column))
        errband = max(errband,err)

        ! write(6,*) "errband", errband
        end do
    end do
    erri = max(erri,errband)


    !if (myid == 0) then
    ! write(6,*) 'error(): size(A%f,1:2)=', size(A%f,1), size(A%f,2)
    ! write(6,*) 'error(): jlim(1:2,pgrid)=', jlim(1,pgrid), jlim(2,pgrid)
    ! write(6,*) 'error(): columns_num(myid)=', columns_num(myid), ' myid=', myid
    !end if


    call MPI_ALLREDUCE(erri,err,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)



  end subroutine


  subroutine finalize(u1,u2,u3,p,div,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!    FINALIZE    !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer iproc,j
  real(8) val
  ! complex(8), intent(inout) :: u1(jlim(1,ugrid):,:), u2(jlim(1,vgrid):,:), u3(jlim(1,ugrid):,:)
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  complex(8) :: div( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )


  iwrite = iwrite+nwrite
  call error(div,myid,ierr)

  u1PL = 0d0
  u2PL = 0d0
  u3PL = 0d0
  ppPL = 0d0
  call modes_to_planes_UVP(u1PL,u1,2,nyu,nyu_LB,myid,status,ierr)
  call modes_to_planes_UVP(u2PL,u2,1,nyv,nyv_LB,myid,status,ierr)
  call modes_to_planes_UVP(u3PL,u3,2,nyu,nyu_LB,myid,status,ierr)
  call modes_to_planes_UVP(ppPL, p,3,nyp,nyp_LB,myid,status,ierr)

  call record_out(u1,myid)

  if (myid/=0) then
      call MPI_SEND(1d0,1,MPI_REAL8,0,101+myid,MPI_COMM_WORLD,ierr)
  else
      do iproc=1,np-1
      call MPI_RECV(val,1,MPI_REAL8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
      write(*,*) 'process',iproc,'over'
      call MPI_SEND(1d0,1,MPI_REAL8,iproc,100,MPI_COMM_WORLD,ierr)
      end do
  end if
  if (myid/=0) then
      call MPI_RECV(val,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  else
      write(*,*) 'process',0,'over'
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid==0) then
      write(*,*) ''
      write(*,*) '---------------LITTLE HARSH---------------'
      write(*,*) "         'That wasn't that harsh'         "
      write(*,*) '------------------------------------------'
      write(*,*) ''
  end if
  call MPI_FINALIZE(ierr)

  end subroutine

  subroutine flowrateIm(Qu,u)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!         FLOW RATE          !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Flow rate: integrating 1st Fourier mode over the whole channel

  use declaration
  implicit none
  integer j
  real(8) Qu
  complex(8) u(nyu_LB:nyu+1)

  Qu = 0d0
  do j = 1,nn+1
      !Qu = Qu + .5d0*(yu(j+1) - yu(j-1))*dreal(u(j) !Old collocated version
      Qu = Qu + (yv(j) - yv(j-1))*dreal(u(j)) !Mixed ugrid vgrid
  end do

  end subroutine
    

  subroutine flowrateRe(Qu,u)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!        FLOW RATE 0         !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j
  real(8) Qu
  real(8) u(1:nn+2)

  Qu = 0d0
  do j = 1,nn+1
      !Qu = Qu + .5d0*(yu(j+1)-yu(j-1))*u(j)
      Qu = Qu + (yv(j) - yv(j-1))*u(j) !Mixed ugrid vgrid
  end do

  end subroutine

  subroutine flowrate_corr(u,mpg,g)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!        FLOW RATE 0         !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  ! use tridLU_3D
  implicit none
  integer j
  real(8) Qx0,Qx1,bdt,g,mpg
  complex(8) u(nyu_LB:nyu+1)

  !!!!!!!!!   dummy u11:  !!!!!!!!!
  ! Super meeeean!!
  bdt = -dRK(kRK)*dt
  u11(0) = 0d0
  do j = 1,nn+1
      u11(j) = bdt
  end do
  u11(nn+2) = 0d0
  call LUsol0(u11,0,nn+2)
  !!!!!!!!! u11 flowrate: !!!!!!!!!
  call flowrateRe(Qx1,u11)
  !!!!! uncorrected flowrate: !!!!!
  call flowrateIm(Qx0,u  )
  !!!!!  pressure correction: !!!!!
  g   = (QxT-Qx0)/Qx1
  mpg = mpg+g
  !!!!!  velocity correction: !!!!!
  do j = 0,nn+2
      u(j) = u(j)+g*u11(j)
  end do

  end subroutine

  subroutine maxvel(u)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!          MAX VEL           !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  complex(8) u(nyu_LB:nyu+1)
  integer j

  Umax = 0d0
  do j = 0,nn+2
      Umax = max(Umax,dreal(u(j)))
  end do

  ! do j= 1, 152
  !   write(6,*) " Umax", Umax, dreal(u(j)), "j", j
  ! end do 

  end subroutine

  subroutine meanflow_ctP(u1,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! MEAN FLOW WITH CONSTANT MPG  !!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer      :: myid
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )

  ! Mean (mode 0,1) is in proc 0 column 1
  if (myid == 0) then 
      ! call flowrateIm(Qx,u1(jlim(1,ugrid),1)) 
      call flowrateIm(Qx,u1(jlim(1,ugrid),1)) 
  end if

  end subroutine

  subroutine meanflow_ctU(u1,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!       MEAN FLOW CORR       !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,myid
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )

  ! Mean (mode 0,1) is in proc 0 column 1
  if (myid == 0) then 
      !call flowrate_corr(u1(jlim(1,ugrid),1),mpgx,dgx)
      call flowrate_corr(u1(jlim(1,ugrid),1),mpgx,dgx)
  end if

  end subroutine

  !subroutine meanpressgrad(mpg,u)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!   MEAN PRESSURE GRADIENT   !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  use declaration
  !  implicit none
  !  real(8) mpg,du0,du1
  !  complex(8) u(N(4,0):N(4,nband)+1)

  !  du0= dreal(u(1) )/(yu(1   )-yu(0 ))
  !  du1=-dreal(u(nn))/(yu(nn+2)-yu(nn+1))

  !  mpg=(du1-du0)/(Ly*Re)

  !end subroutine

  subroutine solveP(p,psi,u1,u2,u3,div,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!     SOLVE P    !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  ! use tridLU_3D
  implicit none
  integer :: i,k,j,column,myid

  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  complex(8) :: psi( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  complex(8) :: div( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )

  real(8) :: dtirk
  
  dtirk = dti/(aRK(kRK)+bRK(kRK))

  call divergence(div,u1,u2,u3,myid)
      
  do column = 1,columns_num(myid)
      do j = jlim(1,pgrid),jlim(2,pgrid)
      psi(j,column) = dtirk*div(j,column) 
      end do
  end do 
  
  !Boundary condition for pressure
  if (myid==0) then
      psi(jlim(1,pgrid),1) = 0d0 !C! Psi (mean) at bottom of channel = 0
  end if

  ! For modified wavenumbers, need BC on last pressure mode (as k2x=k2z=0)
  ! do column = 1,columns_num(iband,myid)
  ! i = columns_i(column,iband,myid)
  ! k = columns_k(column,myid)
  ! if(abs(k1F_x(i)*k1F_x(i))<10e-10.and.abs(k1F_z(k)*k1F_z(k))<10e-10)then
  ! psi(iband)%f(jlim(1,pgrid,iband),column) = 0d0
  ! endif
  ! enddo
      

  call LUsolP(psi,myid,jlim(1,pgrid),jlim(2,pgrid))

  do column = 1,columns_num(myid) 
      do j = jlim(1,pgrid),jlim(2,pgrid)
      p(j,column) = p(j,column)+psi(j,column)
      enddo
  enddo
  
  ! end do

  end subroutine

  subroutine RHS0_u1(du1,u1,Nu1,p,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!   RHS U1 NEW CHECKED  but not the mean p grad bit !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,column,myid
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: du1( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: Nu1(jlim(1,ugrid)+1     :jlim(2,ugrid)-1,     columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  real(8) C1,C2
  complex(8) C3

  !C1 = (aRK(kRK)+bRK(kRK))/Re !For solving for du
  C1 = aRK(kRK)/Re !For solving for u
  C2 = -cRK(kRK)

  ! write(6,*) "C1", C1, "C2", C2, myid

!   if (myid == 0) then
!     write(6,*) 'du1, j=1..10, col=2'
!     do j = 0, 10
!         write(6,*) du1(j,2)
!     end do
!   end if


  call laplacian_U(du1,u1,myid)
  do column = 1,columns_num(myid)
      i = columns_i(column,myid)
      C3 = -dRK(kRK)*k1F_x(i)
      do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      du1(j,column) = C1*du1(j,column) &
  &                              + C2*Nu1(j,column) &
  &                              + C3* p (j,column)
      end do
      
  end do


  !Apply mean pressure gradient int x
  ! mode (0,1)
  ! TODO TRICK: We know that mode(0,1) is the first column of proc 0
  if (myid==0) then
      C1 = -dRK(kRK)*mpgx
      do j=nyu_LB+1,nyu
      !do j=1,nn+1              !       In case you want not to apply pressure gradient in the immersed boundary
      du1(j,1) = du1(j,1) + C1        ! check why they had midband here 
      end do
  end if

    ! if (myid == 0) then
    !     write(6,*) 'du1, j=1..10, col=2'
    !     do j = 0, 10
    !         write(6,*) du1(j,2)
    !     end do
    ! end if


    ! if (myid == 0) then
    !     write(6,*) 'u1, j=1..10, col=2'
    !     do j = 0, 10
    !         write(6,*) u1(j,2)
    !     end do
    ! end if


  end subroutine

  subroutine RHS0_u2(du2,u2,Nu2,p,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!     RHS  U2    !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,column,myid
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: du2( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: Nu2(jlim(1,vgrid)+1:jlim(2,vgrid)-1,columns_num(myid))
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  real(8) C1,C2,C3

  !C1 = (aRK(kRK)+bRK(kRK))/Re !For solving for du
  C1 = aRK(kRK)/Re !For solving for u
  C2 = -cRK(kRK)

  ! write(6,*) "C1", C1, "C2", C2, myid

  
  call laplacian_V(du2,u2,myid)
  do column = 1,columns_num(myid)
      do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
      C3 = -dRK(kRK)*ddthetavi*dthdyv(j)
      du2(j,column) = C1* du2(j  ,column) &
  &                              + C2* Nu2(j  ,column) &
  &                              + C3*( p (j+1,column) & ! centeres to faces --> (j+1)-(j)
  &                                    -p (j,  column))
      end do

      ! if (column == 1 .and. myid == 0) then
      !   write(6,*) 'U2'
      !   write(6,*) 'du1, j=1..10:'
      !   do j = 1, 10
      !     write(6,*) real(du2%f(j,column)), aimag(du2%f(j,column))
      !   end do
      
      ! end if

      ! if (column == 1 .and. myid == 0) then
      !   write(6,*) 'du1, j=1..10:'
      !   do j = 1, 10
      !     write(6,*) 'nu2, p and p '
      !     write(6,*) Nu2%f(j  ,column), p%f(j+1,column),  -p%f(j,  column)
      !   end do

      ! end if

  end do
  ! end do
  
  end subroutine

  subroutine RHS0_u3(du3,u3,Nu3,p,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!     RHS  U3    !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,column,myid
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: du3( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: Nu3(jlim(1,ugrid)+1     : jlim(2,ugrid)-1,    columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  real(8) C1,C2
  complex(8) C3

  !C1 = (aRK(kRK)+bRK(kRK))/Re !For solving for du
  C1 = aRK(kRK)/Re !For solving for u
  C2 = -cRK(kRK)

  call laplacian_U(du3,u3,myid)
  do column = 1,columns_num(myid)
      k = columns_k(column,myid)
      C3 = -dRK(kRK)*k1F_z(k)
      do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      du3(j,column) = C1*du3(j,column) &
  &                              + C2*Nu3(j,column) &
  &                              + C3* p (j,column)
      end do

      ! if (column == 1 .and. myid == 0) then
      !   write(6,*) 'U3'
      !   write(6,*) 'du3, j=1..10:'
      !   do j = 1, 10
      !     write(6,*) real(du3%f(j,column)), aimag(du3%f(j,column))
          
      !   end do
      ! end if

      ! if (column == 1 .and. myid == 0) then
      !   write(6,*) 'du3, j=1..10:'
      !   do j = 1, 10
      !     write(6,*) 'nu2, p'
      !     write(6,*) Nu3%f(j  ,column), p%f(j,  column)
      !   end do

      ! end if
      
  end do

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!  NO FFT BANDS  !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine solveU_W(u,du,w,dw,a_ugrid,myid) !pass (u,du,w,dw)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!    SOLVE U1    !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  ! use tridLU_3D

  implicit none

  integer i,k,j,iband,column,myid
  complex(8) :: u ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: w ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: du( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: dw( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  real(8)  a_ugrid(3,jlim(1,ugrid):jlim(2,ugrid))

      do column = 1,columns_num(myid)
          i = columns_i(column,myid)
          k = columns_k(column,myid)
          du(jlim(1,ugrid),column) = 0d0
          dw(jlim(1,ugrid),column) = 0d0

          do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
          du(j,column) = u(j,column)+dt*(du(j,column)) !For solving for u
          dw(j,column) = w(j,column)+dt*(dw(j,column)) !For solving for w
          end do

          du(jlim(2,ugrid),column) = 0d0
          dw(jlim(2,ugrid),column) = 0d0
      end do

      ! they used to call immersed boundaries here !!!! how do we treat boundaries..? 
      ! idk where to look fo smooth wall

      call LUsolU_W(u,du,w,dw,a_ugrid(1:3,jlim(1,ugrid):jlim(2,ugrid)),ugrid,myid)

  end subroutine

  subroutine solveV(u,du,a_vgrid,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!    SOLVE U1    !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use declaration
  ! use tridLU_3D
  implicit none
  integer i,k,j,iband,column,myid
  complex(8) :: u ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: du( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  real(8)    :: a_vgrid(3,jlim(1,vgrid):jlim(2,vgrid))

      do column = 1,columns_num(myid)
          i = columns_i(column,myid)
          k = columns_k(column,myid)
          du(jlim(1,vgrid),column) = 0d0
          do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
          !         du%f(j,column) = dt*(du%f(j,column)) !For solving for du
          du(j,column) = u(j,column)+dt*(du(j,column)) !For solving for u
          end do
          du(jlim(2,vgrid),column) = 0d0
      end do

      call LUsolV(u,du,a_vgrid(1:3,jlim(1,vgrid):jlim(2,vgrid)),vgrid,myid)

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!   FFT BANDS  (SOLVE AND SOLVEV)  !!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! subroutine solveU(u,du,ufield,myid)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!    SOLVE U1    !!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !   use declaration
  !   implicit none

  !   integer i,k,j,iband,column,myid,ufield
  !   type(cfield)  u
  !   type(cfield) du

  !   !do iband = sband,eband
  !   do column = 1,columns_num(myid)
  !     i = columns_i(column,myid)
  !     k = columns_k(column,myid)
  !     du%f(jlim(1,ugrid),column) = 0d0
  !     do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
  ! !         du%f(j,column) = dt*(du%f(j,column)) !For solving for du
  !       du%f(j,column) = u%f(j,column)+dt*(du%f(j,column)) !For solving for u
  !     end do
  !     du%f(jlim(2,ugrid),column) = 0d0
  !   end do
  !   ! end do

  !   call LUsolV(du,ugrid,ufield,myid)
  ! !call LUsolV(du,ugrid,myid) !Can use for smooth channel

  ! !   do iband = sband,eband
  ! !     do column = 1,columns_num(myid)
  ! !       do j = jlim(1,ugrid),jlim(2,ugrid)
  ! !         u%f(j,column) = u%f(j,column)+du%f(j,column) !For solving for du
  ! !       enddo
  ! !     enddo
  ! !   end do

  !   ! do iband = sband,eband
  !   do column = 1,columns_num(myid)
  !     do j = jlim(1,ugrid),jlim(2,ugrid)
  !       u%f(j,column) = du%f(j,column) !For solving for u
  !     enddo
  !   enddo
  !   ! end do

  ! end subroutine


  ! subroutine solveV(u,du,ufield,myid)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!    SOLVE U1    !!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !   use declaration
  !   implicit none
  !   integer i,k,j,column,myid,ufield
  !   type(cfield)  u
  !   type(cfield) du


  !   ! do iband = sband,eband
  !   do column = 1,columns_num(myid)
  !     i = columns_i(column,myid)
  !     k = columns_k(column,myid)
  !     du%f(jlim(1,vgrid),column) = 0d0
  !     do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
  ! !         du(iband)%f(j,column) = dt*(du(iband)%f(j,column)) !For solving for du
  !       du%f(j,column) = u%f(j,column)+dt*(du%f(j,column)) !For solving for u
  !     end do
  !     du%f(jlim(2,vgrid),column) = 0d0
  !   end do
  !   ! end do

  !   call LUsolV(du,vgrid,ufield,myid)

  ! !   do iband = sband,eband
  ! !     do column = 1,columns_num(iband,myid)
  ! !       do j = jlim(1,vgrid,iband),jlim(2,vgrid,iband)
  ! !         u(iband)%f(j,column) = u(iband)%f(j,column)+du(iband)%f(j,column) !For solving for du
  ! !       enddo
  ! !     enddo
  ! !   end do

  !   ! do iband = sband,eband
  !   do column = 1,columns_num(myid)
  !     do j = jlim(1,vgrid),jlim(2,vgrid)
  !       u%f(j,column) = du%f(j,column) !For solving for u
  !     enddo
  !   enddo
  !   ! end do

  ! end subroutine




  subroutine v_corr(u1,u2,u3,psi,div,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!     V_CORR     !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  
  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,column
  complex(8) kx,kz,dtrk,dtrk_u2
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: div( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  complex(8) :: psi( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  real (8),pointer :: vcorrPL(:,:,:)
  
  real(8) weighting
  
  
  allocate(vcorrPL(Nspec_x+2,Nspec_z,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  !!!!!!!!!      u1:      !!!!!!!!!
  dtrk = dt*(aRK(kRK)+bRK(kRK))
  

  do column = 1,columns_num(myid)
      i  = columns_i(column,myid)
      kx = dtrk*k1F_x(i)
      
      do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
      u1(j,column) = u1(j,column)-kx*psi(j,column)
      end do

  end do


  do column = 1,columns_num(myid)
      u1(jlim(1,ugrid),column) = gridweighting_bc_u1*u1(jlim(1,ugrid)+1,column)
  enddo

  
  ! do iband = 2,3
  do column = 1,columns_num(myid)
      u1(jlim(2,ugrid),column) = gridweighting_bc_u1*u1(jlim(2,ugrid)-1,column)
  enddo
  ! enddo
  
  !!!!!!!!!      u2:      !!!!!!!!!
  dtrk_u2 = dtrk*ddthetavi

  do column = 1,columns_num(myid)
      
      do j = jlim(1,vgrid)+1,jlim(2,vgrid)-1
      u2(j,column) = u2(j,column)-dtrk_u2*dthdyv(j)*(psi(j+1,column)-psi(j,column))
      end do
  !       u2(iband)%f(jlim(1,vgrid),column) = 0d0
  !       u2(iband)%f(jlim(2,vgrid),column) = 0d0

      u2(jlim(1,vgrid),column) = 0d0
      u2(jlim(2,vgrid),column) = 0d0

  end do

  !!!!!!!!!      u3:      !!!!!!!!!

  do column = 1,columns_num(myid)
      k  = columns_k(column,myid)
      kz = dtrk*k1F_z(k)
      
      do j=jlim(1,ugrid)+1,jlim(2,ugrid)-1
      u3(j,column) = u3(j,column)-kz*psi(j,column)
      end do

  end do


  ! do iband = 1,2
  do column = 1,columns_num(myid)
      u3(jlim(1,ugrid),column) = gridweighting_bc_u3*u3(jlim(1,ugrid)+1,column)
  enddo
  ! enddo
  
  ! do iband = 2,3
  do column = 1,columns_num(myid)
      u3(jlim(2,ugrid),column) = gridweighting_bc_u3*u3(jlim(2,ugrid)-1,column)
  enddo

  
  
  !!!!!!!!!  divergence:  !!!!!!!!!

  call divergence(div,u1,u2,u3,myid)

  deallocate(vcorrPL)
  end subroutine