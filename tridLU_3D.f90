! *************************************************************************************************************** !
!
! Last modified by LD 09/01/2026
! - Solve u1 u2 u3 all combined
!
! *************************************************************************************************************** !
!
! Contains: LUsolU_W    - Called from littleharsh for solveU : calls LU_build and LU_dec
!                       - Solves the (1-L) matrix for the U and W velocities
!                       - Decomposed into (1-Lxz)(1-y)
!                       - Contains FFT to solve wall in physical space and interface in Fourier space
!           LUsolV    - Called from littleharsh for solveV : calls LU_build and LU_dec
!                       - Solves the (1-L) matrix for the V velocities
!                       - Decomposed into (1-Lxz)(1-y)
!                       - No FFT all Fourier
!           LUsolP    - Called from littleharsh for solveP : calls LU_buildP and LU_decP
!                       - Solves the lplacian of the pressure (divergance of the gradient)
!           LUsol0    - Called in flowrate_corr (part of const. flow rate) : calls LU_build0 and LU_dec0
!                       - Solves mean mode to calculate pressure gradient for constant mass flow rate
!           LU_build  - Called from Lusol
!                       - Builds the (1-L) matrices. Different for u/w and v
!           LU_buildP - Called from start (with allocation)
!                       - Builds the DG matrix for the pressure
!                       - Only called once as no dt dependence
!           LU_build0 - Called from Lusol0
!                       - Builds LHS for LUsol0
!           LU_dec    - Called from LUsol
!                       - LU Decomposition
!           LU_decP   - Called from LUsolP
!                       - LU Decomposition
!                       - Only called once as no dt dependence
!           LU_dec0   - Called from LUsol0
!                       - LU Decomposition
!
! *************************************************************************************************************** !

  
  subroutine LUsolU(u,rhsu,a,grid,myid)

  !----------------------------------------------------------------------*
    !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
    !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
    ! INPUT m_j number of unknowns
    !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
    !       a(2,:,i,k) inverse of the diagonal of L_j
    !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
    !       (note that the diagonal of U_j is 1)
    !       x(nystart:nyend,i,k)   rhs of the problem
    ! OUTPUT
    !       all untouched but the solution 'x'
    !----------------------------------------------------------------------*

    use declaration
    implicit none
    
    include 'mpif.h'             ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr

    integer    :: i,k,j,column,grid
    integer    :: iopt,myid
    ! complex(8), intent(inout) :: u(jlim(1,grid):,:), w(jlim(1,grid):,:)
    !complex(8), intent(inout) :: rhsu(jlim(1,grid):,:), rhsw(jlim(1,grid):,:)
    complex(8) :: u ( jlim(1,grid)        : jlim(2,grid),      columns_num(myid) )
    complex(8) :: rhsu( jlim(1,grid)      : jlim(2,grid),      columns_num(myid) )
      
    real(8)    :: a(3,jlim(1,grid):jlim(2,grid))


    ! -Original
        do column = 1,columns_num(myid)
          call LU_build(jlim(1,grid),jlim(2,grid),grid,column,myid,a)
          call LU_dec(jlim(1,grid),jlim(2,grid),a)

          u(jlim(1,grid),column)=rhsu(jlim(1,grid),column)*a(2,jlim(1,grid))
          
          do j = jlim(1,grid)+1,jlim(2,grid)
            u(j,column) = (rhsu(j,column)-a(1,j)*u(j-1,column))*a(2,j)
          end do

          do j = jlim(2,grid)-1,jlim(1,grid),-1
            u(j,column) =  u(j,column)-a(3,j)*u(j+1,column)
          end do
        end do

        ! if(myid==0) then 
        !   write(6,*) "solving u/w"
        ! end if 

        ! if(myid==0) then
        !   do column = 1,10
        !     write(6,*) u(2)%f(7,column)
        !   end do
        ! end if 
  end subroutine



  subroutine LUsolP(x,myid,nystart,nyend)
    !----------------------------------------------------------------------*
    !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
    !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
    ! INPUT m_j number of unknowns
    !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
    !       a(2,:,i,k) inverse of the diagonal of L_j
    !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
    !       (note that the diagonal of U_j is 1)
    !       x(nystart:nyend,i,k)   rhs of the problem
    ! OUTPUT
    !       all untouched but the solution 'x'
    !----------------------------------------------------------------------*

    use declaration
    implicit none
    integer i,k,j,column,myid
    integer nystart,nyend
    complex(8) x(nystart:nyend,columns_num(myid))
    
    do column = 1,columns_num(myid)
        i = columns_i(column,myid)
        k = columns_k(column,myid)
        x(nystart,column) = x(nystart,column)*DG(2,nystart,column)
        do j = nystart+1,nyend
          x(j,column) = (x(j,column)-DG(1,j,column)*x(j-1,column))*DG(2,j,column)
        end do
        do j = nyend-1,nystart,-1
          x(j,column) =  x(j,column)-DG(3,j,column)*x(j+1,column)
        end do

    end do


  end subroutine

  subroutine LUsol0(x,nystart,nyend)
    !----------------------------------------------------------------------*
    !                  LUsol for real (0,0) mode
    !----------------------------------------------------------------------*

    use declaration
    implicit none
    integer j
    integer nystart,nyend
    real(8) x(nystart:nyend)
    real(8), allocatable:: a(:,:)
    allocate(a(3,nystart:nyend))

    call LU_build0(nystart,nyend,a)
    call LU_dec  (nystart,nyend,a)

    x(nystart) = x(nystart)*a(2,nystart)
    do j = nystart+1,nyend
        x(j) = (x(j)-a(1,j)*x(j-1))*a(2,j)
    end do
    do j = nyend-1,nystart,-1
        x(j) = x(j)-a(3,j)*x(j+1)
    end do

    deallocate(a)

  end subroutine

  
  subroutine LU_build(nystart,nyend,grid,column,myid,a)
    !-------------------------------------------------------!
    !       specifies original values of a(1:3,j,i,k)       !
    !-------------------------------------------------------!

    use declaration
    implicit none
    integer i,k,j,grid,myid,column
    integer nystart,nyend 
    real(8) k2x,k2z,beta
    real(8) a(3,nystart:nyend)

    beta = bRK(kRK)*dt/Re
    
    if (grid==ugrid) then
    !!!!!!!!!!!!!!!!!!!     u velocity:       !!!!!!!!!!!!!!!!!!!
      i = columns_i(column,myid)
      k = columns_k(column,myid)
      k2x = k2F_x(i)
      k2z = k2F_z(k)

        a(1,nystart) = 0d0
        a(2,nystart) = 1d0  
        a(3,nystart) = gridweighting(1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0

      do j = nystart+1,nyend-1
        a(1,j) =    -beta* dyu2i(1,j)
        a(2,j) = 1d0-beta*(dyu2i(2,j)+k2x+k2z)
        a(3,j) =    -beta* dyu2i(3,j)
      end do

        a(1,nyend) = gridweighting(2)!0d0 !Free-shear -1 no-slip 0
        a(2,nyend) = 1d0
        a(3,nyend) = 0d0

    else if (grid==vgrid) then 
      !!!!!!!!!!!!!!!!!!!     v velocity:       !!!!!!!!!!!!!!!!!!!  
      i = columns_i(column,myid)
      k = columns_k(column,myid)
      k2x = k2F_x(i)
      k2z = k2F_z(k)

        a(1,nystart) = 0d0
        a(2,nystart) = 1d0  
        a(3,nystart) = 0d0

      do j = nystart+1,nyend-1
        a(1,j) =    -beta* dyv2i(1,j)
        a(2,j) = 1d0-beta*(dyv2i(2,j)+k2x+k2z)
        a(3,j) =    -beta* dyv2i(3,j)
      end do

        a(1,nyend) = 0d0
        a(2,nyend) = 1d0
        a(3,nyend) = 0d0
    end if 


  end subroutine


  subroutine LU_buildP(nystart,nyend,myid,a)
    !-------------------------------------------------------!
    !       specifies original values of a(1:3,j,i,k)       !
    !-------------------------------------------------------!

    use declaration
    implicit none
    ! real(8),allocatable:: a(:,:,:)
    real(8) a( 3, jlim(1,pgrid):jlim(2,pgrid), columns_num(myid) )
    integer column,i,k,j,myid
    integer nystart,nyend
    real(8) k2x,k2z

    real(8) D_vec_b(nystart-1:nyend+1)
    real(8) G_vec_b(nystart-1:nyend+1)
    real(8) D_vec_t(nystart-1:nyend+1)
    real(8) G_vec_t(nystart-1:nyend+1)

    do j = nystart,nyend
        D_vec_b(j-1) = -ddthetavi*dthdyu(j)
        G_vec_b(j  ) = -ddthetavi*dthdyv(j)
        D_vec_t(j  ) =  ddthetavi*dthdyu(j)
        G_vec_t(j+1) =  ddthetavi*dthdyv(j)
    end do
    
    
    
    
    do column = 1,columns_num(myid)
        i = columns_i(column,myid)
        k = columns_k(column,myid)
      
        !For exact wavenumbers
        k2x = k2F_x(i)
        k2z = k2F_z(k)
      
        !    !For 2nd order centrered difference wavenumbers
        !    k2x = k1F_x(i)*k1F_x(i)
        !    k2z = k1F_z(k)*k1F_z(k)
    
        a(2,nystart,column) = (D_vec_t(nystart)*G_vec_b(nystart  )) + k2x + k2z
        a(3,nystart,column) = (D_vec_t(nystart)*G_vec_t(nystart+1))
        a(2,nyend,  column) = (D_vec_b(nyend-1)*G_vec_t(nyend    )) + k2x + k2z
        a(1,nyend,  column) = (D_vec_b(nyend-1)*G_vec_b(nyend-1  ))
        do j = nystart+1,nyend-1
          a(2,j,column) = (D_vec_b(j-1)*G_vec_t(j  )) + (D_vec_t(j)*G_vec_b(j)) + k2x + k2z
          a(1,j,column) = (D_vec_b(j-1)*G_vec_b(j-1))
          a(3,j,column) = (D_vec_t(j  )*G_vec_t(j+1))
        end do

    end do
    
    if(myid==0)then
        a(1,nystart,1) = 0d0
        a(2,nystart,1) = 1d0/((yu(1+1)-yu(1)))**2 !C!
        a(3,nystart,1) = 0d0
    end if

    ! For modified wavenumbers, need BC on last pressure mode (as k2x=k2z=0)
    !  do column = 1,columns_num(myid)
    !  i = columns_i(column,myid)
    !  k = columns_k(column,iband,myid)
    !  if(abs(k1F_x(i)*k1F_x(i))<10e-10.and.abs(k1F_z(k)*k1F_z(k))<10e-10)then
    !  a(iband)%f_dg(1,nystart,column) = 0d0
    !  a(iband)%f_dg(2,nystart,column) = 1d0/((yu(1+1)-yu(1)))**2 !C!
    !  a(iband)%f_dg(3,nystart,column) = 0d0
    !  endif
    !  enddo
    

    call LU_decP(nystart,nyend,columns_num(myid),a(:,nystart:nyend,1:columns_num(myid)),myid)

    
  end subroutine

  subroutine LU_build0(nystart,nyend,a)
    !----------------------------------------------------------------------*
    !                  LU_build for real (0,0) mode
    !----------------------------------------------------------------------*

    use declaration
    implicit none
    integer j
    integer nystart,nyend
    real(8) a(3,nystart:nyend)
    real(8) beta

    !!!!!!!!!!!!!!!!!!     velocity     !!!!!!!!!!!!!!!!!!
    beta = bRK(kRK)*dt/Re
    a(1,nystart) = 0d0
    a(2,nystart) = 1d0
    a(3,nystart) = gridweighting(1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
    do j = nystart+1,nyend-1
        a(1,j) =    -beta*dyu2i(1,j)
        a(2,j) = 1d0-beta*dyu2i(2,j)
        a(3,j) =    -beta*dyu2i(3,j)
    end do
    a(1,nyend) = gridweighting(2)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
    a(2,nyend) = 1d0
    a(3,nyend) = 0d0

  end subroutine

  subroutine LU_dec(nystart,nyend,a)
    !----------------------------------------------------------------------*
    !         GIVEN 'a_j' FOR j=nystart:nyend
    !         PERFORMS a_j=L_j*U_j
    !         WITH 'a_j' allmost TRIDIAGONAL
    ! INPUT  nystart:nyend matrix leading dimension
    !        iopt selects original value of 'a'
    !        a_j is given as a tridiagonal with
    !        a(1,nystart:nyend-1) lower diagonal,
    !        a(2,:) main diagonal,
    !        a(3,nystart+1:nyend) upper diagonal
    !        the matrix 'a_j' is assembeld as:
    !  a_j(2,1) a_j(3,1)    0        0         0      .............       0
    !  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
    !     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
    !     0     ...................................................       0
    !     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
    !     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
    !
    ! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
    !         a(2,nystart:nyend) inverse of the diagonal of L_j
    !         a(3,nystart:nyend-1) upper diagonal of U_j
    !         (notice that the diagonal of U_j is 1)
    !----------------------------------------------------------------------*

    use declaration
    implicit none
    integer j,myid
    integer nystart,nyend
    real(8) a(3,nystart:nyend)
  
    a(2,nystart) = 1d0/a(2,nystart)
    do j = nystart+1,nyend
      a(3,j-1) = a(3,j-1)*a(2,j-1)
      a(2,j) = 1d0/(a(2,j)-a(1,j)*a(3,j-1))
    end do

  end subroutine
  !  
  subroutine LU_decP(nystart,nyend,columns,a,myid)
    !----------------------------------------------------------------------*
    !         GIVEN 'a_j' FOR j=nystart:nyend
    !         PERFORMS a_j=L_j*U_j
    !         WITH 'a_j' allmost TRIDIAGONAL
    ! INPUT  nystart:nyend matrix leading dimension
    !        iopt selects original value of 'a'
    !        a_j is given as a tridiagonal with
    !        a(1,nystart:nyend-1) lower diagonal,
    !        a(2,:) main diagonal,
    !        a(3,nystart+1:nyend) upper diagonal
    !        the matrix 'a_j' is assembeld as:
    !  a_j(2,1) a_j(3,1)    0        0         0      .............       0
    !  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
    !     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
    !     0     ...................................................       0
    !     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
    !     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
    !
    ! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
    !         a(2,nystart:nyend) inverse of the diagonal of L_j
    !         a(3,nystart:nyend-1) upper diagonal of U_j
    !         (notice that the diagonal of U_j is 1)
    !----------------------------------------------------------------------*

    use declaration
    implicit none

    integer j,column
    integer columns,myid
    integer nystart,nyend
    real(8) a(3,nystart:nyend,columns)
    
    do column = 1,columns
        a(2,nystart,column) = 1d0/a(2,nystart,column)
        do j = nystart+1,nyend
          a(3,j-1,column) = a(3,j-1,column)*a(2,j-1,column)
          a(2,j  ,column) = 1d0/(a(2,j,column) - a(1,j,column)*a(3,j-1,column))
        end do
    end do

  end subroutine
