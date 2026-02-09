module declaration
  implicit none
  integer, parameter :: Ngal_x  = 6   ! so u has i = 0..Ngal_x-1 = 0..5
  integer, parameter :: Ngal_z  = 4
  integer, parameter :: Nspec_x = 4   ! so Nspec_x/2 = 2
  integer, parameter :: Nspec_z = 2   ! so Nspec_z/2 = 1
end module declaration

program test_der_x_exact
  use declaration
  implicit none

  real(8)    :: u   (0:Ngal_x-1, Ngal_z)
  real(8)    :: dudx(0:Ngal_x-1, Ngal_z)
  complex(8) :: kx(0:Nspec_x/2)

  integer :: i, k

  ! Fill packed-real u with obvious values
  do k = 1, Ngal_z
    do i = 0, Ngal_x-1
      u(i,k) = 100.0d0*k + dble(i)
    end do
  end do

  ! Simple kx values (pure imaginary like i*k), so dudx_c = kx*u_c will rotate pairs
  kx(0) = (0.0d0, 0.0d0)
  kx(1) = (0.0d0, 1.0d0)
  kx(2) = (0.0d0, 2.0d0)

  dudx = -999.0d0

  print *, "==== Input u (packed real) ===="
  do k = 1, Ngal_z
    print *, "k=", k, " u(:,k)=", u(:,k)
  end do

  print *
  print *, "==== kx ===="
  do i = 0, Nspec_x/2
    print *, "i=", i, " kx(i)=", kx(i)
  end do

  call der_x(u, dudx, kx)

  print *
  print *, "==== Output dudx (packed real) ===="
  do k = 1, Ngal_z
    print *, "k=", k, " dudx(:,k)=", dudx(:,k)
  end do

contains

  subroutine der_x(u,dudx,kx)
    use declaration
    implicit none

    integer i,k, ip, kp
    real(8) u   (0:Ngal_x-1,Ngal_z)
    real(8) dudx(0:Ngal_x-1,Ngal_z)

    complex(8) kx(0:Nspec_x/2)

    complex(8), allocatable :: u_c(:,:), dudx_c(:,:)
    allocate(u_c(0:Ngal_x/2,Ngal_z))
    allocate(dudx_c(0:Ngal_x/2,Ngal_z))

    u_c = 0d0
    dudx_c = 0d0

    ! converting real to complex
    do kp = 1,Ngal_z
      do ip = 0,Ngal_x/2
        u_c(ip,kp) = cmplx(u(ip*2,kp), u(ip*2+1,kp))
      end do
    end do

    ! the acc calculation in complex space
    do k = 1,Nspec_z/2
      do i = 0,Nspec_x/2
        dudx_c(i,k) = kx(i)*u_c(i,k)
      end do
      do i = Nspec_x/2+1,Ngal_x/2
        dudx_c(i,k) = 0d0
      end do
    end do

    do k = Nspec_z/2+1,Ngal_z-Nspec_z/2+1
      do i = 0,Ngal_x/2
        dudx_c(i,k) = 0d0
      end do
    end do

    do k = Ngal_z-Nspec_z/2+2,Ngal_z
      do i = 0,Nspec_x/2
        dudx_c(i,k) = kx(i )*u_c(i,k)
      end do
      do i = Nspec_x/2+1,Ngal_x/2
        dudx_c(i,k) = 0d0
      end do
    end do

    ! converting dudx_c back to real
    do kp = 1, Ngal_z
      do ip = 0, Ngal_x/2
        dudx(2*ip,     kp) = dble( dudx_c(ip,kp) )
        dudx(2*ip + 1, kp) = dimag( dudx_c(ip,kp) )
      end do
    end do

  end subroutine der_x

end program test_der_x_exact
