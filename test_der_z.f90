module declaration
  implicit none
  integer, parameter :: Nspec_x = 4
  integer, parameter :: Nspec_z = 2
end module declaration

program test_der_z_exact
  use declaration
  implicit none

  real(8)    :: u   (Nspec_x+2, Nspec_z)
  real(8)    :: dudz(Nspec_x+2, Nspec_z)
  complex(8) :: kz(1:Nspec_z)

  integer :: i, k

  ! Fill u with obvious numbers so you can track packing:
  ! u(1,k), u(2,k) -> mode 0, etc (as YOUR routine currently uses)
  do k = 1, Nspec_z
    do i = 1, Nspec_x+2
      u(i,k) = 10.0d0*k + dble(i)
    end do
  end do

  ! Example kz: pure imaginary i*1, i*2
  kz(1) = (0.0d0, 1.0d0)
  kz(2) = (0.0d0, 2.0d0)

  dudz = -999.0d0

  print *, "==== Input u (real packed) ===="
  do k = 1, Nspec_z
    print *, "k=", k, " u(:,k)=", u(:,k)
  end do

  print *
  print *, "==== kz ===="
  do k = 1, Nspec_z
    print *, "k=", k, " kz=", kz(k)
  end do

  call der_z_N(u, dudz, kz)

  print *
  print *, "==== Output dudz (real packed) ===="
  do k = 1, Nspec_z
    print *, "k=", k, " dudz(:,k)=", dudz(:,k)
  end do

contains

  subroutine der_z_N(u,dudz,kz)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    integer i,k,dk2,k2, ip, kp
    real(8) u   (0:Nspec_x+1,Nspec_z)
    real(8) dudz(0:Nspec_x+1,Nspec_z)
    complex(8) kz(1:Nspec_z)

    complex(8), allocatable :: u_c(:,:), dudz_c(:,:)
    allocate(u_c(0:Nspec_x/2,Nspec_z))
    allocate(dudz_c(0:Nspec_x/2,Nspec_z))

    ! convering real to complex
    do kp = 1,Nspec_z
      do ip = 0,Nspec_x/2
        u_c(ip,kp) = cmplx(u(ip*2,kp), u(ip*2+1,kp))
      end do
    end do

    ! doing the operation in complex space
      do k = 1,Nspec_z
        do i = 0,Nspec_x/2
          dudz_c(i,k) = kz(k)*u_c(i,k)
        end do
      end do

    ! converting dudz_c back to real to be compatible with rest of code
    do kp = 1, Nspec_z
      do ip = 0, Nspec_x/2
        dudz(2*ip,     kp) = dble( dudz_c(ip,kp) )
        dudz(2*ip + 1, kp) = dimag( dudz_c(ip,kp) )
      end do
    end do

  end subroutine der_z_N

end program test_der_z_exact
