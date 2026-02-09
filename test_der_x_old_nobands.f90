module declaration
  implicit none
  integer, parameter :: Ngal_x  = 6   ! i = 0..Ngal_x/2 = 0..3 in complex space
  integer, parameter :: Ngal_z  = 4
  integer, parameter :: Nspec_x = 4   ! kept x modes: 0..Nspec_x/2 = 0..2
  integer, parameter :: Nspec_z = 2   ! kept z band size; Nspec_z/2 = 1
end module declaration

program test_der_x_old_nobands
  use declaration
  implicit none

  complex(8) :: u   (0:Ngal_x/2, Ngal_z)
  complex(8) :: dudx(0:Ngal_x/2, Ngal_z)
  complex(8) :: kx(0:Nspec_x/2)

  integer :: i, k

  ! ---- Fill u with the SAME pattern as your packed-real test ----
  ! In your packed-real test for each k:
  !   mode0 = (100*k + 0) + i(100*k + 1)
  !   mode1 = (100*k + 2) + i(100*k + 3)
  !   mode2 = (100*k + 4) + i(100*k + 5)
  ! We'll build u(i,k) to match those complex values.

  do k = 1, Ngal_z
    u(0,k) = cmplx(100.0d0*k + 0.0d0, 100.0d0*k + 1.0d0, kind=8)
    u(1,k) = cmplx(100.0d0*k + 2.0d0, 100.0d0*k + 3.0d0, kind=8)
    u(2,k) = cmplx(100.0d0*k + 4.0d0, 100.0d0*k + 5.0d0, kind=8)
    u(3,k) = cmplx(-999.0d0, -999.0d0, kind=8)  ! anti-alias x region slot (i=3)
  end do

  ! Same kx as your test: 0, i, 2i
  kx(0) = (0.0d0, 0.0d0)
  kx(1) = (0.0d0, 1.0d0)
  kx(2) = (0.0d0, 2.0d0)

  dudx = cmplx(-999.0d0, -999.0d0, kind=8)

  print *, "==== Input u (complex) ===="
  do k = 1, Ngal_z
    do i = 0, Ngal_x/2
      print *, "k=",k," i=",i," u=",u(i,k)
    end do
  end do

  print *
  print *, "==== kx ===="
  do i = 0, Nspec_x/2
    print *, "i=", i, " kx=", kx(i)
  end do

  call der_x_old_nobands(u, dudx, kx)

  print *
  print *, "==== Output dudx (complex) ===="
  do k = 1, Ngal_z
    do i = 0, Ngal_x/2
      print *, "k=",k," i=",i," dudx=",dudx(i,k)
    end do
  end do

contains

  subroutine der_x_old_nobands(u, dudx, kx)
    use declaration
    implicit none

    integer :: i, k
    complex(8) :: u   (0:Ngal_x/2, Ngal_z)
    complex(8) :: dudx(0:Ngal_x/2, Ngal_z)
    complex(8) :: kx(0:Nspec_x/2)

    ! This matches your ORIGINAL routine structure, but with
    ! N(1,iband) -> Nspec_x,  Ngal(1,iband) -> Ngal_x
    ! N(2,iband) -> Nspec_z,  Ngal(2,iband) -> Ngal_z

    do k = 1, Nspec_z/2
      do i = 0, Nspec_x/2
        dudx(i,k) = kx(i) * u(i,k)
      end do
      do i = Nspec_x/2+1, Ngal_x/2
        dudx(i,k) = (0.0d0, 0.0d0)
      end do
    end do

    do k = Nspec_z/2+1, Ngal_z - Nspec_z/2 + 1
      do i = 0, Ngal_x/2
        dudx(i,k) = (0.0d0, 0.0d0)
      end do
    end do

    do k = Ngal_z - Nspec_z/2 + 2, Ngal_z
      do i = 0, Nspec_x/2
        dudx(i,k) = kx(i) * u(i,k)
      end do
      do i = Nspec_x/2+1, Ngal_x/2
        dudx(i,k) = (0.0d0, 0.0d0)
      end do
    end do

  end subroutine der_x_old_nobands

end program test_der_x_old_nobands
