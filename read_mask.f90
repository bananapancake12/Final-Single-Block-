program read_mask
    implicit none

    integer :: Ny, nx, nz
    integer :: unit
    integer, allocatable :: maskU(:,:,:), maskV(:,:,:), maskW(:,:,:)
    integer :: j

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! U !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    unit = 10

    open(unit=unit, file='maskU.dat', form='unformatted', access='stream', status='old')

    ! Read header
    read(unit) Ny
    read(unit) nx
    read(unit) nz

    print *, "Ny =", Ny
    print *, "nx =", nx
    print *, "nz =", nz

    ! Allocate mask array
    allocate(maskU(nx,nz,Ny))

    ! Read full mask block
    read(unit) maskU

    close(unit)

    print *, "Read U complete."

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! V !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    unit = 10

    open(unit=unit, file='maskV.dat', form='unformatted', access='stream', status='old')

    ! Read header
    read(unit) Ny
    read(unit) nx
    read(unit) nz

    print *, "Ny =", Ny
    print *, "nx =", nx
    print *, "nz =", nz

    ! Allocate mask array
    allocate(maskV(nx,nz,Ny))

    ! Read full mask block
    read(unit) maskV

    close(unit)

    print *, "Read V complete."

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! W !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    unit = 10

    open(unit=unit, file='maskV.dat', form='unformatted', access='stream', status='old')

    ! Read header
    read(unit) Ny
    read(unit) nx
    read(unit) nz

    print *, "Ny =", Ny
    print *, "nx =", nx
    print *, "nz =", nz

    ! Allocate mask array
    allocate(maskW(nx,nz,Ny))

    ! Read full mask block
    read(unit) maskW

    close(unit)

    print *, "Read W complete."



    

    ! Print some sample values
    !print *, "maskU(1,5,34) =", maskU(:,5,34)
    print *, "maskV(1,5,34) =", maskV(:,7,44)
    ! print *, "maskW(1,5,34) =", maskW(:,5,34)

end program read_mask