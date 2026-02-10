! *************************************************************************************************************** !
!
! Last modified by C 12/05/2016
!   - Set up for SHS (fingers crossed)
!   
!   Notes to self: eveyrthing has a commented section below it thats the original code. 
!                  the new section has double comments- this is the stuff that was originally commented
!                  single commments in new secctions are stuff ive removed that arent useful ie textures but left so i know where it is incase  
!
! *************************************************************************************************************** !
!
!  Notes: - when making it single block, getini and mblockini were moved into init_mod.f90
!         - mblock_ini_parabolic_profile has not been edited much- ie it is still multiblock and structs r still there, need to fix inputs to planes to modes etc
!
! Contains:   start - start file - called from littleharsh : calls ygrid, def_k, proc_lims, getbounds
!                     - reads input file, sets up domain, allocates memory, shares info between procs
!             ygrid - Called from start
!                     - Creates (streched) y grid points (seperate for u/w and v points)
!                       - Different grid streching in channel and immersed boundaries
!                     - Calculates dtheta/dy transformation
!                     - Calculates dy2i for all points (2nd derrivative) 
!             def_k - Called from start
!                     - Defines the x and z wavenumbers and their squares
!                       - Also calculates modified wavenumber but not used...
!         getbounds - Called from start
!                     - Calls boundary file and sets up gemomerty for immersed boundaries
! proc_lims_columns - Called from start
!                     - Calculates column lists for each proc and jlims
!                     - Also calculates dk
!                     - Calculates SHS boundary conditions
!  proc_lims_planes - Called from start
!                     - Defiens plane limits to balance between procs
!                     - Defines jgal
!            getini - Called from littleharsh : 
!                     - Grabs the initial conditions from file
!                     - )Or sets up parabolic profile)
!        init_stats - Called from getini : 
!                     - Initialises statistics

!         
!             start  - done, mostly checked, LUBuild etc cold be a bit dodge 
!             ygrid  - done and chcked 
!             def_k  - done and checked 
!         getbounds  - commented out bc i believe we dont use but check...
! proc_lims_columns  - done and checked 
!  proc_lims_planes  - done, check the weighting 
!            getini  - done , needs checking check why flowrate used midband
!        init_stats  - not touched (check what ju is) otherwise all good i think 
!
!
! *************************************************************************************************************** !

subroutine start(myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    START NEW   !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  ! use init_mod
  ! use tridLU_3D
  ! use littleharsh_mod
  implicit none

  include 'mpif.h'            ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer iproc,inputInt(4),k, dimsInt(10)
  integer i
  real(8) inputR(12)           ! Message passing array containing input parameters



  iproc = 0
  do while (pnodes<np)
    iproc  = iproc+1
    pnodes = 2**iproc
  end do


! The master proc is in charge of: 
!  - reading input file
!  - initialise N, Ngal
!  - initialise the y-grid
!  and send this data to the other procs.
  if (myid==0) then
    open(40,file='input.in',form='formatted')
    do i = 1,9
      read(40,10)             ! Input file header
    end do
    read(40,10) Re            ! Reynolds number based on the bulk velocity and the channel half-height 
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10) alp           ! The streamwise periodicity Lx of the box is Lx=2*pi/alp
    read(40,10) bet           ! The spanwise periodicity Lz of the box is Lz=2*pi/bet
    Lx = 2d0*pi/alp
    Lz = 2d0*pi/bet
    Ly = 2d0

    ! Ngal is the physical space 
    ! Nspec is the Fourier space 
    ! They contain the same information except that N is 2/3Ngal in the number of x- and z-grid points,

    ! allocate(Ngal(4,0:4))
    ! Ngal = 0
    ! grid points in the 'x' direction
    ! read(40,20) Ngal(1,nband)
    read(40,20) Ngal_x
    read(40,20)

    ! grid points in the 'z' direction
    ! read(40,20) Ngal(2,nband)
    read(40,20) Ngal_z
    read(40,20)

    ! Ngal(1,nband) = Ngal_x
    ! Ngal(2,nband) = Ngal_z

    ! ! fill cols 1:4 with vals 
    ! Ngal(1,1:nband)= Ngal(1,nband)
    ! Ngal(2,1:nband)= Ngal(2,nband)

    ! grid points in the 'y' direction from 'y=-1' to 'y=1'
    !read(40,20) Ngal(3,nband)
    read(40,20) nyv

    !!!!!!!!!!!!!!!!!!!!!! New variable definitions !!!!!!!!!!!!!!!!!!!!
    nyv = nyv -2
    nyu =  nyv + 1 
    nyp = nyv 

    nyu_LB = 0
    nyv_LB = 0 
    nyp_LB = nyu_LB + 1

    Nspec_x = (2 * Ngal_x)/3
    Nspec_z = (2 * Ngal_z)/3


    write(*,*) 'nyu =', nyu, ' nyv =', nyv, ' nyp =', nyp
    write(*,*) 'nyu_LB =', nyu_LB, ' nyp_LB =', nyp_LB
    write(*,*) 'Nspec_x =', Nspec_x, ' Ngal_x =', Ngal_x
    write(*,*) 'Nspec_z =', Nspec_z, ' Ngal_z =', Ngal_z

    !!!!!!!!!!!!!!!!!!!!!! Creating N for Hist write out !!!!!!!!!!!!!!!!!!!!
    allocate(N(4,0:4))


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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Ngal(3,nband) = nyv
    ! !Ngal(3,nband) = Ngal(3,nband)-2
    ! Ngal(1,nband+1) = -2

    ! allocate(N(4,0:nband+1))
    ! allocate(Ny(3,0:nband+1))

    ! N = Ngal
   
    ! N(1,1:nband) = 2*(Ngal(1,nband)/3)
    ! N(2,1:nband) = 2*(Ngal(2,nband)/3)

    read(40,10) dyq              ! dymax/dymin
    read(40,20) ppp              ! y polinomial exponent
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    ! Initial conditions
    read(40,20) flag_init     ! specifies initial conditions: 1 new, 2 continuing, 3 parabolic
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,30) fnameimb  ! initial conditions file path, if flag_init=1
    read(40,10)               ! - Blank
    read(40,30) dirin         ! directory of the initial conditions
    read(40,30) dirout        ! directory where output files are saved
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    ! Parameters of the numerical method I
    read(40,10) CFL           ! CFL
    read(40,20) nwrite        ! frecuency of output in terminal, meassured in time steps
    read(40,21) flag_ctpress  ! sets constant flow rate (0) or pressure gradient (1)
    read(40,10)               ! - Blank
    if (flag_ctpress==1) then
      read(40,10) mpgx      ! PRESCRIBED MEAN PRESSURE GRADIENT
    else if (flag_ctpress==0) then
      read(40,10)
    else
      write(*,*) 'Wrong presure flag'
      stop
    end if
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    write(6,*) "reading parameters for numerical method"

    ! Parameters of the numerical method II
    read(40,10) !Kib       ! immersed boundaries forcing parameter
    read(40,10) maxerr    ! maximum error for convergence in pressure
    read(40,10) maxt      ! maximum simulation time 't'
    read(40,10) maxA      ! boundary condition stability precission limiter
    maxA = -maxA
    read(40,20) !physlim_bot      ! Defines the last physical space plane (bottom band) for tridiag solve
    read(40,20) !physlim_top      ! Defines the last physical space plane (top band) for tridiag solve
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    write(6,*) "reading surface geometry"

    ! Surface geometry
    read(40,20) geometry_type ! 0 smooth, 1 hydrophobic, 2 roughness
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,20) !ntilex        ! number of tiles in the periodic box
                              ! subroutine boundary checks that 'mod(Ngal(1,1)/ntilex) = 0'
    read(40,20) !ntilez        ! number of tiles in the periodic box
                              ! subroutine boundary checks that 'mod(Ngal(2,1)/ntilez) = 0'
    if (geometry_type /= 0) then
        write(*,*) 'TMust be smooth wall. Moron.'
        stop
    end if 
    read(40,10) !posth         ! riblet height to width ratio
    read(40,20) !npeakx        ! number of points representing the blade tip
    read(40,20) !npeakz        ! number of points representing the blade tip
    read(40,10) !Lfracx        ! number of points representing the blade tip
    read(40,10) !Lfracz        ! number of points representing the blade tip
    read(40,10) !bslpu1        ! number of points representing the blade tip
    read(40,10) !bslpu3        ! number of points representing the blade tip
    
  
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,30) dirlist       ! directory of the nonlinear interaction list
    read(40,30) heading       ! heading of list
    close(40)

    ! write(6,*) "call ygrid"

    ! Creates the geometry in the y-direction
    call ygrid

    ! write(6,*) "ygrid finished"

    allocate(u11(0:nn+2))

    write(ext1,'(i4.4)') Nspec_x       ! Writes the number of point in the x-direction in the first band into ext1
    write(ext2,'(i4.4)') Nspec_z       ! Writes the number of point in the z-direction in the first band into ext2
    write(ext3,'(i4.4)') nn+2         ! Writes the number of point in the y-direction                   into ext3
    write(ext4,'(i5.5)') int(10d0*t)  ! Writes the time                                                 into ext4

    filout     = 'boundary_'
    boundfname = trim(dirout)//trim(filout)//ext1//'x'//ext2//'x'//ext3//'.dat'     ! Example: boundary_192x360x182.dat
                                                                                    !  index('string','t') returns 2.
                                                                                    !  index(filout,' ') returns
                                                                                    !  the value of the first empty cell.
    filout     = 'hist_'
    fnameima   = trim(dirout)//trim(filout)//ext1//'x'//ext2//'x'//ext3//'.dat'     ! Example:     hist_0192x0360x0182.dat

    !open(30,file=fnameima,form='unformatted')
    !write(30) t,Re,alp,bet,mpgx,nband
    !write(30) N

    write(6,*) "writing hist file"
inquire(file=fnameima, exist=exist_file_hist)
if (exist_file_hist) then
open (30,file=fnameima,form='unformatted',status="old", position="append", action="write")
else
open (30,file=fnameima,form='unformatted',status="new", action="write")
write(30) t,Re,alp,bet,mpgx,nband
write(30) N
end if
    

    
    filout     = 'c_four_'
    fnameima   = trim(dirout)//trim(filout)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat' ! Example:   c_four_0192x0360x0182_t01554.dat

! The master thread send to the others all the information they'll need
    inputR( 1) = Re
    inputR( 2) = alp
    inputR( 3) = bet
    inputR( 4) = Lx
    inputR( 5) = Ly
    inputR( 6) = Lz
    inputR( 7) = dtheta
    inputR( 8) = dthetai
    inputR( 9) = CFL
    inputR(10) = maxt
    inputR(11) = mpgx
    inputR(12) = ddthetavi

    inputInt(1) = flag_init
    inputInt(2) = nwrite
    inputInt(3) = nn
    inputInt(4) = flag_ctpress



    ! packing new dimension variables to be sent 
    dimsInt(1)  = Ngal_x
    dimsInt(2)  = Ngal_z
    dimsInt(3)  = Nspec_x
    dimsInt(4)  = Nspec_z
    dimsInt(5)  = nyv
    dimsInt(6)  = nyu
    dimsInt(7)  = nyp
    dimsInt(8)  = nyv_LB
    dimsInt(9)  = nyu_LB
    dimsInt(10) = nyp_LB



    do iproc=1,np-1
      call MPI_SEND(inputR  ,                     12,MPI_REAL8  ,iproc,       iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(inputInt,                      4,MPI_INTEGER,iproc,  1000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dimsInt ,                     10,MPI_INTEGER,iproc,  1500+iproc,MPI_COMM_WORLD,ierr)
      ! call MPI_SEND(N       ,4*(nband+2)            ,MPI_INTEGER,iproc,  2000+iproc,MPI_COMM_WORLD,ierr)
      ! call MPI_SEND(Ngal    ,4*(nband+2)            ,MPI_INTEGER,iproc,  3000+iproc,MPI_COMM_WORLD,ierr)
      ! call MPI_SEND(Ny      ,3*(nband+2)            ,MPI_INTEGER,iproc,  4000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(yu      ,   nyu-nyu_LB+2 ,MPI_REAL8  ,iproc,  7000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dthdyu  ,   nyu-nyu_LB+2 ,MPI_REAL8  ,iproc,  8000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dyu2i   ,3*(nyu-nyu_LB+2),MPI_REAL8  ,iproc,  9000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(yv      ,   nyv-nyv_LB+2 ,MPI_REAL8  ,iproc, 10000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dthdyv  ,   nyv-nyv_LB+2 ,MPI_REAL8  ,iproc, 11000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dyv2i   ,3*(nyv-nyv_LB+2),MPI_REAL8  ,iproc, 12000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dirlist ,          120       , MPI_CHARACTER,iproc, 13000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(heading ,          120       , MPI_CHARACTER,iproc, 14000+iproc,MPI_COMM_WORLD,ierr)
    end do


  else

! All the procs, except the master, receive and store the data
! Every proc knows the physical parameters and the whole geometry describing the problem (including N and Ngal)
    call MPI_RECV(inputR  ,20,MPI_REAL8  ,0,     myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(inputInt,15,MPI_INTEGER,0,1000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dimsInt ,10,MPI_INTEGER,0,1500+myid,MPI_COMM_WORLD,status,ierr)

    Re        = inputR( 1)
    alp       = inputR( 2)
    bet       = inputR( 3)
    Lx        = inputR( 4)
    Ly        = inputR( 5)
    Lz        = inputR( 6)
    dtheta    = inputR( 7)
    dthetai   = inputR( 8)
    CFL       = inputR( 9)
    maxt      = inputR(10)
    mpgx      = inputR(11)
    ddthetavi = inputR(12)

    flag_init    = inputInt(1)
    nwrite       = inputInt(2)
    nn           = inputInt(3)
    flag_ctpress = inputInt(4)


    Ngal_x  = dimsInt(1)
    Ngal_z  = dimsInt(2)
    Nspec_x = dimsInt(3)
    Nspec_z = dimsInt(4)
    nyv     = dimsInt(5)
    nyu     = dimsInt(6)
    nyp     = dimsInt(7)
    nyv_LB  = dimsInt(8)
    nyu_LB  = dimsInt(9)
    nyp_LB  = dimsInt(10)

    ! allocate(N   (4,0:nband+1))
    ! allocate(Ngal(4,0:nband+1))
    ! allocate(Ny  (3,0:nband+1))
    ! call MPI_RECV(N      ,4*(nband+2)            ,MPI_INTEGER,0, 2000+myid,MPI_COMM_WORLD,status,ierr)
    ! call MPI_RECV(Ngal   ,4*(nband+2)            ,MPI_INTEGER,0, 3000+myid,MPI_COMM_WORLD,status,ierr)
    ! call MPI_RECV(Ny     ,3*(nband+2)            ,MPI_INTEGER,0, 4000+myid,MPI_COMM_WORLD,status,ierr)
    allocate(yu    (  nyu_LB:nyu+1))
    allocate(dthdyu(  nyu_LB:nyu+1))
    allocate(dyu2i (3,nyu_LB:nyu+1))
    allocate(yv    (  nyv_LB:nyv+1))
    allocate(dthdyv(  nyv_LB:nyv+1))
    allocate(dyv2i (3,nyv_LB:nyv+1))
    call MPI_RECV(yu     ,   nyu-nyu_LB+2 ,MPI_REAL8  ,0, 7000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dthdyu ,   nyu-nyu_LB+2 ,MPI_REAL8  ,0, 8000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dyu2i  ,3*(nyu-nyu_LB+2),MPI_REAL8  ,0, 9000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(yv     ,   nyv-nyv_LB+2 ,MPI_REAL8  ,0,10000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dthdyv ,   nyv-nyv_LB+2 ,MPI_REAL8  ,0,11000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dyv2i  ,3*(nyv-nyv_LB+2),MPI_REAL8  ,0,12000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dirlist,                120,MPI_CHARACTER  ,0,13000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(heading,                120,MPI_CHARACTER  ,0,14000+myid,MPI_COMM_WORLD,status,ierr)
  end if

! From here on, every proc do everything which follows.

  nstat = min(20,nwrite)

  ! Initialise the FFT

    allocate(buffR_x( 2*Nspec_x+18), &
&            buffC_z(10*Nspec_z+19), &
&            buffRal_x( 2*Ngal_x+18), &
&            buffCal_z(10*Ngal_z+19))

    buffR_x   = 0.0d0
    buffC_z   = 0.0d0
    buffRal_x = 0.0d0
    buffCal_z = 0.0d0

    ! write(*,*) '  size(buffR_x(iband)%b)   = ', size(buffR_x  %b)
    ! write(*,*) '  size(buffC_z(iband)%b)   = ', size(buffC_z  %b)
    ! write(*,*) '  size(buffRal_x(iband)%b) = ', size(buffRal_x%b)
    ! write(*,*) '  size(buffCal_z(iband)%b) = ', size(buffCal_z%b)

    call rfti(Nspec_x,buffR_x )
    call cfti(Nspec_z,buffC_z )
    call rfti(Ngal_x,buffRal_x)
    call cfti(Ngal_z,buffCal_z)


  ! Define wavenumbers k1F_x and k1F_z & look up list for nonlinear
  call def_k

  ! Runge-Kutta coefficients (Le&Moin, 1991)
  aRK(1) =   4d0/15d0    ! aRK = alpha
  aRK(2) =   1d0/15d0    ! aRK = alpha
  aRK(3) =   1d0/ 6d0    ! aRK = alpha
  bRK(1) =   4d0/15d0    ! bRK = beta
  bRK(2) =   1d0/15d0    ! bRK = beta 
  bRK(3) =   1d0/ 6d0    ! bRK = beta
  gRK(1) =   8d0/15d0    ! gRK = gamma
  gRK(2) =   5d0/12d0    ! gRK = gamma
  gRK(3) =   3d0/ 4d0    ! gRK = gamma
  cRK(1) =   0d0         ! cRK = zeta
  cRK(2) = -17d0/60d0    ! cRK = zeta
  cRK(3) = - 5d0/12d0    ! cRK = zeta
  dRK(1) =   8d0/15d0    ! dRK = gamma + zeta
  dRK(2) =   2d0/15d0    ! dRK = gamma + zeta
  dRK(3) =   1d0/ 3d0    ! dRK = gamma + zeta

  ! Runge-Kutta coefficients (Spalart, 1991)
  ! We need to check what is the right coefficient for the pressure
  !aRK(1) =  29d0/ 96d0    ! aRK = alpha
  !aRK(2) = - 3d0/ 40d0    ! aRK = alpha
  !aRK(3) =   1d0/  6d0    ! aRK = alpha
  !bRK(1) =  37d0/160d0    ! bRK = beta
  !bRK(2) =   5d0/ 24d0    ! bRK = beta 
  !bRK(3) =   1d0/  6d0    ! bRK = beta
  !gRK(1) =   8d0/ 15d0    ! gRK = gamma
  !gRK(2) =   5d0/ 12d0    ! gRK = gamma
  !gRK(3) =   3d0/  4d0    ! gRK = gamma
  !cRK(1) =   0d0          ! cRK = zeta
  !cRK(2) = -17d0/ 60d0    ! cRK = zeta
  !cRK(3) = - 5d0/ 12d0    ! cRK = zeta
  !dRK(1) =   8d0/ 15d0    ! dRK = gamma + zeta
  !dRK(2) =   2d0/ 15d0    ! dRK = gamma + zeta
  !dRK(3) =   1d0/  3d0    ! dRK = gamma + zeta

  err = 10d0*maxerr
  !dtv = Re/(-k2F_x(N(1,1)/2)-k2F_z(N(2,1)/2)+4d0/(yu(N(4,0)+1)-yu(N(4,0)))**2) !E! yv or yu??? Smaller timestep?....
! print *, "dtv yv"
  dtv = Re/(-k2F_x(Nspec_x/2)-k2F_z(Nspec_z/2)+4d0/(yv(nyv_LB+1)-yv(nyv_LB))**2) !E! yv or yu??? Smaller timestep?....

  ! write(*,*) 'k2F_x term = ', -k2F_x(N(1,nband)/2)
  ! write(*,*) 'k2F_z term = ', -k2F_z(N(2,nband)/2)
  ! write(*,*) 'yv(N(3,0)+1) = ', yv(N(3,0)+1)
  ! write(*,*) 'yv(N(3,0)  ) = ', yv(N(3,0))
  ! write(*,*) 'dy = yv(N(3,0)+1) - yv(N(3,0)) = ', yv(N(3,0)+1) - yv(N(3,0))
 

10 FORMAT(7X,D10.1)
20 FORMAT(7X,I10)
21 FORMAT(12X,I5)
30 FORMAT(8X,A50)

  ! Calculates the different variables that define planes (phys) and pencils (spectral) for every proc 
  ! Computes planelim, bandPL, crossband and dk 
  ! Computes sband and eband
  !write(6,*) "meow"
  call proc_lims_columns(myid)

  ! !! added nonlin read here so read before allocating proc_lims_planes

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid == 0) then
    write(*,*) 'Getting Weights'
    write(*,*) dirlist
    write(*,*)
    call get_weights
  end if
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! write(*,*) 'finished reading nonlinear interaction list'

  ! write(6,*) "proc lims planes"
  call proc_lims_planes (myid)
  ! write(6,*) "finished proc lims planes"

  
  !write(6,*) "doing grid weighting"

  ! Grid weighting for extrapolating and interpolating the ghost points at the wall 
  !(doesnt look like its used), also eqn doesnt make sense
  allocate(gridweighting(2))
  gridweighting(1) =-(yv(nyv_LB  )-yu(nyu_LB  ))/(yv(nyv_LB  )-yu(nyu_LB+1)) 
  gridweighting(2) = (yu(nyu+1)-yv(nyv+1))/(yv(nyv+1)-yu(nyu  )) 

  ! write(6,*) "gridweighting 1 and 2", gridweighting(1), gridweighting(2)

  gridweighting_bc_u1 = (yu(0) - yv(0)) / (yu(1) - yv(0))
  gridweighting_bc_u3 = (yu(0) - yv(0)) / (yu(1) - yv(0))

  ! write(6,*) 'gridweighting_bc_u1 after = ', gridweighting_bc_u1
  ! write(6,*) 'gridweighting_bc_u3 after = ', gridweighting_bc_u3





  ! Used in interp_v and v_corr
  allocate(gridweighting_interp(2))
  gridweighting_interp(1) = (yu(nyu_LB  )-yu(nyu_LB+1))/(yv(nyv_LB  )-yu(nyu_LB+1))
  gridweighting_interp(2) = (yu(nyu+1)-yu(nyu  ))/(yv(nyv+1)-yu(nyu  ))

  ! write(6,*) "gridweighting_interp", gridweighting_interp(1), gridweighting_interp(2)

  ! PL = PLane
  ! By default variables are stored in columns
  allocate( u1PL       (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u2PL       (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u3PL       (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u1PL_itp   (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u2PL_itp   (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u3PL_itp   (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( ppPL       (igal,kgal,jgal(pgrid,1)-1:jgal(pgrid,2)+1))
  allocate(Nu1PL       (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(Nu2PL       (igal,kgal,jgal(vgrid,1)  :jgal(vgrid,2)  ))
  allocate(Nu3PL       (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(Nu1PL_dy    (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(Nu2PL_dy    (igal,kgal,jgal(vgrid,1)  :jgal(vgrid,2)  ))
  allocate(Nu3PL_dy    (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(uu_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(uv_fPL      (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(uw_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(vu_fPL      (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(vv_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(vw_fPL      (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(wu_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(wv_fPL      (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(ww_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))

  ! du1dy_planes( nx, nz, jplanes of each MPI rank)
  
  allocate(du1dy_planes(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du2dy_planes(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du3dy_planes(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))

  
  allocate(du1dy_planes2(Ngal_x+2,Ngal_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du2dy_planes2(Ngal_x+2,Ngal_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du3dy_planes2(Ngal_x+2,Ngal_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  !allocate(Qcrit(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  !allocate(Qcrit(N(1,2)+2,N(2,2),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate(Qcrit(Ngal_x+2,Ngal_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate( u1PLN(Nspec_x+2,Nspec_z,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u2PLN(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u3PLN(Nspec_x+2,Nspec_z,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u1PL_itpN(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u2PL_itpN(Nspec_x+2,Nspec_z,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u3PL_itpN(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( ppPLN(Nspec_x+2,Nspec_z,jgal(pgrid,1)-1:jgal(pgrid,2)+1))
  
  ! One field per proc (not per band)

  allocate( u1_itp( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )
  allocate( u2_itp( jlim(1,ugrid):jlim(2,ugrid), columns_num(myid) ) )
  allocate( u3_itp( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )

  allocate( Nu1_dy( jlim(1,ugrid)+1:jlim(2,ugrid)-1, columns_num(myid) ) )
  allocate( Nu2_dy( jlim(1,vgrid)+1:jlim(2,vgrid)-1, columns_num(myid) ) )
  allocate( Nu3_dy( jlim(1,ugrid)+1:jlim(2,ugrid)-1, columns_num(myid) ) )

  allocate( uv_f ( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )
  allocate( vv_c ( jlim(1,ugrid):jlim(2,ugrid), columns_num(myid) ) )
  allocate( wv_f ( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )

  allocate( DG( 3, jlim(1,pgrid):jlim(2,pgrid), columns_num(myid) ) )

  allocate( du1dy_columns( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )
  allocate( du2dy_columns( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )
  allocate( du3dy_columns( jlim(1,vgrid):jlim(2,vgrid), columns_num(myid) ) )


  u1_itp        = 0d0
  u2_itp        = 0d0
  u3_itp        = 0d0

  Nu1_dy        = 0d0
  Nu2_dy        = 0d0
  Nu3_dy        = 0d0

  uv_f          = 0d0
  wv_f          = 0d0
  vv_c          = 0d0

  DG            = 0d0

  du1dy_columns = 0d0
  du2dy_columns = 0d0
  du3dy_columns = 0d0



  !C! Matrix fir the laplacian of the pressure
  !Build (and decompose) matirx - only needs to be done once
  
  ! write(6,*) "Calling LU_buildP"

  call LU_buildP(jlim(1,pgrid),jlim(2,pgrid),myid,DG)

  ! write(6,*) "finished LU_buildP"
  allocate(   wx(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1  :jgal(vgrid,2)+1  ))

  allocate(spU (jlim(1,ugrid):jlim(2,ugrid),columns_num(myid)))
  allocate(spV (jlim(1,vgrid):jlim(2,vgrid),columns_num(myid)))
  allocate(spW (jlim(1,ugrid):jlim(2,ugrid),columns_num(myid)))
  allocate(spUV(jlim(1,ugrid):jlim(2,ugrid),columns_num(myid)))
  allocate(spP (jlim(1,pgrid):jlim(2,pgrid),columns_num(myid)))
  spU  = 0d0
  spV  = 0d0
  spW  = 0d0
  spUV = 0d0
  spP  = 0d0

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid == 0) then
    write(*,*) 'Reading in nonlinear interaction list'
    write(*,*) dirlist
    write(*,*)
  end if
  call nonlinRead
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  
end subroutine


subroutine ygrid
! TODO check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      YGRID  new   !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Creates the geometry in the y-direction
! y(j) = a(j-q) + b(j-q)**p

! y is the physical coordinate
! theta is the mapping coordinate (constant increments: dtheta)
!  Since (1) dtheta is constant and (2) the mapping between y and theta is analytical,
!  the second order of the centered finite difference is preserved.

! n   is the total number of point from wall to wall
! q   is the index of the centerline (n/2) (y(q) = 0)
! dyq is dymax/dymin (max stretching)
! p   polinomial exponent
! a and b are obtained such that dy(0)/dy(q)=dyq and y(0)=-1

  use declaration
  implicit none
  integer j,i
  real(8) aaa,bbb,qqq,dy0

  nn  = nyv+1 !v points (as collocated points)

  ! write(6,*) "nn", nn

  qqq = nn/2d0
  bbb = (dyq-1d0)/(qqq*(1d0-qqq**(ppp-1)+dyq*((1d0-qqq)*qqq**(ppp-1)-(1d0-qqq)**ppp)))
  aaa = (1d0-bbb*qqq**ppp)/qqq

  dtheta    = 1d0
  dthetai   = 1d0/dtheta
  ddthetavi = 1.0d0/dtheta

  dy0 = aaa*(1d0-qqq)+bbb*(1d0-qqq)**ppp+1d0

  allocate(yu    (  0:nn+1))
  allocate(dyu2i (3,0:nn+1))
  allocate(dthdyu(  0:nn+1))
  
  allocate(yv    (0:nn))
  allocate(dyv2i (3,0:nn))
  allocate(dthdyv(0:nn))


 !   !'Wall' grid (just smooth channel) ------------------------------------------------------------------------
    
    yu(0)=aaa*(-1d0+.5d0-qqq)+bbb*(-1d0+.5d0-qqq)**ppp  
    do j=0,nn                              ! Flow
      yv(j)=aaa*(j-qqq)+bbb*(j-qqq)**ppp
      yu(j+1)=aaa*(j+.5d0-qqq)+bbb*(j+.5d0-qqq)**ppp  
    end do

  
  ! Analytical expression for dtheta/dy in the lineal regions, and the polinomial one.
    dthdyu(0)=1d0/(aaa+ppp*bbb*(-1d0+.5d0-qqq)**(ppp-1))

    write(6,*) 'dthdyu(0) smooth', dthdyu(0)
    
    do j=0,nn
      dthdyv(j)=1d0/(aaa+ppp*bbb*(j-qqq)**(ppp-1))
      dthdyu(j+1)=1d0/(aaa+ppp*bbb*(j+.5d0-qqq)**(ppp-1))
    end do

! Second-order finite difference coefficient for 1. u(j-1), 2. u(j), 3. u(j+1)
  j=0
  dyu2i(1,j)= 0d0
  dyu2i(2,j)=-2d0/(yu(j+1)-yu(j))**2
  dyu2i(3,j)= 1d0/(yu(j+1)-yu(j))**2

  ! do i = 1,3
  !   write(6,*) "dyu2i", dyu2i(i,j)  ! checked correct
  ! end do
  
  dyv2i(1,j)= 0d0
  dyv2i(2,j)=-2d0/(yv(j+1)-yv(j))**2
  dyv2i(3,j)= 1d0/(yv(j+1)-yv(j))**2


  do j=1,nn-1
    dyu2i(1,j)=2d0/((yu(j)-yu(j-1))*(yu(j+1)-yu(j-1)))
    dyu2i(2,j)=-2d0/((yu(j+1)-yu(j))*(yu(j)-yu(j-1)))
    dyu2i(3,j)=2d0/((yu(j+1)-yu(j))*(yu(j+1)-yu(j-1)))
    
    dyv2i(1,j)=2d0/((yv(j)-yv(j-1))*(yv(j+1)-yv(j-1)))
    dyv2i(2,j)=-2d0/((yv(j+1)-yv(j))*(yv(j)-yv(j-1)))
    dyv2i(3,j)=2d0/((yv(j+1)-yv(j))*(yv(j+1)-yv(j-1)))
  end do

  ! do j= 1,nn-1
  !   write(6,*) dyu2i(1,j), dyu2i(2,j), dyu2i(3,j)
  !   ! write(6,*) dyv2i(1,j), dyv2i(2,j), dyv2i(3,j)
  ! end do

  j=nn
  dyu2i(1,j)=2d0/((yu(j)-yu(j-1))*(yu(j+1)-yu(j-1)))
  dyu2i(2,j)=-2d0/((yu(j+1)-yu(j))*(yu(j)-yu(j-1)))
  dyu2i(3,j)=2d0/((yu(j+1)-yu(j))*(yu(j+1)-yu(j-1)))


  
  dyv2i(1,j)= 1d0/(yv(j)-yv(j-1))**2
  dyv2i(2,j)=-2d0/(yv(j)-yv(j-1))**2
  dyv2i(3,j)= 0d0
  
  ! write(6,*) "j", j, "dyu2i=", dyu2i(1,nn), dyu2i(2,nn), dyu2i(3,nn), dyv2i(1,nn), dyv2i(2,nn), dyv2i(3,nn)

  j=nn+1
  dyu2i(1,j)= 1d0/(yu(j)-yu(j-1))**2
  dyu2i(2,j)=-2d0/(yu(j)-yu(j-1))**2
  dyu2i(3,j)= 0d0

  ! write(6,*) "j", j, "dyu2i=", dyu2i(1,nn+1), dyu2i(2,nn+1), dyu2i(3,nn+1)
  
  nn=nyv

  ! N(3,top_domain)=nn

  ! N(4,nband) = nyv + 1 ! extra poitn for u grid

  ! Ngal(3,3) = nyv
  ! Ngal(4,3) = N(4,3)

  ! write(6,*) "N:"
  ! do i = 1,4
  !   write(6,*) N(i,0:4)
  ! end do

  ! write(6,*) "Ngal:"
  ! do i = 1,4
  !   write(6,*) Ngal(i,0:4)
  ! end do


  ! write(6,*) "Ny:"
  ! do i = 1,3
  !   write(6,*) Ny(i,0:4)
  ! end do
end subroutine 


subroutine def_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    DEFINE K  NEW  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the wavenumbers k1 and their squares k2 = k1**2.
! They are defined for Fourier k1F, which gives an spectral accuracy,
!  which is in principle good. However, the highest frequencies are
!  not filtered and some wiggling may appear.
! In order to filter high frequencies a Differential wavenumber
!  is also defined, k1D, which reproduces the effective k of the
!  central difference scheme.

! ATTENTION: k1D is implemented in the code, but here, after being
!  defined, they are reasigned the same value as k1F.
!  Hence, when k1D is found in the code, it is really using k1F.

! Notes:
! alp = 2*pi/Lx (Lx being the streamwise periodicity of the box)
! bet = 2*pi/Lz (Lx being the spanwise   periodicity of the box)

  use declaration
  implicit none
  integer i,k
!  real(8) dzNgal,dzNgali

  allocate(k1F_x  (0:Nspec_x/2))
  allocate(k2F_x  (0:Nspec_x/2))
  allocate(k1F_z  (1:Nspec_z  ))
  allocate(k2F_z  (1:Nspec_z  ))

 
  !!!!!!   define differential operator eigenvalues !!!!!!
  do i = 0,Nspec_x/2
    k1F_x(i) =  im* alp*i
    k2F_x(i) = -   (alp*i)**2
  end do

  do k = 1,Nspec_z/2
    k1F_z(k) =  im* bet*(k-1)
    k2F_z(k) = -   (bet*(k-1))**2
  end do

  ! do i = 0,Nspec_x/2
  !   write(6,*) "k1F_z", k1F_z(i)
  ! end do 


!    k1F_z(N(2,1)/2+1) = 0d0 !These are now included
!    k2F_z(N(2,1)/2+1) = 0d0 !These are now included
  do k = Nspec_z/2+1,Nspec_z
    k1F_z(k) = -im* bet*(Nspec_z-k+1)
    k2F_z(k) = -   (bet*(Nspec_z-k+1))**2
  end do

  
  ! do i = N(2,3)/2+1,N(2,3)
  !   write(6,*) "k1F_z", k1F_z(i), i 
  ! end do 
  
  allocate(iLkup(-Nspec_x/2:Nspec_x/2), kLkup(-Nspec_z/2:Nspec_z/2), iNeg(-Nspec_x/2:Nspec_x/2))

  ! ! write(*,*) 'size(iLkup) = ', size(iLkup)
  ! write(*,*) 'size(kLkup) = ', size(kLkup)
  ! ! write(*,*) 'size(iNeg)  = ', size(iNeg)


  ! for nonlinear
  !! currently only works for when same discretisation for each band
  do i = 0,Nspec_x/2
    iLkup(i) =  i*2+1
    iNeg(i) = 1
    
  end do
  do i = -Nspec_x/2,-1
    iLkup(i) =  abs(i)*2+1
    iNeg(i) = -1
    !write(6,*) "ilkup", iLkup(i), "iNeg", iNeg(i)
  end do


  do k = 0,Nspec_z/2-1
    kLkup(k) =  k+1
    ! write(6,*) "kLkup", kLkup(k), k 
  end do

  kLkup(Nspec_z/2) = -Nspec_z/2 + 1 + Ngal_z
  do k = -Nspec_z/2,-1
    kLkup(k) = k + 1 + Ngal_z
    ! write(6,*) "kLkup", kLkup(k), k 
  end do
  
  
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    PROC LIMS COL NEW   !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine proc_lims_columns(myid)
  ! Calculates the different variables that define columns (spectral) for every proc 
  ! Creates columns_num, columns_i, columns_k, and jlim
  ! Creating list of columns with x and z coordinates:
  !   columns_i(column_index,myid), columns_k(column_index,myid) and columns_num(myid)
  ! The distribution for x and z is done in spectral space (Ring geometry)
  !    1. Bottom center (short pencils)     | 2 | 3 | 2 |
  !    2. Sides         (long  pencils)   y | 2 |   | 2 |
  !    3. Top    center (short pencils)     | 2 | 1 | 2 |
  !                                               z      

  use declaration
  implicit none

  include 'mpif.h'                         
  ! integer :: ierr

  integer :: iproc, accum, column, i, k, myid
  integer :: columns_num_total
  integer :: max_columns_num
  integer :: columns_num_total_proc_rem, columns_num_total_proc
  integer, allocatable :: columns_total_list_i(:) , columns_total_list_k(:)

  write(6,*) '>>> entered proc_lims_columns, myid =', myid
  ! Creating list of columns with x and z coordinates:
  !   columns_i(column_index,iband,myid), columns_k(column_index,iband,myid) and columns_num(iband,myid)

  ! NOTE 
  !   The columns corresponding to the last Fourier modes are also included.
  !   In the original code they are skipped in z. The * indicates places which need changes in order to skip some columns.

  ! Number of columns to be divided
  columns_num_total = (Nspec_x/2 + 1) * (Nspec_z) !*

  columns_num_total_proc      = columns_num_total /np
  columns_num_total_proc_rem  = columns_num_total  - columns_num_total_proc *np

  ! Array storing the number of columns per proc
  ! Ex. columns_num(myid) contains the number/amount of columns of proc numb "myid"
  allocate(columns_num(0:np-1))

  do iproc = 0,columns_num_total_proc_rem -1
    columns_num(iproc) = columns_num_total_proc  +1
  end do

  do iproc = columns_num_total_proc_rem, np-1
    columns_num(iproc) = columns_num_total_proc
  end do
  
  ! Array with the coordinates of all points
  ! All columns

  allocate(columns_total_list_i(columns_num_total))
  allocate(columns_total_list_k(columns_num_total))
  column = 0
  do   k =  1               ,Nspec_z/2   !do   k =  0               ,N(2,midband)/2-1  
    do i =  0               ,Nspec_x/2
      column = column +1
      columns_total_list_i(column) = i
      columns_total_list_k(column) = k

      !write(*,*) 'column=', column, ' i=', i, ' k=', k
    end do
  end do
  do   k = Nspec_z - Nspec_z/2 + 1, Nspec_z  !do   k = -N(2,midband)/2+1,-1
    do i =  0               ,Nspec_x/2
      column = column +1
      columns_total_list_i(column) = i
      columns_total_list_k(column) = k

      ! write(*,*) 'column=', column, ' i=', i, ' k=', k
    end do
  end do
  
  ! Creating the lists with column coordinates
  ! columns_i(ncolumn,iband,myid) has the x coordiante of column number ncolumn in band iband.
  max_columns_num = maxval(columns_num)
  allocate(columns_i(max_columns_num,0:np-1))
  allocate(columns_k(max_columns_num,0:np-1))


  accum = 0
  do iproc = 0,np-1
    ! write(6,*) "iproc", iproc
    ! write(6,*) "columns_num", columns_num(iproc), iproc
    do column = 1,columns_num(iproc)
      columns_i(column,iproc) = columns_total_list_i(column + accum)
      columns_k(column,iproc) = columns_total_list_k(column + accum)

      ! write(6,*) "columns_i =", columns_i(column, iproc),  &
      !     " columns_k =", columns_k(column, iproc), &
      !     " columns_long_list_i =", columns_total_list_i(column + accum), &
      !     " columns_long_list_k =", columns_total_list_k(column + accum)
    end do
    accum = accum + columns_num(iproc)
  end do
  
  deallocate(columns_total_list_i)
  deallocate(columns_total_list_k)

  allocate(jlim(2,3)) ! jlim(bottom(1)/top(2),grid,iband) 

  ! u and w: Include ghosts points
  jlim(1,ugrid) = nyu_LB
  jlim(2,ugrid) = nyu+1 

  ! write(*,*) 'jlim(1, ugrid)=', N(4,0), ' jlim(2, ugrid)=', N(4,3)+1

  ! v
  jlim(1,vgrid) = nyv_LB
  jlim(2,vgrid) = nyv+1 

  ! write(*,*) 'jlim(1, vgrid)=', N(3,0), ' jlim(2, vgrid)=', N(3,3)+1

  ! p: Ghost points not included
  jlim(1,pgrid) = nyu_LB+1
  jlim(2,pgrid) = nyu
  ! write(*,*) 'jlim(1, pgrid)=',N(4,0)+1, ' jlim(2, pgrid)=', N(4,3)
  ! write(*,*) jlim(1,ugrid), jlim(1,vgrid), 

  ! Used in FOU3D and stats
  ! It's a shift in z, used to align modes in different bands
  allocate(dk(max_columns_num,0:np-1))
  dk = 0
  do iproc = 0,np-1
        do column = 1,columns_num(iproc)
          if (columns_k(column,iproc) > Nspec_z/2) then
            dk(column,iproc) = Nspec_z-Ngal_z
            !write(6,*)'col_k=',columns_k(column,iproc), 'dk=', dk(column,iproc),"iproc", iproc
          end if
    end do
  end do

  !   do iproc = 0,np-1
  !       do column = 1,columns_num(iproc)
  !         if (columns_k(column,iproc) > N(2,3)/2) then
  !           write(6,*)'col_k=',columns_k(column,iproc), 'dk=', dk(column,iproc),"iproc", iproc
  !         end if
  !   end do
  ! end do

  
  allocate(dk_phys(max_columns_num,0:np-1))
  dk_phys = 0
  do iproc = 0,np-1
        do column = 1,columns_num(iproc)
          if (columns_k(column,iproc) > Nspec_z/2) then
            dk_phys(column,iproc) = Nspec_z-Nspec_z
          end if
    end do
  end do

  ! do iproc = 0,np-1
  !   do column = 1,columns_num(iproc)
  !     if (columns_k(column,iproc) > N(2,3)/2) then
  !       write(6,*)'col_k=',columns_k(column,iproc), 'dk_phys=', dk_phys(column,iproc),"iproc", iproc
  !     end if
  !   end do
  ! end do

  
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!    PROC LIMS planes NEW   !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine proc_lims_planes(myid)
  ! Calculates the different variables that define planes (phys) for every proc 
  ! Defines i-k-,j-gal, planelim and bandPL
  ! First computes roughly the proportional load that each proc should handle.
  !  Then it tries to share all the planes, allocating more procs for the finer planes.
  !  Once balanced, it establishes the limit planes (planelim) for the procs.

  use declaration
  use mpi
  implicit none

  integer :: myid, ierr
  ! integer :: iplanes
  integer :: iproc, i 
  integer(8) :: weight_tot, ideal_load
  integer :: jlow, jupp, proc_planes
  integer :: rem_planes

  ! new declarations 
  integer :: j, cut_plane, iplanes
  integer :: startproc, endproc
  real(8) :: cumulative_load, cumulative_load_1, cumulative_load_2
  integer, allocatable :: jmin_plane(:), jmax_plane(:)
  integer, allocatable :: nplanes(:)
  real(8), allocatable :: proc_load(:)

  allocate(planelim  (3,2,0:np-1))
  allocate(limPL_incw(3,2,0:np-1))
  allocate(limPL_excw(3,2,0:np-1))
  allocate(limPL_FFT(3,2,0:np-1))
  allocate(jmin_plane(0:np-1), jmax_plane(0:np-1))
  allocate(proc_load(0:np-1), nplanes(0:np-1))


  jlow = nyu_LB + 1
  jupp = nyu-1 ! jupp in get weights was N(4,3)-1 so weight allocation nneeded to be weight 
                    ! (jlow:jupp+1) but here its jsut jlow to j upp b it included the +1 already

  proc_load = 0.0d0

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DNS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! !-------------------------------------------------------------------------------------
  ! ! Simple equal-plane partition in j: give each proc roughly same no. planes for DNS
  ! !-------------------------------------------------------------------------------------
  
  ! ! total number of active planes in j (without walls)
  ! iplanes = jupp - jlow + 1

  ! ! base number of planes per proc and remainder
  ! proc_planes = iplanes / np
  ! rem_planes  = mod(iplanes, np)

  ! j = jlow
  ! do iproc = 0, np-1

  !    ! first rem_planes procs get one extra plane
  !    if (iproc < rem_planes) then
  !       nplanes(iproc) = proc_planes + 1
  !    else
  !       nplanes(iproc) = proc_planes
  !    end if

  !    jmin_plane(iproc) = j
  !    jmax_plane(iproc) = j + nplanes(iproc) - 1

  !    j = jmax_plane(iproc) + 1

  !    ! if you still want a "load" number, just set it ~ nplanes
  !    proc_load(iproc) = nplanes(iproc)

  ! end do

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!  setting global variables !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! do iproc = 0, np-1

  !   ! default: all grids use full j-range for this proc
  !   planelim(ugrid,1,iproc) = jmin_plane(iproc)
  !   planelim(ugrid,2,iproc) = jmax_plane(iproc)

  !   planelim(vgrid,1,iproc) = jmin_plane(iproc)
  !   planelim(vgrid,2,iproc) = jmax_plane(iproc)

  !   planelim(pgrid,1,iproc) = jmin_plane(iproc)
  !   planelim(pgrid,2,iproc) = jmax_plane(iproc)

  !   ! special case: first proc (shift p-grid start)
  !   if (iproc == 0) then
  !     planelim(pgrid,1,iproc) = jmin_plane(iproc) + 1

  !   ! special case: last proc (shrink v and p by 1 at top)
  !   else if (iproc == np-1) then
  !     planelim(vgrid,2,iproc) = jmax_plane(iproc) - 1
  !     planelim(pgrid,2,iproc) = jmax_plane(iproc) - 1
  !   end if

  !   nplanes(iproc) = planelim(ugrid,2,iproc) - planelim(ugrid,1,iproc) + 1

  !   ! bandPL(iproc) = nband ! keeping this for now bc otherwise code will break but setting at a const

  ! end do


  ! ! limPL_incw: like old code – same as planelim but extend first/last proc by 1 plane
  ! ! limPL_incw is like planelim but including first and last points that were previously removed.
  ! !  It is only used in planes_to_modes_UVP, modes_to_planes_UVP, record_out and stats
  ! limPL_excw = planelim

  ! limPL_incw = planelim
  ! limPL_incw(:,1,0   ) = planelim(:,1,0   ) - 1
  ! limPL_incw(:,2,np-1) = planelim(:,2,np-1) + 1

  ! igal = Ngal_x+2
  ! kgal = Ngal_z


  ! jgal(ugrid,1) = limPL_excw(ugrid,1,myid)
  ! jgal(ugrid,2) = limPL_excw(ugrid,2,myid)
  ! jgal(vgrid,1) = limPL_excw(vgrid,1,myid)
  ! jgal(vgrid,2) = limPL_excw(vgrid,2,myid)
  ! jgal(pgrid,1) = limPL_excw(pgrid,1,myid)
  ! jgal(pgrid,2) = limPL_excw(pgrid,2,myid)

  ! ! write(6,*) "jgal(ugrid,1)", jgal(ugrid,1), "jgal(ugrid,2)", jgal(ugrid,2), myid


  ! if (myid == 0) then
  !   do iproc = 0, np-1
  !     write(6,*) "iproc", iproc, &
  !               "nplanes", nplanes(iproc), &
  !               "proc_load", proc_load(iproc)
  !   end do
  ! end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RNL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Rank 0 computes plane ownership for all ranks (indexed by iproc)
  ! then sends this info to the other procs who can access it 


  if (myid ==0) then 

    weight_tot = sum(weight)  
    ideal_load = weight_tot / np

    write(6,*) "weight_tot",weight_tot, "ideal_load", ideal_load
    ! write(6,*) "weight(151)",weight(151), "weight(1)", weight(1)

    jmin_plane(0) = jlow

    !initialising values 
    cumulative_load = 0.0d0
    iproc = 0
    j = jlow
    ! proc_load = 0
    ! nplanes   = 0

    do while (iproc <= np-1)
      if (j<= jupp) then 
        cumulative_load = cumulative_load + weight(j)
        if (cumulative_load >= ideal_load) then
          ! ! default: cut at j
          ! cut_plane = j
          ! proc_load(iproc) = cumulative_load

          cut_plane = j
          proc_load(iproc) = cumulative_load

          if (j > jmin_plane(iproc)) then
            cumulative_load_2 = cumulative_load - weight(j)     ! load up to j-1
            if (abs(cumulative_load_2 - ideal_load) < abs(proc_load(iproc) - ideal_load)) then
              cut_plane = j-1
              proc_load(iproc) = cumulative_load_2
            end if
          end if


          !define max j fort hat proc and the no. planes 
          jmax_plane(iproc) = cut_plane
          nplanes(iproc) = jmax_plane(iproc) - jmin_plane(iproc) +1
          
          ! if we get to last proc before all j planes are taken, fill the last proc with all remaining planes 
          if (iproc == np-1) then
            jmax_plane(iproc) = jupp
            proc_load(iproc)  = real(weight_tot,8) - sum(proc_load(0:iproc-1))
            exit
          else
            iproc = iproc + 1
            jmin_plane(iproc) = cut_plane + 1
          end if

          cumulative_load = 0.0d0
          j = cut_plane + 1
          endproc = iproc   ! <-- the last filled proc index
        else
          j = j + 1
        end if
      else 
        !if we get to top proc before all planes are filled... 
        endproc = iproc   ! <-- the last filled proc index
        jmax_plane(iproc) = jupp
        proc_load(iproc)  = real(weight_tot,8) - sum(proc_load(0:iproc-1))
        nplanes(iproc) = jmax_plane(iproc) - jmin_plane(iproc) +1
        do i = iproc+1, np-1
          jmin_plane(i) = jupp + 1
          jmax_plane(i) = jupp
          proc_load(i)  = 0.0d0
        end do
        exit
      end if 

    end do

    !!!!!!!!!!!!!!!!!!!!  redistributing planes w nothing.. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! takes one from each of the last n planes and shifts the planes and the weights 
    ! this frees up n planes where one plane will be given to the remaining unloaded procs 
    if (endproc < np-1) then 
      startproc = endproc - (np-1 - endproc) +1
      proc_load(startproc) = proc_load(startproc)-weight(jmax_plane(startproc))
      jmax_plane(startproc) = jmax_plane(startproc) -1
      
      do i = startproc+1, endproc
        proc_load(i) = 0 
        nplanes(i) = nplanes(i)-1   
        jmin_plane(i) = jmax_plane(i-1)+1
        jmax_plane(i) = jmin_plane(i)+ nplanes(i) -1 !this minus 1 bc if start j = 110, then 3 planes is acc 110,111,112
        do j =  jmin_plane(i), jmax_plane(i)
          proc_load(i) = proc_load(i) + weight(j)
        end do 
      end do 

      ! treatmeant of planes w. 0 load

      do i = endproc+1, np-1  
        jmin_plane(i) = jmax_plane(i-1) +1
        jmax_plane(i) = jmin_plane(i)
        proc_load(i) = weight(jmax_plane(i))
      end do 

    end if 

  end if 

  call MPI_BCAST(jmin_plane, np, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(jmax_plane, np, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(proc_load,  np, MPI_DOUBLE_PRECISION,   0, MPI_COMM_WORLD, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!           RNL SPECIFIC         !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do iproc = 0, np-1

    ! default: all grids use full j-range 
    planelim(ugrid,1,iproc) = jmin_plane(iproc)
    planelim(ugrid,2,iproc) = jmax_plane(iproc)

    planelim(vgrid,1,iproc) = jmin_plane(iproc)
    planelim(vgrid,2,iproc) = jmax_plane(iproc)

    planelim(pgrid,1,iproc) = jmin_plane(iproc)
    planelim(pgrid,2,iproc) = jmax_plane(iproc)

    ! special case: first proc (shift p-grid start)
    if (iproc == 0) then
      planelim(pgrid,1,iproc) = jmin_plane(iproc) + 1

    end if

    ! special case: shift final u-grid by 1
    if (iproc == np-1) then
      planelim(ugrid,2,iproc) = jmax_plane(iproc) + 1
    end if

    nplanes(iproc) = planelim(ugrid,2,iproc) - planelim(ugrid,1,iproc) + 1

  end do


  ! limPL_incw: like old code – same as planelim but extend first/last proc by 1 plane
  ! limPL_incw is like planelim but including first and last points that were previously removed.
  !  It is only used in planes_to_modes_UVP, modes_to_planes_UVP, record_out and stats
  limPL_excw = planelim

  limPL_incw = planelim
  limPL_incw(:,1,0   ) = planelim(:,1,0   ) - 1
  limPL_incw(:,2,np-1) = planelim(:,2,np-1) + 1

  igal = Ngal_x+2
  kgal = Ngal_z


  jgal(ugrid,1) = limPL_excw(ugrid,1,myid)
  jgal(ugrid,2) = limPL_excw(ugrid,2,myid)
  jgal(vgrid,1) = limPL_excw(vgrid,1,myid)
  jgal(vgrid,2) = limPL_excw(vgrid,2,myid)
  jgal(pgrid,1) = limPL_excw(pgrid,1,myid)
  jgal(pgrid,2) = limPL_excw(pgrid,2,myid)

  !write(6,*) "jgal(ugrid,1)", jgal(ugrid,1), "jgal(ugrid,2)", jgal(ugrid,2), myid


  if (myid == 0) then
    do iproc = 0, np-1
      write(6,*) "iproc", iproc, &
                "nplanes", nplanes(iproc), &
                "proc_load", proc_load(iproc)
    end do
  end if

  deallocate(jmin_plane, jmax_plane, nplanes, proc_load)



end subroutine

subroutine getini(u1,u2,u3,p,div,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     GETINI  NEW   !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Grab the initial condition and initialize some variables

  use declaration
  ! use littleharsh_mod
  ! use transpose

  implicit none

  include 'mpif.h'                                  ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables
  !complex(8) :: u1(jlim(1,ugrid):,:), u2(jlim(1,vgrid):,:), u3(jlim(1,ugrid):,:), p(:,:)
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
  ! complex(8) :: div(:,:)
  complex(8) :: div( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )

  ! write(*,*) "call shape(u1)=", shape(u1), "lb=", lbound(u1), "ub=", ubound(u1)


  if (myid==0) then
    write(*,*) 'Launching...'
  end if

  if (flag_init==1) then       ! Initial conditions borrowed from another 
    if (myid==0) then
      write(*,*) 'starting from multi-block, flat channel?'
      write(*,*) 'start file: ',trim(dirin)//trim(fnameimb)
    end if
    call mblock_ini(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
    if (myid==0) then
      ! call flowrateIm(Qx,u1(nyu_LB,1))        !check why flowrate used midband ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      call flowrateIm(Qx,u1(:,1))
      !write(6,*) "flowrateIm", Qx
      if (flag_ctpress/=1) then
        QxT = Qx
      end if
    end if
    iter0 = 0
    mpgz  = 0d0
    t     = 0d0
  else if (flag_init==2) then  ! Continuing simulation
    if (myid==0) then
      write(*,*) 'continuing simulation'
      write(*,*) 'start file:',fnameimb
    end if
    call mblock_ini(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
    if (myid==0) then
      !call flowrateIm(Qx,u1(nyu_LB,1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      call flowrateIm(Qx,u1(:,1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      ! TODO change this condition for 'flag_ctpress == 0'
      if (flag_ctpress/=1) then  ! Constant flow rate (flag_ctpress == 1, constant pressure gradient)
        QxT = Qx
      end if
    end if
    mpgz = 0d0
  else if (flag_init==3) then  ! Parabolic profile
    if (myid==0) then
      write(*,*) 'parabolic profile'
    end if
    !call mblock_ini_parabolic_profile(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
    if (myid==0) then
      call flowrateIm(Qx,u1(:,1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      ! TODO change this condition for 'flag_ctpress == 0'
      if (flag_ctpress/=1) then  ! Constant flow rate (flag_ctpress == 1, constant pressure gradient)
        QxT = Qx
      end if
    end if
    iter0 = 0
    mpgz  = 0d0
    t     = 0d0
  else 
    write(*,*) 'INITIAL CONDITIONS NOT IMPLEMENTED YET'
    call MPI_FINALIZE(ierr)
    stop
  end if

  ! write(6,*) "check1 ", myid

  ! Broadcast the time step and time.
  ! If it's initialized from a previous simulation this value is already known, otherwise it's set to 0
  call MPI_BCAST(iter0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(t    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
  iter   = iter0
  !  iter0=iter-nstat
  iwrite = iter

  ! 'Probably' this is used to initialize the divergence in the case of a new simulation
  ! write(6,*) "call divergence ", myid
  call divergence(div,u1,u2,u3,myid)

  ! write(6,*) "call init stats", myid
  call init_stats(myid)
  ! write(6,*) "finished init_stats", myid
  
  ! write(6,*) "call init_sl stats", myid
  ! call init_sl_stats(myid)

  if (myid==0) then
    write(*,*) 'Lx    ',Lx
    write(*,*) 'Ly    ',Ly
    write(*,*) 'Lz    ',Lz
    write(*,*) 'Nx    ',Nspec_x
    write(*,*) 'Nz    ',Nspec_z
    write(*,*) 'Nyv   ',Nyv
    write(*,*) 'Nyu   ',nyu
    write(*,*) 'Ngalx ',Ngal_x
    write(*,*) 'Ngalz ',Ngal_z
    write(*,*) ''
    write(*,*) 'dymin ',yu(1)-yu(0)
    write(*,*) 'dymax ',yu((nyu+1)/2+1)-yu((nyu+1)/2)
    write(*,*) ''
    write(*,*) 'Re ',Re
    write(*,*) 't  ',t
  end if

end subroutine


subroutine mblock_ini(u1,u2,u3,p,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!   MBLOCK INI   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  ! use transpose
  implicit none

  include 'mpif.h'                                  ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables

  integer j,i,k,column
  ! real(8) sigmaz,z,sigmax,x,fact
  ! complex(8), intent(inout) :: u1(:,:), u2(:,:), u3(:,:), p(:,:)
  complex(8) :: u1 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: u2 ( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) )
  complex(8) :: u3 ( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
  complex(8) :: p  ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )

  u1PL = 0d0
  u2PL = 0d0
  u3PL = 0d0
  ppPL = 0d0

  ! write(*,*) "shape(u1)=", shape(u1), "lb=", lbound(u1), "ub=", ubound(u1)



  u1 = 0d0
  u2 = 0d0
  u3 = 0d0
  p = 0d0

  write(6,*) " Calling read in"
  call read_in(myid)
  
  ! write(6,*) " start planes to modes"
  call planes_to_modes_UVP(u1,u1PL,2,nyu,nyu_LB,myid,status,ierr)
  call planes_to_modes_UVP(u2,u2PL,1,nyv,nyv_LB,myid,status,ierr)
  call planes_to_modes_UVP(u3,u3PL,2,nyu,nyu_LB,myid,status,ierr)
  call planes_to_modes_UVP(p ,ppPL,3,nyp,nyp_LB,myid,status,ierr)
  ! write(6,*) " finished pplanes to modes"

  if (recv_flg ==1) then
    call reconstruct_u3(u3,u1,u2,myid)
  end if 

  ! if(myid ==0) then 
  !   do j= 1, 22
  !     write(6,*) p%f(j,1)
  !   end do
  ! end if 

end subroutine


subroutine mblock_ini_parabolic_profile(u1,u2,u3,p,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! MBLOCK INI PARABOLIC PROFILE NEW !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 

  use declaration
  ! use transpose
  implicit none

  include 'mpif.h'                                  ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables
  ! type(cfield) u1,u2,u3,p
  complex(8), intent(in) :: u1(jlim(1,ugrid):,:), u2(jlim(1,vgrid):,:), u3(jlim(1,ugrid):,:)
  complex(8), intent(in) :: p(:,:)
  integer nx,nz
  integer j,ju1,ju2,jv1,jv2,jp1,jp2, iband
  integer, allocatable:: nxxu(:),nzzu(:),nxxv(:),nzzv(:)
  real(8), allocatable:: buffSR(:,:)

  integer column,i,k

  u1PL = 0d0
  u2PL = 0d0
  u3PL = 0d0
  ppPL = 0d0

  ! uncomment, remove

  ! u1%f = 0d0
  ! u2%f = 0d0
  ! u3%f = 0d0
  ! p %f = 0d0


  ! ! Function to create parabolic profile
  ! allocate(nxxu(nyu_LB:nyu+1),nzzu(nyu_LB:nyu+1))
  ! allocate(nxxv(nyv_LB:nyv+1),nzzv(nyv_LB:nyv+1))
  ! nxxu(nyu_LB)=N(1,1)+2
  ! nzzu(nyu_LB)=N(2,1)
  ! nxxv(nyv_LB)=N(1,1)+2
  ! nzzv(nyv_LB)=N(2,1)
  ! do iband=1,3
  !   do j=N(4,iband-1)+1,N(4,iband)
  !     nxxu(j)=N(1,iband)+2
  !     nzzu(j)=N(2,iband)
  !   end do
  !   do j=N(3,iband-1)+1,N(3,iband)
  !     nxxv(j)=N(1,iband)+2
  !     nzzv(j)=N(2,iband)
  !   end do
  ! end do
  ! nxxv(nyv+1)=Nspec_x+2
  ! nzzv(nyv+1)=Nspec_z
  ! nxxu(nyu+1)=Nspec_x+2
  ! nzzu(nyu+1)=Nspec_z
  ! if (myid==0) then
  !  ju1=jgal(2,1)-1
  !  ju2=jgal(2,2)
  !  ju1=max(ju1,nyu_LB)
  ! else
  !  ju1=jgal(2,1)
  !  ju2=jgal(2,2)
  !  if (jgal(2,2)==nyu) then
  !    ju2=jgal(2,2)+1
  !  end if
  !  ju1=max(ju1,nyu_LB)
  !  ju2=min(ju2,nyu+1)
  ! end if
  ! if (myid==0) then
  !  jv1=jgal(1,1)-1
  !  jv2=jgal(1,2)
  !  jv1=max(jv1,nyv_LB)
  ! else
  !  jv1=jgal(1,1)
  !  jv2=jgal(1,2)
  !  if (jgal(1,2)==nyv) then
  !    jv2=jgal(1,2)+1
  !  end if
  !  jv1=max(jv1,nyv_LB)
  !  jv2=min(jv2,nyv+1)
  ! end if
  ! if (myid==0) then
  !  jp1=jgal(3,1)-1
  !  jp2=jgal(3,2)
  !  jp1=max(jp1,nyu_LB+1)
  ! else
  !  jp1=jgal(3,1)
  !  jp2=jgal(3,2)
  !  if (jgal(3,2)==nyu-1) then
  !    jp2=jgal(3,2)+1
  !  end if
  !  jp1=max(jp1,nyu_LB+1)
  !  jp2=min(jp2,nyu+1-1)
  ! end if
  ! !!!!!!!!!!!!!!    u1    !!!!!!!!!!!!!!
  ! do j=ju1,ju2
  !   nx=nxxu(j)
  !   nz=nzzu(j)
  !   allocate(buffSR(nx,nz))
  !   buffSR     =0d0
  !   buffSR(1,1)=(0.5d0)*Re*mpgx*(yu(j)**2-1)                    ! Parabolic profile: 1/2*Re*dp/dx*(y^2-1) (Reynolds bulk)
  !   call buff_to_u(u1PL(1,1,j),buffSR,nx,nz,igal,kgal)
  !   deallocate(buffSR)
  ! end do
  ! !!!!!!!!!!!!!!    u2    !!!!!!!!!!!!!!
  ! do j=jv1,jv2
  !   nx=nxxv(j)
  !   nz=nzzv(j)
  !   allocate(buffSR(nx,nz))
  !   buffSR     =0
  !   buffSR(1,1)=1e-6                                          ! Some noise
  !   call buff_to_u(u2PL(1,1,j),buffSR,nx,nz,igal,kgal)
  !   deallocate(buffSR)
  ! end do
  ! !!!!!!!!!!!!!!    u3    !!!!!!!!!!!!!!
  ! do j=ju1,ju2
  !   nx=nxxu(j)
  !   nz=nzzu(j)
  !   allocate(buffSR(nx,nz))
  !   buffSR     =0d0
  !   buffSR(1,1)=1e-5                                          ! Some noise 
  !   call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
  !   deallocate(buffSR)
  ! end do
  ! !!!!!!!!!!!!!!    p     !!!!!!!!!!!!!!
  ! do j=jp1,jp2
  !   nx=nxxu(j)
  !   nz=nzzu(j)
  !   allocate(buffSR(nx,nz))
  !   buffSR     =0d0
  !   buffSR(1,1)=1d0                                           ! Some noise 
  !   call buff_to_u(ppPL(1,1,j),buffSR,nx,nz,igal,kgal)
  !   deallocate(buffSR)
  ! end do
  ! deallocate(nxxu,nzzu,nxxv,nzzv)

  ! call planes_to_modes_UVP(u1,u1PL,ugrid,nyu,nyu_LB,myid,status,ierr)
  ! call planes_to_modes_UVP(u2,u2PL,vgrid,nyv,nyv_LB,myid,status,ierr)
  ! call planes_to_modes_UVP(u3,u3PL,ugrid,nyu,nyu_LB,myid,status,ierr)
  ! call planes_to_modes_UVP(p ,ppPL,pgrid,nyp,nyp_LB,myid,status,ierr)

end subroutine


subroutine read_in(myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    READ IN  NEW   !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reads the file with the inital condition

! The variables u1PL, u2PL, u3PL and ppPL are sent to every proc
! The procs only stores the planes they have to compute
!  (jgal is the local name of planelim for j)

  use declaration
  ! use FOU3D_mod
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid
  integer nx,nz,nxin,nzin,nband2
  integer j,jin,iproc,dummI,ju1,ju2,jv1,jv2,jp1,jp2,i
  real(8) dummRe,Re2,alp2,bet2
  real(8), allocatable:: buffSR(:,:),dumm_y(:)
  integer, allocatable:: dummint(:),N2(:,:)
  integer, allocatable:: nxxu(:),nzzu(:),nxxv(:),nzzv(:),nxxp(:),nzzp(:)

  filout   = fnameimb(3:index(fnameimb,' ')-1)

  if (myid==0) then
    ! Read the file with the initial conditions
    fnameimb = trim(dirin)//'/u1'//filout
    open(10,file=fnameimb,form='unformatted')
    allocate(dummint(88))
    read(10) t,Re2,alp2,bet2,dummRe,nband2,iter0,dummint
!    if (flag_ctpress==1) then
!      mpgx=dummRe
!    end if
    deallocate(dummint)

    ! Checks that the size of the variables in the file matches with the setup of actual problem
    ! Checks nband, alp, bet, N
    if (nband2/=3) then
      write(*,*) ''
      write(*,*) 'wrong startfile?'
      write(*,*) 'nband,nband_old=',nband,nband2
      stop
    else if (alp2/=alp .or. bet2/=bet) then
      write(*,*) ''
      write(*,*) 'wrong startfile?'
      write(*,*) 'alp,alp_old=',alp,alp2
      write(*,*) 'bet,bet_old=',bet,bet2
      write(*,*) ''
    end if
    write(*,*) 'Re_old,Re',Re2,Re
    allocate(N2(4,0:nband2+1))
    read(10) N2
    close(10)


    if(myid ==0) then 
      write(6,*) "N2"
      do i= 1,4
        write(6,*) N2(i, 0:4)
      end do 
    end if 

    dummI=0

    if (Nspec_x/=N2(1,3)) then
      dummI=2
    end if


    if (nyv_LB/=N2(3,0)) then
      dummI=2
    else if (nyv/=N2(3,3)) then
      dummI=2
    end if

    if (dummI==1) then
      write(*,*) ''
      write(*,*) 'wrong startfile?',dummI
      write(*,*) 'N_new='
      write(*,*) Nspec_x
      write(*,*) Nspec_z
      write(*,*) nyv
      write(*,*) nyu
      write(*,*) ''
      write(*,*) 'N_old='
      write(*,*) N2(1,0:4)
      write(*,*) N2(2,0:4)
      write(*,*) N2(3,0:4)
      write(*,*) N2(4,0:4)
      stop
    else if (dummI==2) then
      write(*,*) ''
      write(*,*) 'WARNING startfile size',dummI
      write(*,*) 'N_new='
      write(*,*) Nspec_x
      write(*,*) Nspec_z
      write(*,*) nyv
      write(*,*) nyu
      write(*,*) ''
      write(*,*) 'N_old='
      write(*,*) N2(1,0:4)
      write(*,*) N2(2,0:4)
      write(*,*) N2(3,0:4)
      write(*,*) N2(4,0:4)
      write(*,*) ''
    end if

    ! Despite possible unmatches, it sends N to the procs
    do iproc=1,np-1
      call MPI_SEND(N2,4*(nband2+2),MPI_INTEGER,iproc,121*iproc,MPI_COMM_WORLD,ierr)
    end do

    ! nxx(j) and nzz(j) store the number of points at a plane j
    allocate(nxxu(N2(4,0):N2(4,3)+1),nzzu(N2(4,0):N2(4,3)+1))
    allocate(nxxv(N2(3,0):N2(3,3)+1),nzzv(N2(3,0):N2(3,3)+1))
    allocate(nxxp(N2(4,0)+1:N2(4,3)+1-1),nzzp(N2(4,0)+1:N2(4,3)+1-1))
    nxxu(N2(4,0))=N2(1,3)+2
    nzzu(N2(4,0))=N2(2,3)
    nxxv(N2(3,0))=N2(1,3)+2
    nzzv(N2(3,0))=N2(2,3)
    nxxp(N2(4,0)+1)=N2(1,3)+2
    nzzp(N2(4,0)+1)=N2(2,3)

    do j=N2(4,0)+1,N2(4,3)
      nxxu(j)=N2(1,3)+2
      nzzu(j)=N2(2,3)
    end do
    do j=N2(3,0)+1,N2(3,3)
      nxxv(j)=N2(1,3)+2
      nzzv(j)=N2(2,3)
    end do
    do j= N2(4,0)+1,N2(4,3)-1
      nxxp(j)=N2(1,3)+2
      nzzp(j)=N2(2,3)
    end do

    ! nxxp(N2(4,1))=nxxp(N2(4,1)+1)
    ! nzzp(N2(4,1))=nzzp(N2(4,1)+1)
    ! nxxp(N2(4,2))=nxxp(N2(4,2)-1)
    ! nzzp(N2(4,2))=nzzp(N2(4,2)-1)
    ! nxxp(N2(4,2)+1)=nxxp(N2(4,2)-1)
    ! nzzp(N2(4,2)+1)=nzzp(N2(4,2)-1)
    ! end do


    nxxu(N2(4,3)+1)=N2(1,3)+2
    nzzu(N2(4,3)+1)=N2(2,3)
    nxxv(N2(3,3)+1)=N2(1,3)+2
    nzzv(N2(3,3)+1)=N2(2,3)
    nxxp(N2(4,3)+1-1)=N2(1,3)+2
    nzzp(N2(4,3)+1-1)=N2(2,3)
    !filout=fnameimb(3:index(fnameimb,' ')-1)
    !filout=fnameimb(10:index(fnameimb,' ')-1) !22 as now in subfolder

    write(*,*) 'getting'

    ! Reads u1, u2, u3 and p,  and send them to the procs

    !!!!!!!!!!!!!!    u1    !!!!!!!!!!!!!!
    fnameimb = trim(dirin)//'/u1'//filout
    write(*,*) 'u1 from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    ju1=jgal(2,1)-1
    ju2=jgal(2,2)
    ju1=max(ju1,N2(4,0))
    do j=N2(4,0),ju1-1
      read(10)
    end do
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(u1PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      ju1=planelim(2,1,iproc)
      ju2=planelim(2,2,iproc)
      if (planelim(2,2,iproc)==nyu) then
        ju2=planelim(2,2,iproc)+1
      end if
      ju1=max(ju1,N2(4,0))
      ju2=min(ju2,N2(4,3)+1)
      do j=ju1,ju2
        nx=nxxu(j)
        nz=nzzu(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,123*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)

    !!!!!!!!!!!!!!    u2    !!!!!!!!!!!!!!
    fnameimb = trim(dirin)//'/u2'//filout
    write(*,*) 'u2 from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    jv1=jgal(1,1)-1
    jv2=jgal(1,2)
    jv1=max(jv1,N2(3,0))
    do j=N2(3,0),jv1-1
      read(10)
    end do
    do j=jv1,jv2
      nx=nxxv(j)
      nz=nzzv(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(u2PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      jv1=planelim(1,1,iproc)
      jv2=planelim(1,2,iproc)
      if (planelim(1,2,iproc)==nyv.and.iproc==np-1) then
        jv2=planelim(1,2,iproc)+1
      end if
      jv1=max(jv1,N2(3,0))
      jv2=min(jv2,N2(3,3)+1)
!jv2=min(jv2,N2(3,nband))
      do j=jv1,jv2
        nx=nxxv(j)
        nz=nzzv(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,124*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)

    ! !!!!!!!!!!!!!!    u3    !!!!!!!!!!!!!!
    ! fnameimb = trim(dirin)//'/u3'//filout
    ! write(*,*) 'u3 from file ',trim(fnameimb),','
    ! open(10,file=fnameimb,form='unformatted')
    ! read(10)
    ! read(10)
    ! read(10)
    ! ju1=jgal(2,1)-1
    ! ju2=jgal(2,2)
    ! ju1=max(ju1,N2(4,0))
    ! do j=N2(4,0),ju1-1
    !   read(10)
    ! end do
    ! do j=ju1,ju2
    !   nx=nxxu(j)
    !   nz=nzzu(j)
    !   allocate(buffSR(nx,nz))
    !   read(10) jin,dummI,nxin,nzin,dummRe,buffSR
    !   call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
    !   deallocate(buffSR)
    !   if (nx/=nxin .or. nz/=nzin) then
    !     write(*,*) 'WARNING!: unexpected size of plane',j
    !   end if
    ! end do
    ! do iproc=1,np-1
    !   ju1=planelim(2,1,iproc)
    !   ju2=planelim(2,2,iproc)
    !   if (planelim(2,2,iproc)==nyu) then
    !     ju2=planelim(2,2,iproc)+1
    !   end if
    !   ju1=max(ju1,N2(4,0))
    !   ju2=min(ju2,N2(4,3)+1)
    !   do j=ju1,ju2
    !     nx=nxxu(j)
    !     nz=nzzu(j)
    !     allocate(buffSR(nx,nz))
    !     read(10) jin,dummI,nxin,nzin,dummRe,buffSR
    !     call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,ierr)
    !     deallocate(buffSR)
    !     if (nx/=nxin .or. nz/=nzin) then
    !       write(*,*) 'WARNING!: unexpected size of plane',j
    !     end if
    !   end do
    ! end do
    ! close(10)
    
    !!!!!!!!!!!!!!    u3  NEW  !!!!!!!!!!!!!!
    fnameimb = trim(dirin)//'/u3'//filout
    write(*,*) 'u3 from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    ju1=jgal(2,1)-1
    ju2=jgal(2,2)
    ju1=max(ju1,N2(4,0))
    do j=N2(4,0),ju1-1
      read(10)
    end do
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      buffSR(:,:) = 0
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR(:,1)

      if (j == ju1) then
        recv_flg = merge(1, 0, nzin == 1)
        call MPI_BCAST(recv_flg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if

      call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    if (nzin == 1)  then
      do iproc=1,np-1
        ju1=planelim(2,1,iproc)
        ju2=planelim(2,2,iproc)
        if (planelim(2,2,iproc)==nyu) then
          ju2=planelim(2,2,iproc)+1
        end if
        ju1=max(ju1,N2(4,0))
        ju2=min(ju2,N2(4,3)+1)
        do j=ju1,ju2
          nx=nxxu(j)
          nz=nzzu(j)
          allocate(buffSR(nx,nz))
          buffSR(:,:) = 0
          read(10) jin,dummI,nxin,nzin,dummRe,buffSR(:,1)
          ! call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,ierr)
          call MPI_SEND(buffSR(:,1), nx, MPI_REAL8, iproc, 125*iproc, MPI_COMM_WORLD, ierr)
          deallocate(buffSR)
          if (nx/=nxin .or. nzin/=1) then
            write(*,*) 'WARNING!: unexpected size of plane',j
          end if
        end do
      end do
    else 
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
      do iproc=1,np-1
        ju1=planelim(2,1,iproc)
        ju2=planelim(2,2,iproc)
        if (planelim(2,2,iproc)==nyu) then
          ju2=planelim(2,2,iproc)+1
        end if
        ju1=max(ju1,N2(4,0))
        ju2=min(ju2,N2(4,3)+1)
        do j=ju1,ju2
          nx=nxxu(j)
          nz=nzzu(j)
          allocate(buffSR(nx,nz))
          read(10) jin,dummI,nxin,nzin,dummRe,buffSR
          call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,ierr)
          deallocate(buffSR)
          if (nx/=nxin .or. nz/=nzin) then
            write(*,*) 'WARNING!: unexpected size of plane',j
          end if
        end do
      end do
    end if 
    close(10)
    


    !!!!!!!!!!!!!!    p     !!!!!!!!!!!!!!

    fnameimb = trim(dirin)//'/p'//filout
    write(*,*) 'and p from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    jp1=jgal(3,1)-1
    jp2=jgal(3,2)

    
    !if(myid ==0) then 
    ! write (6,*) "jp1", jp1, "jp2", jp2, "jgal(3,1)", jgal(3,1), myid 
    !end if 

    jp1=max(jp1,N2(4,0)+1)
    do j=N2(4,0)+1,jp1-1
      read(10)
    end do
    do j=jp1,jp2
      nx=nxxp(j)
      nz=nzzp(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(ppPL(1,1,j),buffSR,nx,nz,igal,kgal)
      !write(6,*) "ppPL", ppPL(1,1,j), j 
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      jp1=planelim(3,1,iproc)
      jp2=planelim(3,2,iproc)
      if (planelim(3,2,iproc)==nyu-1.and.iproc==np-1) then
        jp2=planelim(3,2,iproc)+1
      end if
      jp1=max(jp1,N2(4,0)+1)
      jp2=min(jp2,N2(4,3)+1-1)
!jp2=min(jp2,N2(4,nband)+1-1-1)
      do j=jp1,jp2
        nx=nxxp(j)
        nz=nzzp(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,126*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)
    write(*,*) ''

    deallocate(nxxu,nzzu,nxxv,nzzv,nxxp,nzzp)
    deallocate(N2)

  else

    ! The procs receive u1, u2, u3 and p
    ! Those variables are stored in u1PL, u2PL, u3PL and ppPL
    ! The procs only stores the planes they have to compute
    !  (jgal is the local name of planelim for j)

    allocate(N2(4,0:4))
    call MPI_RECV(N2,4*(3+2),MPI_INTEGER,0,121*myid,MPI_COMM_WORLD,status,ierr)
    allocate(nxxu(N2(4,0):N2(4,3)+1),nzzu(N2(4,0):N2(4,3)+1))
    allocate(nxxv(N2(3,0):N2(3,3)+1),nzzv(N2(3,0):N2(3,3)+1))
    nxxu(N2(4,0))=N2(1,3)+2
    nzzu(N2(4,0))=N2(2,3)
    nxxv(N2(3,0))=N2(1,3)+2
    nzzv(N2(3,0))=N2(2,3)

    do j=N2(4,0)+1,N2(4,3)
      nxxu(j)=N2(1,3)+2
      nzzu(j)=N2(2,3)
    end do
    do j=N2(3,0)+1,N2(3,3)
      nxxv(j)=N2(1,3)+2
      nzzv(j)=N2(2,3)
    end do
    ! end do
    nxxu(N2(4,3)+1)=N2(1,3)+2
    nzzu(N2(4,3)+1)=N2(2,3)
    nxxv(N2(3,3)+1)=N2(1,3)+2
    nzzv(N2(3,3)+1)=N2(2,3)
   
    ju1=jgal(2,1)
    ju2=jgal(2,2)
    if (jgal(2,2)==nyu) then
      ju2=jgal(2,2)+1
    end if
    ju1=max(ju1,N2(4,0))
    ju2=min(ju2,N2(4,3)+1)
    jv1=jgal(1,1)
    jv2=jgal(1,2)
    if (jgal(1,2)==nyv.and.myid==np-1) then
      jv2=jgal(1,2)+1
    end if
    jv1=max(jv1,N2(3,0))
    jv2=min(jv2,N2(3,3)+1)
!jv2=min(jv2,N2(3,3)+1-1)
    jp1=jgal(3,1)
    jp2=jgal(3,2)
    if (jgal(3,2)==nyu-1.and.myid==np-1) then
      jp2=jgal(3,2)+1
    end if
    jp1=max(jp1,N2(4,0)+1)
    jp2=min(jp2,N2(4,3)+1-1)
!jp2=min(jp2,N2(4,nband)+1-1-1)
    !!!!!!!!!!!!!!    u1    !!!!!!!!!!!!!!
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,123*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(u1PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    !!!!!!!!!!!!!!    u2    !!!!!!!!!!!!!!
    do j=jv1,jv2
      nx=nxxv(j)
      nz=nzzv(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,124*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(u2PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    !!!!!!!!!!!!!!    u3    !!!!!!!!!!!!!!
    call MPI_BCAST(recv_flg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))

      if (recv_flg == 1) then
        buffSR(:,:) = 0d0
        call MPI_RECV(buffSR(:,1), nx, MPI_REAL8, 0, 125*myid, MPI_COMM_WORLD, status, ierr)
        call buff_to_u(u3PL(1,1,j), buffSR, nx, nz, igal, kgal)
      else 
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,125*myid,MPI_COMM_WORLD,status,ierr)
        call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
      end if 
  
      deallocate(buffSR)
    end do
    !!!!!!!!!!!!!!    p     !!!!!!!!!!!!!!
    do j=jp1,jp2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,126*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(ppPL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    deallocate(nxxu,nzzu,nxxv,nzzv)
    deallocate(N2)
  end if
  !write(*,*) 'finished read in'

end subroutine

subroutine reconstruct_u3(u3,u1,u2,myid)
  use declaration
  implicit none

  integer :: myid
  integer :: column, i, k, j
  complex(8) :: kx, kzF
  complex(8), allocatable :: du2dy(:,:)
  complex(8) :: u1( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )
  complex(8) :: u2( jlim(1,vgrid) : jlim(2,vgrid), columns_num(myid) )
  complex(8) :: u3( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )

  allocate(du2dy(jlim(1,ugrid)+1:jlim(2,ugrid)-1,columns_num(myid)))

  ! reconstruct u3 from continuity, each processor reconstructs their own planes
  ! reconstructed is only calculated from div.u =0 on the internal grid points 
  call der_yu_h(du2dy,u2,myid)  

  do column = 1,columns_num(myid)
    i   = columns_i(column,myid)
    k   = columns_k(column,myid)
    kx  = k1F_x(i)  ! im*kx
    kzF = k1F_z(k)  ! im*kz

        ! reconstructing internal points
    if (abs(kzF) > 1d-12) then
    ! if (kzF /= 0) then
      do j = jlim(1,ugrid)+1,jlim(2,ugrid)-1
        u3(j,column) = (-kx*u1(j,column) - du2dy(j,column)) / kzF
      end do
    end if

    ! reconstructing ghost points
    u3(jlim(1,ugrid),column) = -gridweighting(1) * u3(jlim(1,ugrid)+1,column)
    u3(jlim(2,ugrid),column) = -gridweighting(2) * u3(jlim(2,ugrid)-1,column)
  end do 

end subroutine 



subroutine init_stats(myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    init stats  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer ju1,ju2,jv1,jv2,jp1,jp2,myid, i 

  allocate(Um (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)),U2m (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(Vm (limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)),V2m (limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)))
  allocate(Wm (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)),W2m (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(Pm (limPL_incw(pgrid,1,myid):limPL_incw(pgrid,2,myid)),P2m (limPL_incw(pgrid,1,myid):limPL_incw(pgrid,2,myid)))
  allocate(wxm(limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)),wx2m(limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)))
  allocate(UVm(limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(UWm(limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(VWm(limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)))

  ju1 = -limPL_incw(ugrid,2,myid)+nyu+1
  ju2 = -limPL_incw(ugrid,1,myid)+nyu+1
  jv1 = -limPL_incw(vgrid,2,myid)+nyv+1
  jv2 = -limPL_incw(vgrid,1,myid)+nyv+1
  jp1 = -limPL_incw(pgrid,2,myid)+nyu
  jp2 = -limPL_incw(pgrid,1,myid)+nyu

  !write(6,*) "finished if statement", myid

  Um    = 0d0
  U2m   = 0d0
  Vm    = 0d0
  V2m   = 0d0
  Wm    = 0d0
  W2m   = 0d0
  Pm    = 0d0
  P2m   = 0d0
  UVm   = 0d0
  UWm   = 0d0
  VWm   = 0d0
  wxm   = 0d0
  wx2m  = 0d0
  istat = 0

  !write(6,*) "finished setting to 0", myid
  

 if (myid==0) then
   write(ext4,'(i5.5)') int(t)
   fnameimb = 'stats_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
 end if

end subroutine



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!    nonlin read old    !!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! this is the slightly adapted version of nonlinRead from before the weights
! subroutine nonlinRead(myid)
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !!!!!!!!!!!!    read in important interactions  !!!!!!!!!!!!!!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   use declaration
!     implicit none

!     integer j, jread, x, type, length, jlow, jupp, len, i, myid
!     character*2 extj
!     integer, allocatable :: tmpInt(:)

!     jlow = limPL_excw(ugrid,1,myid)
!     jupp = limPL_excw(ugrid,2,myid)

!     ! write(6,*) "jlow", jlow, "jupp", jupp
!     write(*,*) 'rank', myid, 'limPL_excw ugrid=', limPL_excw(ugrid,1,myid), limPL_excw(ugrid,2,myid)



!     if (jlow > (Ngal(3,nband)+1)/2) then !! bad :/
!       jlow = jlow-1
!     end if
!     jupp = min(jupp,Ngal(3,nband))

!     allocate(nonlin(jlow:jupp,9))

!     do j = jlow,jupp
!       if (j > (Ngal(3,nband)+1)/2) then
!         jread = Ngal(3,nband)+1-j
!       else 
!         jread = j
!         write(6,*) "jread", jread
!       end if
!       write(extj,'(i2.2)') jread

!       fnameima = trim(dirlist)//trim(heading)//extj//'.dat'

!       ! write(*,*) 'DEBUG rank', myid, ' j=', j, ' jread=', jread, ' extj=[', extj, ']'

      
!       open(40, file=fnameima, form='unformatted',access='stream', status='old')

!       allocate(tmpInt(3))
!       read(40) tmpInt
!       if (tmpInt(1)/= N(1,nband)/2-1) then
!         write(*,*) "nx number does not agree between list and simulation"
!         stop
!       elseif (tmpInt(2)/= N(2,nband)/2-1) then
!         write(*,*) "nz number does not agree between list and simulation"
!         stop
!       elseif (tmpInt(3)/= jread) then
!         write(*,*) "something wrong with the j index of the list"
!         stop
!       end if
!       deallocate(tmpInt)
      
!       allocate(tmpInt(2))
      

!       ! weight(j)=0
!       do x = 1,9
!         read(40) tmpInt
!         type = tmpInt(1)
!         length = tmpInt(2)

!         ! write(6,*) "length", length
!         allocate(nonlin(j,type)%list(length,4))
!         read(40) nonlin(j,type)%list

!         ! write(6,*) "read done rank", myid, " j=", j
!         ! do len = 1,length
!         !   read(40) nonlin(j,type)%list(len,:)
!         ! end do 
!         ! write(6,*) 'rank', myid
!       end do

!       ! write(6,*) "j=", j, "jlow", jlow, "jupp", jupp 
!       ! write(6,*) 'rank', myid, 'j=', j, 'jlow=', jlow, 'jupp=', jupp

!       deallocate(tmpInt)
!       close(40)
!     end do


!   end subroutine
! 

! this is joys original version of nonlin read  
subroutine nonlinRead
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!    read in important interactions  !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! currently quite bad, only works for bands of the same size
  ! also file read with access = stream is bad
  ! also some files are read multiple times due to upper and lower channel
  ! also the load in each processor is not perfectly balanced

  use declaration
  implicit none

  integer j, jread, x, type, length, jlow, jupp, len, i
  character*2 extj
  integer, allocatable :: tmpInt(:)

  jlow = min(jgal(vgrid, 1), jgal(ugrid, 1))
  jupp = max(jgal(vgrid, 2), jgal(ugrid, 2))

  !write(6,*) "jgal", size(jgal,1), size(jgal,2)

  ! ! jgal is (3,2)
  ! write(6,*) "jgal ="
  ! do i = 1, 3
  !     write(6,*) jgal(i,1), jgal(i,2)
  ! end do

  ! write(6,*) "ngal", size(ngal,1), size(ngal,2)
  ! write(6,*) "ngal:"
  ! do i = 1,4
  !   write(6,*) ngal(i,0:4)
  ! end do


  if (jlow > (nyv+1)/2) then !! bad :/
    jlow = jlow-1
  end if
  jupp = min(jupp,nyv)

  allocate(nonlin(jlow:jupp,9))

  do j = jlow,jupp
    if (j > (nyv+1)/2) then
      jread = nyv+1-j
    else
      jread = j
    end if
    write(extj,'(i2.2)') jread

    fnameima = trim(dirlist)//trim(heading)//extj//'.dat'

    open(40, file=fnameima, form='unformatted',access='stream', status='old')

    allocate(tmpInt(3))
    read(40) tmpInt
    if (tmpInt(1)/= Nspec_x/2-1) then
      write(*,*) "nx number does not agree between list and simulation"
      stop
    elseif (tmpInt(2)/= Nspec_z/2-1) then
      write(*,*) "nz number does not agree between list and simulation"
      stop
    elseif (tmpInt(3)/= jread) then
      write(*,*) "something wrong with the j index of the list"
      stop
    end if
    deallocate(tmpInt)

    allocate(tmpInt(2))

    do x = 1,9
      read(40) tmpInt
      type = tmpInt(1)
      length = tmpInt(2)
      allocate(nonlin(j,type)%list(length,4))
      read(40) nonlin(j,type)%list
      ! do len = 1,length
      !   read(40) nonlin(j,type)%list(len,:)
      ! end do
    end do

    deallocate(tmpInt)
    close(40)
  end do


end subroutine



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!    nonlin read but w weights    !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_weights
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!    read in important interactions  !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! currently quite bad, only works for bands of the same size
  ! also file read with access = stream is bad
  ! also some files are read multiple times due to upper and lower channel
  ! also the load in each processor is not perfectly balanced

  use declaration
  implicit none
  


  integer j, jread, x, type, length, jlow, jupp, len, i
  character*2 extj
  integer, allocatable :: tmpInt(:)

  integer, parameter :: CHUNK = 1000000          ! 1e6 ints ~ 4 MB buffer
  integer :: nrem, nread
  integer, allocatable:: dump(:)

  if (.not. allocated(dump)) then
  allocate(dump(CHUNK))
  end if

  jlow = nyu_LB +1 
  jupp = nyu-1


  if (allocated(weight)) deallocate(weight)
  allocate(weight(jlow:jupp+1))
  weight(:) = 0   ! initialise


  do j = jlow,jupp
    if (j > (nyv+1)/2) then
      jread = nyv+1-j
    else 
      jread = j
    end if
    write(extj,'(i2.2)') jread

    fnameima = trim(dirlist)//trim(heading)//extj//'.dat'
    
    open(40, file=fnameima, form='unformatted',access='stream', status='old')

    allocate(tmpInt(3))
    read(40) tmpInt
    if (tmpInt(1)/= Nspec_x/2-1) then
      write(*,*) "nx number does not agree between list and simulation"
      stop
    elseif (tmpInt(2)/= Nspec_z/2-1) then
      write(*,*) "nz number does not agree between list and simulation"
      stop
    elseif (tmpInt(3)/= jread) then
      write(*,*) "something wrong with the j index of the list"
      stop
    end if
    deallocate(tmpInt)
    
    allocate(tmpInt(2))

    ! write(6,*) " check 1"

    weight(j)=0
    do x = 1,9
      read(40) tmpInt
      type = tmpInt(1)
      length = tmpInt(2)
      ! write(6,*) "length", length
      weight(j) = weight(j) + length
      ! write(6,*) "weight", weight(j), j 

      ! if (allocated(nonlin(j,type)%list)) deallocate(nonlin(j,type)%list)
      ! allocate(nonlin(j,type)%list(length,4))

      ! read(40) nonlin(j,type)%list
      ! write(6,*) "read done"
      ! do len = 1,length
      !   read(40) nonlin(j,type)%list(len,:)
      ! end do 

      ! Skip the list data: list is (length,4) integers in the file
      nrem = length
      do while (nrem > 0)
        nread = min(nrem, CHUNK)
        read(40) dump(1:nread)
        nrem = nrem - nread
      end do
      
    end do

    write(6,*) "j=", j, "weight", weight(j)
    

    ! write(6,*) "weight(j)", weight(:)

    deallocate(tmpInt)
    close(40)

  end do

  ! do i= jlow, jupp+1
  !   write(6,*) "j=", j, "weight", weight(j)
  ! end do 

  write(6,*) "check 4"
end subroutine
