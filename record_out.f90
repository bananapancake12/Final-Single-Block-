! module rec_out
!   use declaration
!   implicit none

! contains

    
    ! subroutine u_to_buff(buffSR,u,nx,nz,igal,kgal)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!    U to BUFF   !!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ! Rearrange z-modes before writing and after reading

    !     implicit none

    !     integer nx,nz,igal,kgal
    !     real(8) u(igal,kgal)
    !     real(8) buffSR(nx,nz)
    !     integer i,k,dkk

    !     do k = 1,nz/2
    !     do i = 1,nx
    !         buffSR(i,k) = u(i,k    )
    !     end do
    !     end do
    !     do i = 1,nx
    !     !buffSR(i,nz/2+1) = 0d0
    !     end do
    !     do k = nz/2+1,nz
    !     dkk = kgal-nz
    !     do i = 1,nx
    !         buffSR(i,k) = u(i,k+dkk)
    !     end do
    !     end do

    ! end subroutine

  
    ! subroutine flowrateIm(Qu,u)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!         FLOW RATE          !!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ! Flow rate: integrating 1st Fourier mode over the whole channel

    ! use declaration
    ! implicit none
    ! integer j
    ! real(8) Qu
    ! complex(8) u(nyu_LB:nyu+1)

    ! Qu = 0d0
    ! do j = 1,nn+1
    !     !Qu = Qu + .5d0*(yu(j+1) - yu(j-1))*dreal(u(j) !Old collocated version
    !     Qu = Qu + (yv(j) - yv(j-1))*dreal(u(j)) !Mixed ugrid vgrid
    ! end do

    ! end subroutine
    
    ! subroutine maxvel(u)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!          MAX VEL           !!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! use declaration
    ! implicit none

    ! complex(8) u(nyu_LB:nyu+1)
    ! integer j

    ! Umax = 0d0
    ! do j = 0,nn+2
    !     Umax = max(Umax,dreal(u(j)))
    ! end do

    ! ! do j= 1, 152
    ! !   write(6,*) " Umax", Umax, dreal(u(j)), "j", j
    ! ! end do 

    ! end subroutine

  
    ! subroutine record_out(u1,myid)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!   RECORD OUT   !!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! use declaration
    ! ! use littleharsh_mod
    ! implicit none

    ! include 'mpif.h'             ! MPI variables
    ! integer status(MPI_STATUS_SIZE),ierr,myid

    ! complex(8)   u1( jlim(1,ugrid) : jlim(2,ugrid), columns_num(myid) )
    ! integer nx,nz, i 
    ! integer j,jmax,iproc
    ! real(8), allocatable:: buffSR(:,:)
    ! integer, allocatable:: dummint(:)
    ! real(8) Uslip  

    ! ! Rebuilding N 
    ! if (.not. allocated(N)) allocate(N(4,0:4))


    ! N= 0 
    ! N(1,1:3) = Nspec_x
    ! N(1,4) = -2
    ! N(2,1:3) = Nspec_z
    ! N(3,3) = nyv
    ! N(4,3) = nyu

    ! if (myid == 0) then
    !     write(6,*) "N:"
    !     do i = 1,4
    !     write(6,*) N(i,0:4)
    !     end do
    ! end if 

    ! if (myid/=0) then
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
    !     call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
    !     call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,123*myid,MPI_COMM_WORLD,ierr)
    !     end do
    !     do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
    !     call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
    !     call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,124*myid,MPI_COMM_WORLD,ierr)
    !     end do
    !     do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
    !     call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
    !     call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,125*myid,MPI_COMM_WORLD,ierr)
    !     end do
    !     do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
    !     call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
    !     call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,126*myid,MPI_COMM_WORLD,ierr)
    !     end do

    !     deallocate(buffSR)
    ! else
    !     write(ext4,'(i5.5)') int(10d0*(t))!int(t)!
    !     allocate(dummint(88))
    !     dummint = 0
    !     !!!!!!!!!!!!!    u1    !!!!!!!!!!!!!
    !     fnameima = 'output/u1_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    !     open(10,file=fnameima,form='unformatted')
    !     write(10) t,Re,alp,bet,mpgx,nband,iter,dummint 
    !     write(10) N
    !     write(10) yu,dthetavi,dthdyu
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
    !     call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
    !     write(10) j,1,nx,nz,yu(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     do iproc = 1,np-1
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
    !         call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,123*iproc,MPI_COMM_WORLD,status,ierr)
    !         write(10) j,1,nx,nz,yu(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     end do
    !     close(10)
    !     !!!!!!!!!!!!!    u2    !!!!!!!!!!!!!
    !     fnameima='output/u2_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    !     open(10,file=fnameima,form='unformatted')
    !     write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    !     write(10) N
    !     write(10) yv,dthetavi,dthdyv
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
    !     call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
    !     write(10) j,2,nx,nz,yv(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     do iproc = 1,np-1
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
    !         call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,124*iproc,MPI_COMM_WORLD,status,ierr)
    !         write(10) j,2,nx,nz,yv(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     end do
    !     close(10)
    !     !!!!!!!!!!!!!    u3    !!!!!!!!!!!!!
    !     fnameima = 'output/u3_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    !     open(10,file=fnameima,form='unformatted')
    !     write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    !     write(10) N
    !     write(10) yu,dthetavi,dthdyu
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
    !     call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
    !     write(10) j,3,nx,nz,yu(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     do iproc = 1,np-1
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
    !         call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,status,ierr)
    !         write(10) j,3,nx,nz,yu(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     end do
    !     close(10)
    !     !!!!!!!!!!!!!    p     !!!!!!!!!!!!!
    !     fnameima = 'output/p_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    !     open(10,file=fnameima,form='unformatted')
    !     write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    !     write(10) N
    !     write(10) yu,dthetavi,dthdyu
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
    !     call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
    !     write(10) j,4,nx,nz,yu(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     do iproc = 1,np-1
    !     nx = Nspec_x+2
    !     nz = Nspec_z
    !     allocate(buffSR(nx,nz))
    !     do j = limPL_incw(pgrid,1,iproc),limPL_incw(pgrid,2,iproc)
    !         call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,126*iproc,MPI_COMM_WORLD,status,ierr)
    !         write(10) j,4,nx,nz,yu(j),buffSR
    !     end do
    !     deallocate(buffSR)
    !     end do
    !     close(10)
    !     deallocate(dummint)
    ! end if

    ! if (myid==0) then

    !     Uslip = ((u1(1,1))*(-1d0-yu(0))+(u1(0,1))*(yu(1)+1d0))/(yu(1)-yu(0))

    !     ! call flowrateIm(Qx,u1(nyu_LB,1))
    !     call flowrateIm(Qx,u1(:,1))
    !     ! call maxvel(u1(nyu_LB,1))
    !     call maxvel(u1(:,1))
    !     write(*,*) ''
    !     write(*,*) 'iter',iter
    !     write(*,*) 't   ',t
    !     write(*,*) 'dtv ',dtv
    !     write(*,*) 'dtc ',dtc
    !     write(*,*) 'dt  ',dt
    !     write(*,*) 'err ',err
    !     write(*,*) 'Qx  ',Qx
    !     if (flag_ctpress==0) then
    !     write(*,*) 'QxT ',QxT
    !     write(*,*) 'mpgx',mpgx
    !     write(*,*) 'dpgx',dgx
    !     else
    !     write(*,*) 'mpgx',mpgx
    !     end if
    !     write(*,*) 'Umax',Umax
    !     write(*,*) 'Uslp',Uslip
    ! !    write(*,*) 'utau',utau
    ! ! Save to history file
    !     if (flag_ctpress==0) then
    !     write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,QxT,mpgx,dgx,Umax,Uslip
    !     else
    !     write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,    mpgx,    Umax,Uslip
    !     end if
    !     flush(30)
    ! end if

    ! end subroutine

! end module rec_out