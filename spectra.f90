
! module spectra_mod
!   use declaration
!   implicit none
! contains  
  
  
  subroutine spectra(u1,u2,u2_itp2,u3,p,myid)

    use declaration
    implicit none

    integer myid
    ! complex(8), intent(in) :: u1(jlim(1,ugrid):,:), u2(jlim(1,vgrid):,:), u3(jlim(1,ugrid):,:)
    complex(8) :: u1( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) )
    complex(8) :: u2( jlim(1,vgrid)      : jlim(2,vgrid),      columns_num(myid) ) 
    complex(8) :: u3( jlim(1,ugrid)      : jlim(2,ugrid),      columns_num(myid) ) 
    complex(8) :: p ( jlim(1,pgrid)      : jlim(2,pgrid),      columns_num(myid) )
    complex(8) :: u2_itp2( jlim(1,ugrid):jlim(2,ugrid), columns_num(myid) )

    call buildsp  (spU ,u1,                   ugrid,myid)
    call buildsp  (spV ,u2,                   vgrid,myid)
    call buildsp  (spW ,u3,                   ugrid,myid)
    call buildspUV(spUV,u1,u2_itp2,                 myid)
    call buildsp  (spP ,p ,                   pgrid,myid)

  end subroutine

  subroutine buildsp(sp,u,grid,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!    buildD   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none
    integer j,column,grid,myid
    real(8)    sp(jlim(1,grid):jlim(2,grid),columns_num(myid))
    complex(8) :: u (jlim(1,grid):jlim(2,grid),columns_num(myid))

    ! write(6,*) "meow check"

    do column = 1,columns_num(myid)
      do j = jlim(1,grid),jlim(2,grid)
        sp(j,column) = sp(j,column)+dreal(u(j,column))**2+dimag(u(j,column))**2
      end do
    end do

  end subroutine

  subroutine buildspUV(sp,u,v,myid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!    buildD   !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none
    integer j,column,grid,myid
    real(8)    sp(jlim(1,ugrid):jlim(2,ugrid),columns_num(myid))
    complex(8) u (jlim(1,ugrid):jlim(2,ugrid),columns_num(myid))
    complex(8) v (jlim(1,ugrid):jlim(2,ugrid),columns_num(myid))

    do column = 1,columns_num(myid)
      do j = jlim(1,ugrid),jlim(2,ugrid)
          sp(j,column) = sp(j,column)+dreal(u(j,column))*dreal(v(j,column)) &
  &                                  +dimag(u(j,column))*dimag(v(j,column))
      end do
    end do

  end subroutine

  subroutine write_spect(myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!! WRITE SPECT !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    include 'mpif.h'             ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr,myid

    real(8), allocatable:: buffSPu(:,:,:),buffSPv(:,:,:),buffSPp(:,:,:)
    integer msizeu,msizev,msizep

    if (myid/=0) then
      msizeu = (jlim(2,ugrid)-jlim(1,ugrid)+1)*columns_num(myid)
      msizev = (jlim(2,vgrid)-jlim(1,vgrid)+1)*columns_num(myid)
      msizep = (jlim(2,pgrid)-jlim(1,pgrid)+1)*columns_num(myid)
      call MPI_SEND(spU ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(spV ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(spW ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(spUV,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(spP ,msizep,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
    else
      allocate(buffSPu(0:Nspec_x/2,1:Nspec_z/2+1,jlim(1,ugrid):jlim(2,ugrid)))
      allocate(buffSPv(0:Nspec_x/2,1:Nspec_z/2+1,jlim(1,vgrid):jlim(2,vgrid)))
      allocate(buffSPp(0:Nspec_x/2,1:Nspec_z/2+1,jlim(1,pgrid):jlim(2,pgrid)))

      fnameimb='output/spU_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
      call recvwrspec(spU ,buffSPu,ugrid,myid,status,ierr)
      fnameimb='output/spV_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
      call recvwrspec(spV ,buffSPv,vgrid,myid,status,ierr)
      fnameimb='output/spW_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
      call recvwrspec(spW ,buffSPu,ugrid,myid,status,ierr)
      fnameimb='output/spUV_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
      call recvwrspec(spUV,buffSPu,ugrid,myid,status,ierr)
      fnameimb='output/spP_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
      call recvwrspec(spP ,buffSPp,pgrid,myid,status,ierr)
      deallocate(buffSPu)
      deallocate(buffSPv)
      deallocate(buffSPp)
    end if

    spU  = 0d0
    spV  = 0d0
    spW  = 0d0
    spUV = 0d0
    spP  = 0d0

  end subroutine

  subroutine recvwrspec(spX,buffSP,grid,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!    recv+write spec    !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    include 'mpif.h'             ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr

    integer i,k,kk,j,iproc,column,grid,myid
    integer msize
    real(8)    spX(jlim(1,grid):jlim(2,grid),columns_num(myid))
    !The following +1 in 1:N(2,midband)/2+1 shouldn't be there but is consistent with the postprocessing....
    real(8) buffSP(0:Nspec_x/2,1:Nspec_z/2+1,jlim(1,grid):jlim(2,grid))
    real(8), allocatable:: buffrecv(:,:)
    integer, allocatable:: dummint(:)

    buffSP = 0d0

    do j = jlim(1,grid),jlim(2,grid)
      do column = 1,columns_num(myid)
        i = columns_i(column,myid)
        k = columns_k(column,myid)
        if (k > Nspec_z/2) then
          kk = Nspec_z+2-k
        else
          kk = k
        end if
        ! write(6,*) "kk", kk, j 
        buffSP(i,kk,j) = spX(j,column)
      end do
    end do

    do iproc = 1,np-1
      msize = (jlim(2,grid)-jlim(1,grid)+1)*columns_num(iproc)
      allocate(buffrecv(jlim(1,grid):jlim(2,grid),columns_num(iproc)))
      call MPI_RECV(buffrecv,msize,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,status,ierr)
      do j = jlim(1,grid),jlim(2,grid)
        do column = 1,columns_num(iproc)
          i = columns_i(column,iproc)
          k = columns_k(column,iproc)
          if (k > Nspec_z/2) then
            kk = Nspec_z+2-k
          else
            kk = k
          end if
          buffSP(i,kk,j) = buffSP(i,kk,j)+buffrecv(j,column)
        end do
      end do
      deallocate(buffrecv)
    end do

    do j = jlim(1,grid),jlim(2,grid)
      do k = 1,Nspec_z/2
        do i = 1,Nspec_x/2
          buffSP(i,k,j) = 2d0*buffSP(i,k,j)
        end do
      end do
    end do

    open(10,file=fnameimb,form='unformatted')
    allocate(dummint(88))
    dummint = 0
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    deallocate(dummint)
    write(10) N
    write(10) yu,dthetai,dthdyu
    write(10) yv,dthetai,dthdyv
    write(10) istat
    write(10) buffSP
    close(10)

  end subroutine


! end module spectra_mod