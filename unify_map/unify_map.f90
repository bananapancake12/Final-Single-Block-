program unify_map

  use declaration
  implicit none

  integer :: jpl
  integer :: isamp

  ! Read setup + allocate arrays
  call start

  ! ------------------------------------------------------------
  ! Loop over all heights
  ! ------------------------------------------------------------
  do jpl = 1, PLoINum

    write(*,*) 'Processing height ', jpl, ' of ', PLoINum

    ! Reset accumulators for this height
    convs_uu = 0d0
    convs_uv = 0d0
    convs_uw = 0d0

    convs_vu = 0d0
    convs_vv = 0d0
    convs_vw = 0d0

    convs_wu = 0d0
    convs_wv = 0d0
    convs_ww = 0d0

    nsamp_total = 0

    ! ----------------------------------------------------------
    ! Loop over all samples
    ! ----------------------------------------------------------
    do isamp = 1, nsamp
      call read_map_file(jpl,isamp)

    end do

    ! ----------------------------------------------------------
    ! Write merged file for this height
    ! ----------------------------------------------------------
    call write_map_file(jpl)

  end do

  write(*,*) ''
  write(*,*) 'DONE'
  write(*,*) ''

end program



subroutine start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    START    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer :: i
  real(8) :: dum1, dum2
  character(len=200) :: fn

  pi = 4d0*datan(1d0)

  ! ------------------------------------------------------------
  ! Read unify_map.input
  ! ------------------------------------------------------------
  open(40,file='unify_map.input',form='formatted')
  read(40,31) filelist
  read(40,20) firstgood
  read(40,31) fildir
  read(40,20) nx
  read(40,20) nz
  read(40,20) ny
  close(40)

  ! ------------------------------------------------------------
  ! Read list of sample suffixes
  ! ------------------------------------------------------------
  open(40,file=filelist,form='formatted')
  read(40,37) nsamp
  allocate(extL(nsamp))
  do i = 1, nsamp
    read(40,34) extL(i)
  end do
  close(40)

  ! ------------------------------------------------------------
  ! Read PlaneOfInterest.txt
  ! ------------------------------------------------------------
  fn = trim(fildir)//'postPars/PlaneOfInterest.txt'
  open(40,file=fn,form='formatted')

  read(40,*) PLoINum
  allocate(NYoI(PLoINum))
  allocate(yPLoI(PLoINum))

  do i = 1, PLoINum
    read(40,*) NYoI(i), yPLoI(i)
  end do

  close(40)

  ! ------------------------------------------------------------
  ! Read ModesOfInterest.txt
  ! ------------------------------------------------------------
  fn = trim(fildir)//'postPars/ModesOfInterest.txt'
  open(40,file=fn,form='formatted')

  read(40,*) MoINumX
  allocate(NXfib(MoINumX))
  allocate(NXoI(MoINumX))

  do i = 1, MoINumX
      read(40,*) dum1, dum2, dum1, NXoI(i)
  end do

  read(40,*)
  read(40,*) MoINumZ
  allocate(NZfib(MoINumZ))
  allocate(NZoI(MoINumZ))

  do i = 1, MoINumZ
      read(40,*) dum1, dum2, dum1, NZoI(i)
  end do

  close(40)

  MoINumt = MoINumX * MoINumZ

  ! ------------------------------------------------------------
  ! Allocate one-height-at-a-time accumulators
  ! ------------------------------------------------------------
  allocate(convs_uu(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_uv(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_uw(4,4,MoINumX,MoINumZ,MoINumt))

  allocate(convs_vu(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_vv(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_vw(4,4,MoINumX,MoINumZ,MoINumt))

  allocate(convs_wu(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_wv(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_ww(4,4,MoINumX,MoINumZ,MoINumt))

  allocate(convs_uu_tmp(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_uv_tmp(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_uw_tmp(4,4,MoINumX,MoINumZ,MoINumt))

  allocate(convs_vu_tmp(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_vv_tmp(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_vw_tmp(4,4,MoINumX,MoINumZ,MoINumt))

  allocate(convs_wu_tmp(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_wv_tmp(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(convs_ww_tmp(4,4,MoINumX,MoINumZ,MoINumt))

  convs_uu = 0d0
  convs_uv = 0d0
  convs_uw = 0d0

  convs_vu = 0d0
  convs_vv = 0d0
  convs_vw = 0d0

  convs_wu = 0d0
  convs_wv = 0d0
  convs_ww = 0d0

  nsamp_total = 0

  write(*,*) 'PLoINum ', PLoINum
  write(*,*) 'MoINumX ', MoINumX
  write(*,*) 'MoINumZ ', MoINumZ
  write(*,*) 'MoINumt ', MoINumt
  write(*,*) 'nsamp   ', nsamp
  write(*,*) ''

20 FORMAT(10X,I10)
31 FORMAT(10X,A60)
34 FORMAT(A5)
37 FORMAT(I5)

end subroutine

subroutine read_map_file(jpl,isamp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  READ MAP  !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer :: jpl, isamp
  integer :: i, k, m, s
  integer :: nsamp_file
  integer :: nxt_file, nzt_file
  integer :: NYoI_file
  real(8) :: yPLoI_file

  convs_uu_tmp = 0d0
  convs_uv_tmp = 0d0
  convs_uw_tmp = 0d0

  convs_vu_tmp = 0d0
  convs_vv_tmp = 0d0
  convs_vw_tmp = 0d0

  convs_wu_tmp = 0d0
  convs_wv_tmp = 0d0
  convs_ww_tmp = 0d0

  ! ------------------------------------------------------------
  ! Build filename
  ! ------------------------------------------------------------
  write(extnkx,'(i3.3)') MoINumX
  write(extnkz,'(i3.3)') MoINumZ
  write(extny ,'(i3.3)') int(yPLoI(jpl))

  fnameima = trim(fildir)//'TRD_'// &
             extnkx//'_'//extnkz//'_'//extny//'_t'// &
             trim(extL(isamp))//'.dat'

  open(20,file=fnameima,form='unformatted')
  read(20) Re, alp, bet, mpgx, Ntx, Ntz
  read(20) NYoI_file, yPLoI_file
  read(20) NXfib, NZfib
  read(20) NXoI, NZoI
  read(20) nsamp_file

  nsamp_total = nsamp_total + nsamp_file

  ! ------------------------------------------------------------

  do i = 1, MoINumX
    do k = 1, MoINumZ

    read(20)
    read(20) nxt_file, nzt_file

    if (nxt_file /= NXoI(i) .or. nzt_file /= NZoI(k)) then
      write(*,*) 'Target mode mismatch in file ', trim(fnameima)
      write(*,*) 'At i,k = ', i, k
      write(*,*) 'Expected NXoI,NZoI = ', NXoI(i), NZoI(k)
      write(*,*) 'Read     nxt,nzt   = ', nxt_file, nzt_file
      stop
    end if

      do m = 1, 4
        do s = 1, 4

          ! -------------------------
          ! uu uv uw
          ! -------------------------
          read(20) &
            convs_uu_tmp(m,s,i,k,:), &
            convs_uv_tmp(m,s,i,k,:), &
            convs_uw_tmp(m,s,i,k,:)

          ! -------------------------
          ! vu vv vw
          ! -------------------------
          read(20) &
            convs_vu_tmp(m,s,i,k,:), &
            convs_vv_tmp(m,s,i,k,:), &
            convs_vw_tmp(m,s,i,k,:)

          ! -------------------------
          ! wu wv ww
          ! -------------------------
          read(20) &
            convs_wu_tmp(m,s,i,k,:), &
            convs_wv_tmp(m,s,i,k,:), &
            convs_ww_tmp(m,s,i,k,:)

        end do
      end do

    end do
  end do

  close(20)

  ! --------------------   Accumulate    ---------------------


  convs_uu = convs_uu + convs_uu_tmp
  convs_uv = convs_uv + convs_uv_tmp
  convs_uw = convs_uw + convs_uw_tmp

  convs_vu = convs_vu + convs_vu_tmp
  convs_vv = convs_vv + convs_vv_tmp
  convs_vw = convs_vw + convs_vw_tmp

  convs_wu = convs_wu + convs_wu_tmp
  convs_wv = convs_wv + convs_wv_tmp
  convs_ww = convs_ww + convs_ww_tmp

end subroutine


subroutine write_map_file(jpl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  WRITE MAP  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer :: jpl
  integer :: i, k, f, s
  character(len=256) :: fnameList

  write(extnkx,'(i3.3)') MoINumX
  write(extnkz,'(i3.3)') MoINumZ
  write(extny ,'(i3.3)') int(yPLoI(jpl))

  fnameList = 'TRD_'// &
              extnkx//'_'//extnkz//'_'//extny//'.dat'

  open(unit=50, file=fnameList, form='unformatted')

  write(50) Re, alp, bet, mpgx, Ntx, Ntz
  write(50) NYoI(jpl), yPLoI(jpl)
  write(50) NXfib(:), NZfib(:)
  write(50) NXoI(:), NZoI(:)
  write(50) nsamp_total

  do i = 1, MoINumX
    do k = 1, MoINumZ

      write(50)
      write(50) NXoI(i), NZoI(k)

      do f = 1, 4
        do s = 1, 4

          write(50) &
            convs_uu(f,s,i,k,:), &
            convs_uv(f,s,i,k,:), &
            convs_uw(f,s,i,k,:)

          write(50) &
            convs_vu(f,s,i,k,:), &
            convs_vv(f,s,i,k,:), &
            convs_vw(f,s,i,k,:)

          write(50) &
            convs_wu(f,s,i,k,:), &
            convs_wv(f,s,i,k,:), &
            convs_ww(f,s,i,k,:)

        end do
      end do

    end do
  end do

  close(50)

end subroutine