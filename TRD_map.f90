
subroutine map(myid)

  use declaration
  implicit none 
  include 'mpif.h'
  integer status(MPI_STATUS_SIZE),ierr,myid
  
  integer :: j, i, k, NRplxz


  real(8), allocatable :: u1_ind(:,:), u2_ind(:,:), u3_ind(:,:)
  real(8), allocatable :: u1PL_map(:,:,:),u2PL_map(:,:,:), u3PL_map(:,:,:)
  integer, parameter :: LowChn = 1, UppChn = 2
  integer :: whiChn

  real(8), allocatable :: convs_uu_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_uv_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_uw_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_vu_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_vv_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_vw_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_wu_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_wv_tmp(:,:,:,:,:,:)
  real(8), allocatable :: convs_ww_tmp(:,:,:,:,:,:)

  integer :: it_moi, jpl, jbf, idx
  integer :: k_ind, IkkcNeg
  integer :: ric, ncs
  integer :: Ntx, Ntz, f, s, jL, jU

  integer ::  n_here


  integer :: ii, kk, ind, sigCase
  integer :: NxUpp, NzUpp, jplex
  integer :: buffInt   ! if you use merge() with 1d0 below
  

  integer :: xlim(4), zlim(4)

  real(8) :: kc_x, kc_z
  real(8) :: u1C_Re, u1C_Im
  real(8) :: u2C_Re, u2C_Im
  real(8) :: u3C_Re, u3C_Im


  integer :: icsub, iia, kka, iib, kkb, dk2, di2, nx, nz
  integer :: IkkNeg, IiiNeg
  character(len=3)  :: extnkx, extnkz, extny, extmod
  character(len=256):: fnameList
  character(len=256):: output, map_output

  ! integer :: n_planes, idx
  ! integer, allocatable :: jpl_listL(:)

  integer :: max_iia, max_kka

  real(8), allocatable:: buffSR(:,:)
  integer :: planeSize, src, dst

  real(8), allocatable :: planeBuf(:,:,:,:,:), recvPlane(:,:,:,:,:)

  allocate( convs_uu_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_uv_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_uw_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_vu_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_vv_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_vw_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_wu_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_wv_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )
  allocate( convs_ww_tmp(4,4,PLoINum, MoINumX, MoINumZ, MoINumt) )

  convs_uu_tmp = 0.0d0
  convs_uv_tmp = 0.0d0
  convs_uw_tmp = 0.0d0
  convs_vu_tmp = 0.0d0
  convs_vv_tmp = 0.0d0
  convs_vw_tmp = 0.0d0
  convs_wu_tmp = 0.0d0
  convs_wv_tmp = 0.0d0
  convs_ww_tmp = 0.0d0


  max_iia = 0
  max_kka = 0



  NRplxz = (Nspec_x+2) * Nspec_z
  Ntx = Nspec_x
  Ntz = Nspec_z


  allocate(u1_ind(NRplxz, jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(u2_ind(NRplxz, jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(u3_ind(NRplxz, jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(u1PL_map(Nspec_x+2,Nspec_z,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(u2PL_map(Nspec_x+2,Nspec_z,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(u3PL_map(Nspec_x+2,Nspec_z,jgal(ugrid,1)-1:jgal(ugrid,2)+1))

  nsamp = nsamp + 1


  ! removing the band of zeroes in the antialiasing region
  nx = Nspec_x+2
  nz = Nspec_z

  allocate(buffSR(nx,nz))

  do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
    call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
    u1PL_map(:,:,j) = buffSR
  end do

  do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
    call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
    u2PL_map(:,:,j) = buffSR
  end do

  do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
    call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
    u3PL_map(:,:,j) = buffSR
  end do

  deallocate(buffSR)

  ! flattening the array to match the indexing of post processing script 
  
  do j = jgal(ugrid,1)-1, jgal(ugrid,2)+1
    do k = 1, Nspec_z
      do i = 1, Nspec_x+2
        u1_ind( (Nspec_x+2)*(k-1) + i , j ) = u1PL_map(i,k,j)
        u3_ind( (Nspec_x+2)*(k-1) + i , j ) = u3PL_map(i,k,j)
      end do
    end do
  end do

  do j = jgal(vgrid,1)-1, jgal(vgrid,2)+1
    do k = 1, Nspec_z
      do i = 1, Nspec_x+2
        u2_ind( (Nspec_x+2)*(k-1) + i , j ) = u2PL_map(i,k,j)
      end do
    end do
  end do



  ! ! ----- Building lists of what PLoI each rank owns ----- !
  ! n_planes = 0 
  ! do jpl = 1, PLoINum
  !   j   = NYoI(jpl)
  !   if (j >= jgal(ugrid,1) .and. j <= jgal(ugrid,2)) then
  !     n_planes = n_planes + 1
  !   end if 
  ! end do 

  ! allocate(jpl_listL(n_planes))

  ! idx = 0
  ! do jpl = 1, PLoINum
  !   j   = NYoI(jpl)
  !   if (j >= jgal(ugrid,1) .and. j <= jgal(ugrid,2)) then 
  !     idx = idx + 1
  !     jpl_listL(idx) = jpl
  !   end if 
  ! end do 

  ! ! --- Upper list: decide ownership using mirrored physical plane jU
  ! n_planesU = 0
  ! do jpl = 1, PLoINum
  !   jU = nyf - NYoI(jpl) -1
  !   ! write(6,*) "JU", jU, "nyf", nyf
  !   if (jU >= jgal(ugrid,1) .and. jU <= jgal(ugrid,2)) n_planesU = n_planesU + 1
  ! end do
  
  ! if (allocated(jpl_listU)) deallocate(jpl_listU)
  ! allocate(jpl_listU(n_planesU))

  ! idx = 0
  ! do jpl = 1, PLoINum
  !   jU = nyf - NYoI(jpl) -1
  !   if (jU >= jgal(ugrid,1) .and. jU <= jgal(ugrid,2)) then
  !     idx = idx + 1
  !     jpl_listU(idx) = jpl
  !     ! write(6,*) "jpl_listU(idx)", jpl_listU(idx)
  !   end if
  ! end do


  ! do i = 0, np-1
  !   if (myid == i) then
  !     write(6,*) "Rank", myid, "owns LOW planes:", jpl_listL !, jgal(ugrid,1), jgal(ugrid,2)
  !     write(6,*) "Rank", myid, "owns UPP planes:", jpl_listU !, jgal(ugrid,1), jgal(ugrid,2)
  !   end if
  ! end do

  ! After jpl_list is built, may need to send and recive any planes we dont already own...
  ! ---- ugrid halo exchange for u1,u3 (always, if neighbour exists) ----
  if (myid /= np-1) then
    call MPI_SENDRECV( u1_ind(:, jgal(ugrid,2)),   NRplxz, MPI_REAL8, myid+1, 101, &
                      u1_ind(:, jgal(ugrid,2)+1), NRplxz, MPI_REAL8, myid+1, 102, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
    call MPI_SENDRECV( u3_ind(:, jgal(ugrid,2)),   NRplxz, MPI_REAL8, myid+1, 201, &
                      u3_ind(:, jgal(ugrid,2)+1), NRplxz, MPI_REAL8, myid+1, 202, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
  end if

  if (myid /= 0) then
    call MPI_SENDRECV( u1_ind(:, jgal(ugrid,1)),   NRplxz, MPI_REAL8, myid-1, 102, &
                      u1_ind(:, jgal(ugrid,1)-1), NRplxz, MPI_REAL8, myid-1, 101, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
    call MPI_SENDRECV( u3_ind(:, jgal(ugrid,1)),   NRplxz, MPI_REAL8, myid-1, 202, &
                      u3_ind(:, jgal(ugrid,1)-1), NRplxz, MPI_REAL8, myid-1, 201, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
  end if

  if (myid /= np-1) then
    call MPI_SENDRECV( u2_ind(:, jgal(vgrid,2)),   NRplxz, MPI_REAL8, myid+1, 301, &
                      u2_ind(:, jgal(vgrid,2)+1), NRplxz, MPI_REAL8, myid+1, 302, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
  end if

  if (myid /= 0) then
    call MPI_SENDRECV( u2_ind(:, jgal(vgrid,1)),   NRplxz, MPI_REAL8, myid-1, 302, &
                      u2_ind(:, jgal(vgrid,1)-1), NRplxz, MPI_REAL8, myid-1, 301, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
  end if

  RBf_u1 = 0.0d0
  RBf_u2 = 0.0d0
  RBf_u3 = 0.0d0

  ! if(myid == 0) then
  !   do i = 1, PLoINumEx
  !     write(6,*) "jList_Buff", jList_Buff(i)
  !   end do
  ! end if

  do whiChn = LowChn, UppChn
    do jplex = 1, PLoINumEx
      jL = jList_Buff(jplex)

      ! --- mirror indices that match TRDv9 ---
      if (whiChn == LowChn) then
        jU = jL                 ! u1/u3 index
        j  = jL                 ! u2 index (same on low side)
      else
        jU = nyf - jL           ! u1/u3 mirror (TRDv9)
        j  = (nyf - 1) - jL      ! u2 mirror uses nnc=nyf-1 (TRDv9)
        ! write(6,*) "j", j, "jL", jL-1
      end if

      ! u1/u3 live on ugrid -> use jU
      if (jU >= jgal(ugrid,1)-1 .and. jU <= jgal(ugrid,2)+1) then
        RBf_u1(:, jplex, whiChn) = u1_ind(:, jU)
        RBf_u3(:, jplex, whiChn) = u3_ind(:, jU)
      end if

      ! u2 lives on vgrid -> use j (nnc-mirror in upper)
      if (jplex <= PLoINumEx-1) then
        if (j >= jgal(vgrid,1)-1 .and. j <= jgal(vgrid,2)+1) then
          RBf_u2(:, jplex, whiChn) = u2_ind(:, j)
        end if
      end if
    end do
  end do

  RBf_u2(:,:,UppChn) = -RBf_u2(:,:,UppChn)

  ! if( myid == 0 ) then
  !   ! write(*,*) "RBf_u1", RBf_u1( 1:100, 1, LowChn)
  !   ! write(6,*) "RBf_u2", RBf_u2( 1:100, 1, LowChn)
  !   write(6,*) "RBf_u3", RBf_u3( 1:100, 1, LowChn)
  ! end if 


  ! if( myid == 0 ) then

  !   write(6,*) "RBf_u2(:,jplex,LowChn)", RBf_u2(1:100,1,LowChn)
  !   write(6,*) "RBf_u3(:,jplex,LowChn)", RBf_u3(1:100,1,LowChn)
  !   write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  !   write(6,*) "RBf_u3(:,jplex,UppChn)", RBf_u3(1:100,1,UppChn)
  ! end if 

  ! write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  ! write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  ! write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  ! write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  ! if( myid == 7 ) then
  ! !   write(6,*) "RBf_u2(:,jplex,LowChn)", RBf_u2(1:100,1,LowChn)
  !   write(6,*) "RBf_u3(:,jplex,LowChn)", RBf_u3(1:100,1,LowChn)
  !   write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  !   write(6,*) "RBf_u3(:,jplex,UppChn)", RBf_u3(1:100,1,UppChn)
  ! end if 
  ! ---- Calculating and writing ----
  if (myid ==0 ) then
    write(*,*) 'Calculating'
  end if 

  Do whiChn = LowChn, UppChn
    it_moi = 0
    
    if (whiChn == LowChn) then
      n_here = n_planesL
    else
      n_here = n_planesU
    end if

    do idx = 1, n_here

      if (whiChn == LowChn) then
        jpl = jpl_listL(idx)
      else
        jpl = jpl_listU(idx)
      end if

      j   = NYoI(jpl)
      jbf = buffIndj(jpl)

      ! if (whiChn == LowChn) then
      !     write(*,"(4X,A22,I3,1X,A2,1X,I3)") 'Lower channel: Plane #', jpl, 'of', PLoINum
      ! else
      !     write(*,"(4X,A22,I3,1X,A2,1X,I3)") 'Upper channel: Plane #', jpl, 'of', PLoINum
      ! end if


      !Storing u, v, w and the d/dy for that plane
      u1pl_tmp =   RBf_u1(:,jbf  ,whiChn)
      s1pl_tmp = ( RBf_u1(:,jbf+1,whiChn) - RBf_u1(:,jbf-1,whiChn) ) / 2.d0 * (dthdyu(j)*ddthetavi)
      u2pl_tmp = ( RBf_u2(:,jbf  ,whiChn) * (yu(j)-yv(j-1)) &
                & + RBf_u2(:,jbf-1,whiChn) * (yv(j  )-yu(j)) ) / (yv(j)-yv(j-1))
      s2pl_tmp = ( RBf_u2(:,jbf  ,whiChn) - RBf_u2(:,jbf-1,whiChn) )        * (dthdyu(j)*ddthetavi)
      u3pl_tmp =   RBf_u3(:,jbf  ,whiChn)
      s3pl_tmp = ( RBf_u3(:,jbf+1,whiChn) - RBf_u3(:,jbf-1,whiChn) ) / 2.d0 * (dthdyu(j)*ddthetavi)


      do i = 1, MoINumX
        do k = 1, 2*MoINumZ
          ! if (mod(k, 2) == 1) cycle !!!!!!!!!!!!!!!!!!!!!!!!!
          
          it_moi = it_moi + 1
          k_ind = (k + 1) / 2
          IkkcNeg = 1-mod(k,2)
          ! write(*,"(8X,A6,I3,1X,A2,1X,I3)") 'mode #', (i-1)*MoINumZ*2+k, 'of', MoINumX*2*MoINumZ
          ric    = indMoI(i,k)%SubMatA%RIndC
          kc_x   = indMoI(i,k)%NXoI * alp
          kc_z   = indMoI(i,k)%NZoI * bet
          
          u1C_Re = u1pl_tmp(ric  )
          u1C_Im = u1pl_tmp(ric+1)
          u2C_Re = u2pl_tmp(ric  )
          u2C_Im = u2pl_tmp(ric+1)
          u3C_Re = u3pl_tmp(ric  )
          u3C_Im = u3pl_tmp(ric+1)

          xlim(1:4) = [Ntx/2-1, Ntx/2-1-i, Ntx/2-1-i, Ntx/2-1]
          zlim(1:4) = [Ntz/2-1, Ntz/2-1, Ntz/2-1 - k_ind,Ntz/2-1 - k_ind]

          ncs           =           indMoI(i,k)%SubMatA%nsub
          ka_x(  1:ncs) = alp *     indMoI(i,k)%SubMatA%ni(  :)
          ka_z(  1:ncs) = bet *     indMoI(i,k)%SubMatA%nk(  :)
          u1A_Re(1:ncs) = u1pl_tmp( indMoI(i,k)%SubMatA%RInd(:)   )
          u1A_Im(1:ncs) = u1pl_tmp( indMoI(i,k)%SubMatA%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatA%CC(  :) ! -kx terms: conjugate of +kx
          s1A_Re(1:ncs) = s1pl_tmp( indMoI(i,k)%SubMatA%RInd(:)   )
          s1A_Im(1:ncs) = s1pl_tmp( indMoI(i,k)%SubMatA%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatA%CC(  :)
          u2A_Re(1:ncs) = u2pl_tmp( indMoI(i,k)%SubMatA%RInd(:)   )
          u2A_Im(1:ncs) = u2pl_tmp( indMoI(i,k)%SubMatA%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatA%CC(  :)
          s2A_Re(1:ncs) = s2pl_tmp( indMoI(i,k)%SubMatA%RInd(:)   )
          s2A_Im(1:ncs) = s2pl_tmp( indMoI(i,k)%SubMatA%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatA%CC(  :)
          u3A_Re(1:ncs) = u3pl_tmp( indMoI(i,k)%SubMatA%RInd(:)   )
          u3A_Im(1:ncs) = u3pl_tmp( indMoI(i,k)%SubMatA%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatA%CC(  :)
          s3A_Re(1:ncs) = s3pl_tmp( indMoI(i,k)%SubMatA%RInd(:)   )
          s3A_Im(1:ncs) = s3pl_tmp( indMoI(i,k)%SubMatA%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatA%CC(  :)
      
          u1B_Re(1:ncs) = u1pl_tmp( indMoI(i,k)%SubMatB%RInd(:)   )
          u1B_Im(1:ncs) = u1pl_tmp( indMoI(i,k)%SubMatB%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatB%CC(  :)
          u2B_Re(1:ncs) = u2pl_tmp( indMoI(i,k)%SubMatB%RInd(:)   )
          u2B_Im(1:ncs) = u2pl_tmp( indMoI(i,k)%SubMatB%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatB%CC(  :)
          s2B_Re(1:ncs) = s2pl_tmp( indMoI(i,k)%SubMatB%RInd(:)   )
          s2B_Im(1:ncs) = s2pl_tmp( indMoI(i,k)%SubMatB%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatB%CC(  :)
          u3B_Re(1:ncs) = u3pl_tmp( indMoI(i,k)%SubMatB%RInd(:)   )
          u3B_Im(1:ncs) = u3pl_tmp( indMoI(i,k)%SubMatB%RInd(:)+1 ) &
                      & *           indMoI(i,k)%SubMatB%CC(  :)

          allocate(buff_EP(3,2,ncs), buff_Re(3,ncs), buff_Im(3,ncs))


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate advection term for uu, uv, uw !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! life would be better if people labelled things consistently- uu, vu, wu in the eqn -__- 
          buff_Re(1,1:ncs) = ( u1B_Im(1:ncs)*u1A_Re(1:ncs) + u1B_Re(1:ncs)*u1A_Im(1:ncs) ) * kc_x
          buff_Im(1,1:ncs) = ( u1B_Im(1:ncs)*u1A_Im(1:ncs) - u1B_Re(1:ncs)*u1A_Re(1:ncs) ) * kc_x
          buff_Re(2,1:ncs) = - u2B_Re(1:ncs)*s1A_Re(1:ncs) + u2B_Im(1:ncs)*s1A_Im(1:ncs) &
                              & - s2B_Re(1:ncs)*u1A_Re(1:ncs) + s2B_Im(1:ncs)*u1A_Im(1:ncs)
          buff_Im(2,1:ncs) = - u2B_Re(1:ncs)*s1A_Im(1:ncs) - u2B_Im(1:ncs)*s1A_Re(1:ncs) &
                              & - s2B_Re(1:ncs)*u1A_Im(1:ncs) - s2B_Im(1:ncs)*u1A_Re(1:ncs)
          buff_Re(3,1:ncs) = ( u3B_Im(1:ncs)*u1A_Re(1:ncs) + u3B_Re(1:ncs)*u1A_Im(1:ncs) ) * kc_z
          buff_Im(3,1:ncs) = ( u3B_Im(1:ncs)*u1A_Im(1:ncs) - u3B_Re(1:ncs)*u1A_Re(1:ncs) ) * kc_z



          buff_EP(:,1,1:ncs) = buff_Re(:,1:ncs)*u1C_Re                           + buff_Im(:,1:ncs)*u1C_Im
          buff_EP(:,2,1:ncs) = buff_Re(:,1:ncs)*u1C_Re/sqrt(u1C_Re**2+u1C_Im**2) + &
                              & buff_Im(:,1:ncs)*u1C_Im/sqrt(u1C_Re**2+u1C_Im**2)
      

          buff_fold = 0

          do icsub = 1, ncs
              iia = indMoI(i,k)%SubMatA%ni(  icsub) 
              kka = indMoI(i,k)%SubMatA%nk(  icsub) 
              kka = kka*(1-2*IkkcNeg)
              IkkNeg = (1 - max(isign(1,kka),0))*2
              IiiNeg = 1 - max(isign(1,iia),0)

              sigCase = 1+IkkNeg + IiiNeg ! 1: (+,+), 2: (-,+), 3: (+,-), 4: (-,-)
              buff_fold(:,1,sigCase, abs(iia),abs(kka)) = buff_fold(:,1,sigCase, abs(iia),abs(kka)) + buff_EP(:,1, icsub)
              buff_fold(:,3,sigCase, abs(iia),abs(kka)) = buff_fold(:,3,sigCase, abs(iia),abs(kka)) + buff_EP(:,2, icsub)
              ! if (mod(icsub,12)==0) then
              !     write(*,*) iia, kka, sigCase
              ! end if



              iib = indMoI(i,k)%SubMatB%ni(  icsub) 
              kkb = indMoI(i,k)%SubMatB%nk(  icsub) 
              kkb = kkb*(1-2*IkkcNeg)
              IkkNeg = (1 - max(isign(1,kkb),0))*2
              IiiNeg = 1 - max(isign(1,iib),0)
              sigCase = 1+IkkNeg + IiiNeg
              buff_fold(:,2,sigCase, abs(iib),abs(kkb)) = buff_fold(:,2,sigCase, abs(iib),abs(kkb)) + buff_EP(:,1, icsub)
              buff_fold(:,4,sigCase, abs(iib),abs(kkb)) = buff_fold(:,4,sigCase, abs(iib),abs(kkb)) + buff_EP(:,2, icsub)       
          end do
          buff_fold(:,:,1,NXoI(i), NZoI(k_ind)) = 0.0d0

          do ii = 1, MoINumX
              do kk = 1, MoINumZ
                  ind = (kk-1)*MoINumX + ii
                  do sigCase = 1,4
                      NxUpp = merge(xlim(sigCase), NXlim(ii,2), xlim(sigCase) >= NXlim(ii,1) .and. xlim(sigCase) <= NXlim(ii,2))
                      NzUpp = merge(zlim(sigCase), NZlim(kk,2), zlim(sigCase) >= NZlim(kk,1) .and. zlim(sigCase) <= NZlim(kk,2))
                      buffInt = merge(1, 0, NXoI(i) >= NXlim(ii,1) .and. NXoI(i) <= NxUpp .and. NZoI(k_ind) >= NZlim(kk,1) .and. NZoI(k_ind) <= NzUpp)
                      buffInt = max(((NxUpp-NXlim(ii,1)+1)*(NzUpp-NZlim(kk,1)+1)-buffInt),1)
                      buff_fib(:,1:4, sigCase,ind) = sum( sum(buff_fold(:,1:4, sigCase,NXlim(ii,1):NxUpp, NZlim(kk,1):NzUpp), dim=4), dim=3 ) / buffInt
                  end do
              end do
          end do
          convs_uu_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_uu_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(1, 1:4,1:4, :)
          convs_uv_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_uv_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(2, 1:4,1:4, :) 
          convs_uw_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_uw_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(3, 1:4,1:4, :) 


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate advection term for vu, vv, vw !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          buff_Re(1,1:ncs) = ( u1B_Im(1:ncs)*u2A_Re(1:ncs) + u1B_Re(1:ncs)*u2A_Im(1:ncs) ) * kc_x
          buff_Im(1,1:ncs) = ( u1B_Im(1:ncs)*u2A_Im(1:ncs) - u1B_Re(1:ncs)*u2A_Re(1:ncs) ) * kc_x
          buff_Re(2,1:ncs) = - u2B_Re(1:ncs)*s2A_Re(1:ncs) + u2B_Im(1:ncs)*s2A_Im(1:ncs) &
                                  & - s2B_Re(1:ncs)*u2A_Re(1:ncs) + s2B_Im(1:ncs)*u2A_Im(1:ncs)
          buff_Im(2,1:ncs) = - u2B_Re(1:ncs)*s2A_Im(1:ncs) - u2B_Im(1:ncs)*s2A_Re(1:ncs) &
                                  & - s2B_Re(1:ncs)*u2A_Im(1:ncs) - s2B_Im(1:ncs)*u2A_Re(1:ncs)
          buff_Re(3,1:ncs) = ( u3B_Im(1:ncs)*u2A_Re(1:ncs) + u3B_Re(1:ncs)*u2A_Im(1:ncs) ) * kc_z
          buff_Im(3,1:ncs) = ( u3B_Im(1:ncs)*u2A_Im(1:ncs) - u3B_Re(1:ncs)*u2A_Re(1:ncs) ) * kc_z


          buff_EP(:,1,1:ncs) = buff_Re(:,1:ncs)*u2C_Re                           + buff_Im(:,1:ncs)*u2C_Im
          buff_EP(:,2,1:ncs) = buff_Re(:,1:ncs)*u2C_Re/sqrt(u2C_Re**2+u2C_Im**2) &
                            & + buff_Im(:,1:ncs)*u2C_Im/sqrt(u2C_Re**2+u2C_Im**2)

          buff_fold = 0
          ! 1: (+,+), 2: (-,+), 3: (+,-), 4: (-,-)
          do icsub = 1, ncs
              iia = indMoI(i,k)%SubMatA%ni(  icsub) 
              kka = indMoI(i,k)%SubMatA%nk(  icsub) 
              kka = kka*(1-2*IkkcNeg)
              IkkNeg = (1 - max(isign(1,kka),0))*2
              IiiNeg = 1 - max(isign(1,iia),0)
              sigCase = 1+IkkNeg + IiiNeg ! 1: (+,+), 2: (-,+), 3: (+,-), 4: (-,-)
              buff_fold(:,1,sigCase, abs(iia),abs(kka)) = buff_fold(:,1,sigCase, abs(iia),abs(kka)) + buff_EP(:,1, icsub)
              buff_fold(:,3,sigCase, abs(iia),abs(kka)) = buff_fold(:,3,sigCase, abs(iia),abs(kka)) + buff_EP(:,2, icsub)

              iib = indMoI(i,k)%SubMatB%ni(  icsub) 
              kkb = indMoI(i,k)%SubMatB%nk(  icsub)
              kkb = kkb*(1-2*IkkcNeg) 
              IkkNeg = (1 - max(isign(1,kkb),0))*2
              IiiNeg = 1 - max(isign(1,iib),0)
              sigCase = 1+IkkNeg + IiiNeg
              buff_fold(:,2,sigCase, abs(iib),abs(kkb)) = buff_fold(:,2,sigCase, abs(iib),abs(kkb)) + buff_EP(:,1, icsub)
              buff_fold(:,4,sigCase, abs(iib),abs(kkb)) = buff_fold(:,4,sigCase, abs(iib),abs(kkb)) + buff_EP(:,2, icsub)
          end do

          buff_fold(:,:,1,NXoI(i), NZoI(k_ind)) = 0.0d0

          do ii = 1, MoINumX
              do kk = 1, MoINumZ
                  ind = (kk-1)*MoINumX + ii
                  do sigCase = 1,4
                      !Using the fibonacci again
                      NxUpp = merge(xlim(sigCase), NXlim(ii,2), xlim(sigCase) >= NXlim(ii,1) .and. xlim(sigCase) <= NXlim(ii,2))
                      NzUpp = merge(zlim(sigCase), NZlim(kk,2), zlim(sigCase) >= NZlim(kk,1) .and. zlim(sigCase) <= NZlim(kk,2))
                      buffInt = merge(1, 0, NXoI(i) >= NXlim(ii,1) .and. NXoI(i) <= NxUpp .and. NZoI(k_ind) >= NZlim(kk,1) .and. NZoI(k_ind) <= NzUpp)
                      buffInt = max(((NxUpp-NXlim(ii,1)+1)*(NzUpp-NZlim(kk,1)+1)-buffInt),1)
                      !Averaging buff_fold within the fibonacci regions from before, and storing in buff_fib
                      buff_fib(:,1:4, sigCase,ind) = sum( sum(buff_fold(:,1:4, sigCase,NXlim(ii,1):NxUpp, NZlim(kk,1):NzUpp), dim=4), dim=3 ) / buffInt
                  end do
              end do
          end do

          convs_vu_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_vu_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(1, 1:4,1:4, :)
          convs_vv_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_vv_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(2, 1:4,1:4, :) 
          convs_vw_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_vw_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(3, 1:4,1:4, :)    

          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculate advection term for wu, wv, ww !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          buff_Re(1,1:ncs) = ( u1B_Im(1:ncs)*u3A_Re(1:ncs) + u1B_Re(1:ncs)*u3A_Im(1:ncs) ) * kc_x  !this is the -ikx*uu
          buff_Im(1,1:ncs) = ( u1B_Im(1:ncs)*u3A_Im(1:ncs) - u1B_Re(1:ncs)*u3A_Re(1:ncs) ) * kc_x
          buff_Re(2,1:ncs) = - u2B_Re(1:ncs)*s3A_Re(1:ncs) + u2B_Im(1:ncs)*s3A_Im(1:ncs) &
                                  & - s2B_Re(1:ncs)*u3A_Re(1:ncs) + s2B_Im(1:ncs)*u3A_Im(1:ncs)
          buff_Im(2,1:ncs) = - u2B_Re(1:ncs)*s3A_Im(1:ncs) - u2B_Im(1:ncs)*s3A_Re(1:ncs) &
                                  & - s2B_Re(1:ncs)*u3A_Im(1:ncs) - s2B_Im(1:ncs)*u3A_Re(1:ncs)
          buff_Re(3,1:ncs) = ( u3B_Im(1:ncs)*u3A_Re(1:ncs) + u3B_Re(1:ncs)*u3A_Im(1:ncs) ) * kc_z
          buff_Im(3,1:ncs) = ( u3B_Im(1:ncs)*u3A_Im(1:ncs) - u3B_Re(1:ncs)*u3A_Re(1:ncs) ) * kc_z


          buff_EP(:,1,1:ncs) = buff_Re(:,1:ncs)*u3C_Re                           + buff_Im(:,1:ncs)*u3C_Im
          buff_EP(:,2,1:ncs) = buff_Re(:,1:ncs)*u3C_Re/sqrt(u3C_Re**2+u3C_Im**2) & 
                            & + buff_Im(:,1:ncs)*u3C_Im/sqrt(u3C_Re**2+u3C_Im**2)
          buff_fold = 0
          do icsub = 1, ncs
              iia = indMoI(i,k)%SubMatA%ni(  icsub) 
              kka = indMoI(i,k)%SubMatA%nk(  icsub) 
              kka = kka*(1-2*IkkcNeg)
              IkkNeg = (1 - max(isign(1,kka),0))*2
              IiiNeg = 1 - max(isign(1,iia),0)
              sigCase = 1+IkkNeg + IiiNeg ! 1: (+,+), 2: (-,+), 3: (+,-), 4: (-,-)
              buff_fold(:,1,sigCase, abs(iia),abs(kka)) = buff_fold(:,1,sigCase, abs(iia),abs(kka)) + buff_EP(:,1, icsub)
              buff_fold(:,3,sigCase, abs(iia),abs(kka)) = buff_fold(:,3,sigCase, abs(iia),abs(kka)) + buff_EP(:,2, icsub)

              iib = indMoI(i,k)%SubMatB%ni(  icsub) 
              kkb = indMoI(i,k)%SubMatB%nk(  icsub) 
              kkb = kkb*(1-2*IkkcNeg)
              IkkNeg = (1 - max(isign(1,kkb),0))*2
              IiiNeg = 1 - max(isign(1,iib),0)
              sigCase = 1+IkkNeg + IiiNeg
              buff_fold(:,2,sigCase, abs(iib),abs(kkb)) = buff_fold(:,2,sigCase, abs(iib),abs(kkb)) + buff_EP(:,1, icsub)
              buff_fold(:,4,sigCase, abs(iib),abs(kkb)) = buff_fold(:,4,sigCase, abs(iib),abs(kkb)) + buff_EP(:,2, icsub)
          end do
      
          buff_fold(:,:,1,NXoI(i), NZoI(k_ind)) = 0.0d0

          do ii = 1, MoINumX
              do kk = 1, MoINumZ
                  ind = (kk-1)*MoINumX + ii
                  do sigCase = 1,4
                      NxUpp = merge(xlim(sigCase), NXlim(ii,2), xlim(sigCase) >= NXlim(ii,1) .and. xlim(sigCase) <= NXlim(ii,2))
                      NzUpp = merge(zlim(sigCase), NZlim(kk,2), zlim(sigCase) >= NZlim(kk,1) .and. zlim(sigCase) <= NZlim(kk,2))
                      buffInt = merge(1, 0, NXoI(i) >= NXlim(ii,1) .and. NXoI(i) <= NxUpp .and. NZoI(k_ind) >= NZlim(kk,1) .and. NZoI(k_ind) <= NzUpp)
                      buffInt = max(((NxUpp-NXlim(ii,1)+1)*(NzUpp-NZlim(kk,1)+1)-buffInt),1)
                      buff_fib(:,1:4, sigCase,ind) = sum( sum(buff_fold(:,1:4, sigCase,NXlim(ii,1):NxUpp, NZlim(kk,1):NzUpp), dim=4), dim=3 ) / buffInt
                  end do
              end do
          end do

          convs_wu_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_wu_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(1, 1:4,1:4, :)
          convs_wv_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_wv_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(2, 1:4,1:4, :) 
          convs_ww_tmp(1:4,1:4, jpl, i, k_ind, :) = convs_ww_tmp(1:4,1:4, jpl, i, k_ind, :) + buff_fib(3, 1:4,1:4, :)
          deallocate(buff_EP, buff_Re, buff_Im)

          


        end do
      end do
    end do
  End Do

  ! ! MPI reduce to combine upper and lower channel data!!
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_uu_tmp, size(convs_uu_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_uv_tmp, size(convs_uv_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_uw_tmp, size(convs_uw_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_vu_tmp, size(convs_vu_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_vv_tmp, size(convs_vv_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_vw_tmp, size(convs_vw_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_wu_tmp, size(convs_wu_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_wv_tmp, size(convs_wv_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, convs_ww_tmp, size(convs_ww_tmp), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  planeSize = 4*4*MoINumX*MoINumZ*MoINumt
  allocate(planeBuf(4,4,MoINumX,MoINumZ,MoINumt))
  allocate(recvPlane(4,4,MoINumX,MoINumZ,MoINumt))

  ! -------------------------
  ! lower owners send only their own planes
  ! -------------------------
  do idx = 1, n_planesL
    jpl = jpl_listL(idx)

    if (ownerL(jpl) /= ownerU(jpl)) then
      dst = ownerU(jpl)

      planeBuf = convs_uu_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 1, MPI_COMM_WORLD, ierr)

      planeBuf = convs_uv_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 2, MPI_COMM_WORLD, ierr)

      planeBuf = convs_uw_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 3, MPI_COMM_WORLD, ierr)

      planeBuf = convs_vu_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 4, MPI_COMM_WORLD, ierr)

      planeBuf = convs_vv_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 5, MPI_COMM_WORLD, ierr)

      planeBuf = convs_vw_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 6, MPI_COMM_WORLD, ierr)

      planeBuf = convs_wu_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 7, MPI_COMM_WORLD, ierr)

      planeBuf = convs_wv_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 8, MPI_COMM_WORLD, ierr)

      planeBuf = convs_ww_tmp(:,:,jpl,:,:,:)
      call MPI_SEND(planeBuf, planeSize, MPI_REAL8, dst, 1000 + 10*jpl + 9, MPI_COMM_WORLD, ierr)
    end if
  end do


  ! -------------------------
  ! upper owners receive only their own planes
  ! -------------------------
  do idx = 1, n_planesU
    jpl = jpl_listU(idx)

    if (ownerL(jpl) /= ownerU(jpl)) then
      src = ownerL(jpl)

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 1, MPI_COMM_WORLD, status, ierr)
      convs_uu_tmp(:,:,jpl,:,:,:) = convs_uu_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 2, MPI_COMM_WORLD, status, ierr)
      convs_uv_tmp(:,:,jpl,:,:,:) = convs_uv_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 3, MPI_COMM_WORLD, status, ierr)
      convs_uw_tmp(:,:,jpl,:,:,:) = convs_uw_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 4, MPI_COMM_WORLD, status, ierr)
      convs_vu_tmp(:,:,jpl,:,:,:) = convs_vu_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 5, MPI_COMM_WORLD, status, ierr)
      convs_vv_tmp(:,:,jpl,:,:,:) = convs_vv_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 6, MPI_COMM_WORLD, status, ierr)
      convs_vw_tmp(:,:,jpl,:,:,:) = convs_vw_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 7, MPI_COMM_WORLD, status, ierr)
      convs_wu_tmp(:,:,jpl,:,:,:) = convs_wu_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 8, MPI_COMM_WORLD, status, ierr)
      convs_wv_tmp(:,:,jpl,:,:,:) = convs_wv_tmp(:,:,jpl,:,:,:) + recvPlane

      call MPI_RECV(recvPlane, planeSize, MPI_REAL8, src, 1000 + 10*jpl + 9, MPI_COMM_WORLD, status, ierr)
      convs_ww_tmp(:,:,jpl,:,:,:) = convs_ww_tmp(:,:,jpl,:,:,:) + recvPlane

    end if
  end do


  convs_uu = convs_uu + convs_uu_tmp
  convs_uv = convs_uv + convs_uv_tmp
  convs_uw = convs_uw + convs_uw_tmp

  convs_vu = convs_vu + convs_vu_tmp
  convs_vv = convs_vv + convs_vv_tmp
  convs_vw = convs_vw + convs_vw_tmp

  convs_wu = convs_wu + convs_wu_tmp
  convs_wv = convs_wv + convs_wv_tmp
  convs_ww = convs_ww + convs_ww_tmp

  if (myid == 7) then
    write(6,*) "convs_uw", convs_uw(1,1,1, 5, 10, 1:100)
  end if 

  deallocate(convs_uu_tmp, convs_uv_tmp, convs_uw_tmp)
  deallocate(convs_vu_tmp, convs_vv_tmp, convs_vw_tmp)
  deallocate(convs_wu_tmp, convs_wv_tmp, convs_ww_tmp)

  deallocate(u1_ind, u2_ind, u3_ind)
  deallocate(u1PL_map, u2PL_map, u3PL_map)


    

end subroutine 

subroutine write_map(myid)
  
  use declaration
  implicit none
  include 'mpif.h'

  integer :: myid, Ntx, Ntz, ierr
  integer :: idx
  integer :: jpl,j,i,k,f,s

  character(len=3)  :: extnkx, extnkz, extny, extmod
  character(len=256):: output, map_output
  character(len=256):: fnameList
  
  Ntx = Nspec_x
  Ntz = Nspec_z


  write(*,*) ''

  ! print *, 'Begin writing'

  ! do jpl = 1,PLoINum

  do idx = 1, n_planesU
    jpl = jpl_listU(idx)
    j   = NYoI(jpl)
    write(extnkx,'(i3.3)') MoINumX
    write(extnkz,'(i3.3)') MoINumZ
    write(extny ,'(i3.3)') int(yPLoi(jpl))

    output = 'output'
    map_output = 'map_output'

    write(ext4,'(i5.5)') int(100d0*(t))!int(t)!

    fnameList = trim(output)//'/'//trim(map_output)//'/TRD_'//extnkx//'_'//extnkz//'_'//extny//'_t'//ext4//'.dat'

    open(unit=50, file=fnameList, form='unformatted')
    write(50) Re, alp, bet, mpgx, Ntx, Ntz
    write(50) NYoI(jpl), yPLoi(jpl)
    write(50) NXfib(:), NZfib(:)
    write(50) NXoI(:), NZoI(:)
    write(50) nsamp
    do i = 1,MoINumX
      do k = 1, MoINumZ
        write(50)
        write(50) NXoI(i), NZoI(k)
        do f = 1, 4
          do s = 1,4
          write(50) convs_uu(f, s, jpl, i, k, :), convs_uv(f, s, jpl, i, k, :), convs_uw(f, s, jpl, i, k, :)
          write(50) convs_vu(f, s, jpl, i, k, :), convs_vv(f, s, jpl, i, k, :), convs_vw(f, s, jpl, i, k, :)
          write(50) convs_wu(f, s, jpl, i, k, :), convs_wv(f, s, jpl, i, k, :), convs_ww(f, s, jpl, i, k, :)

              ! if (i==7 .and. k==7) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !     write(*,*) convs_uu(1, 1, jpl, i, k, 50:52)
              ! end if
          end do
        end do
      end do
    end do
    close(50)
  end do

  convs_uu = 0.0d0
  convs_uv = 0.0d0
  convs_uw = 0.0d0
  convs_vu = 0.0d0
  convs_vv = 0.0d0
  convs_vw = 0.0d0
  convs_wu = 0.0d0
  convs_wv = 0.0d0
  convs_ww = 0.0d0
  nsamp = 0

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

end subroutine