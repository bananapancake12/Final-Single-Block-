module declaration

  implicit none

  real(8) :: pi
  real(8) :: t, Re, alp, bet, mpgx

  integer :: nx, nz, ny
  integer :: Ntx, Ntz
  integer :: firstgood
  integer :: nsamp
  integer :: nsamp_total

  integer :: MoINumX, MoINumZ, MoINumt
  integer :: PLoINum

  integer :: NYoI_jpl
  real(8) :: yPLoI_jpl

  character(len=80)  :: filelist
  character(len=100) :: fildir
  character(len=100) :: filout
  character(len=100) :: fnameima

  character(len=5), allocatable :: extL(:)

  character(len=3) :: extnkx
  character(len=3) :: extnkz
  character(len=3) :: extny
  character(len=5) :: ext4

  integer, allocatable :: NYoI(:)
  real(8), allocatable :: yPLoI(:)

  real(8), allocatable :: NXfib(:)
  real(8), allocatable :: NZfib(:)
  integer, allocatable :: NXoI(:)
  integer, allocatable :: NZoI(:)

  real(8), allocatable :: convs_uu(:,:,:,:,:)
  real(8), allocatable :: convs_uv(:,:,:,:,:)
  real(8), allocatable :: convs_uw(:,:,:,:,:)

  real(8), allocatable :: convs_vu(:,:,:,:,:)
  real(8), allocatable :: convs_vv(:,:,:,:,:)
  real(8), allocatable :: convs_vw(:,:,:,:,:)

  real(8), allocatable :: convs_wu(:,:,:,:,:)
  real(8), allocatable :: convs_wv(:,:,:,:,:)
  real(8), allocatable :: convs_ww(:,:,:,:,:)

  real(8), allocatable :: buu(:), buv(:), buw(:)
  real(8), allocatable :: bvu(:), bvv(:), bvw(:)
  real(8), allocatable :: bwu(:), bwv(:), bww(:)

end module