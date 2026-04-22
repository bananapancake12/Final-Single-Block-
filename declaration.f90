module declaration

  ! Static variables
  
  real(8), parameter :: pi = 4d0*datan(1d0)
  complex, parameter :: im = dcmplx(0d0,1d0)
  
  integer, parameter :: nband = 3
  
  integer, parameter :: ugrid = 2 ! centres     including ghost points: u and w
  integer, parameter :: vgrid = 1 ! faces: v
  integer, parameter :: pgrid = 3 ! centres not including ghost points: p


  ! NEW AND CHANGED VARIABLES (no more N and Ngal)
  
  integer :: Ngal_x, Ngal_z
  integer :: Nspec_x, Nspec_z
  integer :: nyu, nyv, nyp, nyu_LB, nyv_LB, nyp_LB

  
  ! All variables

  integer bandit(3)
  integer np,pnodes
  integer,allocatable:: procs(:)
  integer,allocatable:: N(:,:),Ngal(:,:),Ny(:,:)
  integer,allocatable:: planelim(:,:,:)
  integer,allocatable:: limPL_incw(:,:,:),limPL_excw(:,:,:)
  integer,allocatable:: limPL_FFT(:,:,:)
  ! integer,allocatable:: bandPL(:)
  ! integer,allocatable:: bandPL_FFT(:)
  integer            :: jgal(3,2),igal,kgal 
  integer, allocatable :: columns_num(:)
  integer, allocatable :: columns_i(:,:)
  integer, allocatable :: columns_k(:,:)
  integer, allocatable :: jlim(:,:)
  integer, allocatable :: dk(:,:)
  integer, allocatable :: dk_phys(:,:)

  integer(8), allocatable :: weight(:)
  integer recv_flg

  real(8) :: t1

    
  integer physlim_bot
  integer physlim_top
  
  integer nn
  integer iter,iter0,nwrite,iwrite, nsamp,iter0mp, iwrite_map, nwrite_map
  integer nstat,istat, nmap
  integer flag_init,flag_ctpress
  real(8) nextqt
  integer geometry_type
  integer kRK
  real(8) t,dt,CFL,maxt,dtv,dtc,dti,Re
  real(8) gridweighting_bc_u1,gridweighting_bc_u3
  real(8) alp,bet             ! Wavelengths
  real(8) Lx,Ly,Lz            ! Size of the computational box
  real(8) Kib

  ! Grid
  real(8),pointer:: yu(:),dyu2i(:,:),dthdyu(:)
  real(8),pointer:: yv(:),dyv2i(:,:),dthdyv(:)
  real(8) dtheta,dthetai
  real(8) dthetavi,ddthetavi
  real(8),pointer:: gridweighting(:)
  real(8),pointer:: gridweighting_interp(:)
  integer ppp
  real(8) dyq

  ! Variables in planes
  real(8), allocatable :: u1PL(:,:,:), u2PL(:,:,:), u3PL(:,:,:)
  real(8), allocatable :: u1PL_itp(:,:,:), u2PL_itp(:,:,:), u3PL_itp(:,:,:)
  real(8), allocatable :: Nu1PL(:,:,:), Nu2PL(:,:,:), Nu3PL(:,:,:)
  ! real(8), allocatable :: du1PL(:,:,:), du2PL(:,:,:), du3PL(:,:,:)
  real(8), allocatable :: wx(:,:,:), ppPL(:,:,:)
  real(8), allocatable :: Qcrit(:,:,:)
  real(8), allocatable :: u1PLN(:,:,:), u2PLN(:,:,:), u3PLN(:,:,:), ppPLN(:,:,:)
  real(8), allocatable :: u1PL_itpN(:,:,:), u2PL_itpN(:,:,:), u3PL_itpN(:,:,:)

  ! Cross products in planes
  real(8), allocatable :: uu_cPL(:,:,:), uv_fPL(:,:,:), uw_cPL(:,:,:)
  real(8), allocatable :: vu_fPL(:,:,:), vv_cPL(:,:,:), vw_fPL(:,:,:)
  real(8), allocatable :: wu_cPL(:,:,:), wv_fPL(:,:,:), ww_cPL(:,:,:)

  real(8), allocatable :: Nu1PL_dy(:,:,:), Nu2PL_dy(:,:,:), Nu3PL_dy(:,:,:)


  ! Spectra
  real(8),pointer::  spU(:,:), spV(:,:), spW(:,:)
  real(8),pointer:: spUV(:,:), spP(:,:)

  ! Statistics. Mean statistics
  real(8), allocatable :: Um(:), U2m(:)
  real(8), allocatable :: Vm(:), V2m(:)
  real(8), allocatable :: Wm(:), W2m(:)
  real(8), allocatable :: Pm(:), P2m(:)

  real(8), allocatable :: UVm(:)
  real(8), allocatable :: UWm(:)
  real(8), allocatable :: VWm(:)

  real(8), allocatable :: wxm(:), wx2m(:)


  real(8),pointer:: u11(:)
  
  complex(8),allocatable:: k1F_x(:),k1F_z(:)
  real(8)   ,allocatable:: k2F_x(:),k2F_z(:)

  ! Runge-Kutta coefficients
  real(8) aRK(3),bRK(3),gRK(3),cRK(3),dRK(3)

  real(8) err,maxerr,maxA!,lambda,lambdaQ,Re_div,iRediv
  real(8) mpgx,mpgz,dgx,dgz,QxT,Qx,Qz,Umax,utau
  character*120 fnameima,fnameimb,fnameimc,boundfname,filout,directory
  logical exist_file_hist
  character*120 dirin,dirout, dirlist, heading
  character*4 ext1,ext2,ext3
  character*5 ext4


  real(8), allocatable :: buffR_x(:), buffC_z(:), buffRal_x(:), buffCal_z(:)

  complex(8), allocatable :: u1_itp(:,:),u2_itp(:,:),u3_itp(:,:)
  complex(8), allocatable :: Nu1_dy(:,:),Nu2_dy(:,:),Nu3_dy(:,:)
  complex(8), allocatable :: uv_f(:,:), wv_f(:,:), vv_c(:,:)

  
  ! Omega x
  real(8),      allocatable :: du1dy_planes(:,:,:)
  real(8),      allocatable :: du2dy_planes(:,:,:)
  real(8),      allocatable :: du3dy_planes(:,:,:)
  
  real(8),      allocatable :: du1dy_planes2(:,:,:)
  real(8),      allocatable :: du2dy_planes2(:,:,:)
  real(8),      allocatable :: du3dy_planes2(:,:,:)
  
  complex(8), allocatable :: du1dy_columns(:,:), du2dy_columns(:,:), du3dy_columns(:,:)
  
  real(8), allocatable :: DG(:,:,:)
  
  !real(8) bslip
  
  ! for nonlinear interaction list (added by JC)
  integer, parameter :: int1 = selected_int_kind(2)  ! at least 2 decimal digits → 1 byte
  type :: nonlinList
    integer(kind=int1), allocatable :: list(:,:)
  end type nonlinList
  type(nonlinList), allocatable :: nonlin(:,:)

  integer, allocatable:: iLkup(:), kLkup(:), iNeg(:)



  !!!!!!!!!!!!!!      TRD_v2 New variables     !!!!!!!!!!!!!!!!!
  integer :: nxf, nzf, nyf                 
  integer :: nx_trd, nz_trd, nnc                   
  integer :: nfib                          

  integer :: MoINumX, MoINumZ, MoINumt     ! sizes

  integer, allocatable :: NXlim(:,:), NZlim(:,:)   
  real(8), allocatable :: NXfib(:), NZfib(:)       
  integer, allocatable :: NXoI(:),  NZoI(:)        

  integer :: PLoINum, PLoINumEx, inoutunit
  integer, allocatable :: NYoI(:)                 
  real(8), allocatable :: yPLoi(:)               
  integer, allocatable :: jList_Buff(:)     
  integer, allocatable :: buffIndj(:)   
  real(8) :: Re_tau


 type PLj
      integer :: j
      type(PLj), pointer :: next
  end type PLj

  type CtrMod
    integer :: RIndC
    integer :: nsub
    integer, allocatable :: ni(:), nk(:)
    integer, allocatable :: CC(:)
    integer, allocatable :: RInd(:)
  end type CtrMod

  type Mxz
    integer :: NXoI, NZoI
    type(CtrMod), allocatable :: SubMatA, SubMatB
  end type Mxz

  type(Mxz), allocatable :: indMoI(:,:)   ! (MoINumX, 2*MoINumZ)
  !real(8), allocatable :: u1_ind(:,:), u2_ind(:,:), u3_ind(:,:)


  ! === TRD work arrays ===
  real(8), allocatable :: RBf_u1(:,:,:), RBf_u2(:,:,:), RBf_u3(:,:,:), RBf_pr(:,:,:)
  real(8), allocatable :: u1pl_tmp(:), u2pl_tmp(:), u3pl_tmp(:)
  real(8), allocatable :: s1pl_tmp(:), s2pl_tmp(:), s3pl_tmp(:)


  real(8), allocatable :: u1A_Re(:), u1A_Im(:), s1A_Re(:), s1A_Im(:)
  real(8), allocatable :: u2A_Re(:), u2A_Im(:), s2A_Re(:), s2A_Im(:)
  real(8), allocatable :: u3A_Re(:), u3A_Im(:), s3A_Re(:), s3A_Im(:)

  real(8), allocatable :: u1B_Re(:), u1B_Im(:)
  real(8), allocatable :: u2B_Re(:), u2B_Im(:), s2B_Re(:), s2B_Im(:)
  real(8), allocatable :: u3B_Re(:), u3B_Im(:)

  real(8), allocatable :: ka_x(:), ka_z(:)

  real(8), allocatable :: convs_uu(:,:,:,:,:,:), convs_uv(:,:,:,:,:,:), convs_uw(:,:,:,:,:,:)
  real(8), allocatable :: convs_vu(:,:,:,:,:,:), convs_vv(:,:,:,:,:,:), convs_vw(:,:,:,:,:,:)
  real(8), allocatable :: convs_wu(:,:,:,:,:,:), convs_wv(:,:,:,:,:,:), convs_ww(:,:,:,:,:,:)

  real(8), allocatable :: buff_fold(:,:,:,:,:), buff_fib(:,:,:,:)

  real(8), allocatable :: buff_Re(:,:), buff_Im(:,:), buff_EP(:,:,:)

  integer :: n_planesU
  integer, allocatable :: jpl_listU(:)





end module
