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

  real(8) :: t1

    
  integer physlim_bot
  integer physlim_top
  
  integer nn
  integer iter,iter0,nwrite,iwrite !itersl,nstatsl
  integer nstat,istat
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

end module
