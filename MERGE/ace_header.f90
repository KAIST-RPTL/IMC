module ace_header

implicit none
character(100):: library_path    !> Path of ACE library


type TabularDataForm
  integer :: NR        !> number of interpolation regions
  real(8), allocatable :: NBT(:) !> ENDF interpolation parameters, NBT(I), I=1,NR
  real(8), allocatable :: INT(:) !> ENDF interpolation parameters, INT(I), I=1,NR
  integer :: NE        !> number of energies
  real(8), allocatable :: E(:)   !> number of energies, E(I), I=1,NE
  real(8), allocatable :: F(:)   !> energy dependent data, F(I), 1=1,NE
end type 



type, extends (TabularDataForm) :: NuTotDataForm 
  integer :: NC                !> number of coefficients
  real(8), allocatable :: C(:) !> Coefficients C(I), I=1,NC
end type


type CrossSectionDataForm
  integer :: IE        !> first energy grid index correspond to E(:) in ESZ Block
  integer :: NE        !> number of consecutive entries
  real(8), allocatable :: cx(:) !> cross sections, sig(I), I=IE,IE+NE-1
end type



type AngularDistDataForm
  !> case( dist_flag(IE) = 0 ) :: isotropic distribution, no data is needed. 
  !> case( dist_flag(IE) < 0 ) :: tabular angular distribution
  !>   LDAT( 1 : 2 + 3*NP )
  !>   LDAT(1) = JJ, interpolation flag
  !>     LDAT(1)=1 :: histogram
  !>     LDAT(1)=2 :: lin-lin
  !>   LDAT(2) = NP, number of points in the distribution
  !>   LDAT( 3      : 3+NP-1   ) = CSOUT(1:NP), cosine scattering angular grid
  !>   LDAT( 3+NP   : 3+2*NP-1 ) = PDF(1:NP), probability density function 
  !>   LDAT( 3+2*NP : 3+3*NP-1 ) = CDF(1:NP), cumulative density function
  !> case( dist_flag(IE) > 0 ) :: 32 equiprobable cosine bins for scattering 
  !real(8), allocatable :: LDAT(:) 
  
  !> 32 equiprobable distribution 
  !real(8), allocatable :: P(:) 
  
  !> tabular angular distribution
  !integer :: JJ = 0, NP = 0
  !real(8), allocatable :: CSOUT(:) 
  !real(8), allocatable :: PDF(:) 
  !real(8), allocatable :: CDF(:) 
  
  real(8), allocatable :: LDAT(:)
end type


type AngularDist
  integer :: flag !> Reaction MT dependent flag !> flag = 0  :: no angular distribution data (isotropic distribution)
                                                !> flag = -1 :: angular distribution data are given in DLW Block (LAW=44)
                                                !> flag > 0  :: angular distribution data are given
  integer :: NE     !> number of energies at which angular distributions are tabulated
  real(8), allocatable :: E(:)         !> energy grid, E(IE), IE=1,NE
  integer, allocatable :: dist_flag(:) !> distribution flag, dist_flag(IE), IE=1,NE
                                       !> dist_flag(IE) = 1  :: 32 equiprobable bin distribution
                                       !> dist_flag(IE) = -1 :: tabular angular distribution  
                                       !> dist_flag(IE) = 0  :: isotropic
  type (AngularDistDataForm), allocatable :: dist(:) !> angular distribution array, dist(IE), IE=1,NE
end type


type, extends (TabularDataForm) :: EnergyDistDataForm
  integer :: law
  integer :: IDAT
  real(8), allocatable :: LDAT(:)
end type


type EnergyDist
  integer :: nlaw
  type (EnergyDistDataForm), allocatable :: dist(:) !> energy distribution array, dist(I), I=1,nlaw
end type


type, extends (TabularDataForm) :: PrecursorDataForm
  real(8) :: decay_const
end type

type UNRtype 
    logical :: URES = .false.    !> Logical indicator that the XS has probability table
    integer :: N        !> Number of incident energies where there is a probability table
    integer :: M        !> Length of table; i.e., number of probabilities, typically 20
    integer :: INT      !> Interpolation parameter between tables =2 lin-lin; =5 log-log
    integer :: ILF, ILFidx       !> Inelastic competition flag
    integer :: IOA, IOAidx       !> Other absorption flag
    integer :: IFF      !> Factors flag
    real(8), allocatable :: E(:), P(:,:,:)
    real(8) :: Emin, Emax
    integer :: unridx
endtype 


real(8) :: EUmin, EUmax !> Min/Max Energy for whole UNR

integer :: n_unr = 0
integer, allocatable :: uresiso(:)
integer :: nueg  = 0
real(8) :: UEGMAX = 3D1

!> 22/08/05 ~ UNIONIZED ENERGY GRID
type UNITED
    integer, allocatable ::  Egrid(:) 
    real(8), allocatable ::   sigt(:) !> Total XS in united grid
    real(8), allocatable ::  sigel(:) !> Elastic XS in ...
    real(8), allocatable ::  sig2n(:) !> (n,2n) XS in ...
    real(8), allocatable ::  sig3n(:) !> (n,3n) XS in ...
    real(8), allocatable ::  sig4n(:) !> (n,4n) XS in ...
    real(8), allocatable ::   sigd(:) !> (n,gamma) XS
    real(8), allocatable ::   sigf(:) !> Fission XS in ...
    real(8), allocatable :: signuf(:) !> Fission Neutron Production
    real(8), allocatable ::  sigqf(:) !> Fission heating ..
    real(8), allocatable ::  sigal(:) !> (n,Alpha) XS
    real(8), allocatable ::   sigp(:) !> (n,p) XS
endtype


!Nuclear data library in ace format
type AceFormat
  character(200) :: library       !> name of library for each isotope
  integer :: ZAID                !> ZAID number
  logical :: excited             !> True if excited
  integer :: NXS(1:16)           !> number array in ace format
  integer :: JXS(1:32)           !> pointer array in ace format
  real(8) :: temp                !> temperature in [MeV]
  real(8) :: atn                 !> ratio of atomic mass to neutron mass
  integer :: sab_iso = 0         !> not zero, which isotope considered S(a,b)
  integer :: resonant = 0        !> resonant isotope?
  character(20) :: xslib         !> xslib indicator (.80c, ...)

  !Data blocks in ace format
  !> ESZ_Block // FIS_Block
  real(8), allocatable :: E(:)     !> energies, E(I), I=1,NXS(3)
  real(8), allocatable :: sigt(:)  !> total cross sections, E(I), I=1,NXS(3)
  real(8), allocatable :: sigd(:)  !> disappearance cross sections, sigd(I), I=1,NXS(3)
  real(8), allocatable :: sigel(:) !> elastic cross sections, sigel(I), I=1,NXS(3)
  real(8), allocatable :: sigf(:)  !> fission cross sections, sigf(I), I=IE, IE+NE-1
  real(8), allocatable :: sigqf(:) !> Q-fission cross sections
  real(8), allocatable :: siga(:)  !> (ADDITIONAL) siga = sigf + sigd
  real(8), allocatable :: H(:)     !> average heating numbers, H(I), I=1,NXS(3)

  real(8) :: qval
  
  !> NU_block
  logical :: nu_block_exist = .true.      !> existence of NU block
  logical :: nu_del_block_exist = .true.      !> existence of NU DEL block
  integer :: nu_tot_flag         !> 1 = polynomial function flag, 2 = tabular data flag
  integer :: nu_del_flag         !> 1 = polynomial function flag, 2 = tabular data flag
  type(NuTotDataForm) :: nu_tot 
  type(TabularDataForm) :: nu_del                       !> exist if JXS(24) > 0 
  type(PrecursorDataForm), allocatable :: prcr(:)       !> exist if JXS(24) > 0

  !> MTR Block
  integer, allocatable :: MT(:) !> ENDF MT numbers, MT(I), I=1,NXS(4)

  !> LQR Block
  real(8), allocatable :: Q(:) !> Q-value of reaction MT, Q(I), I=1,NXS(4)

  !> TYR Block
  integer, allocatable :: TY(:) !> Neutron release for reaction MT, TY(I), I=1,NXS(4)

  !> SIG Block
  type (CrossSectionDataForm), allocatable :: sig_MT(:) !> cross section arrays for reaction MT(I), sig_MT(I), I=1,NXS(4)

  !> AND Block
  integer, allocatable :: ang_flag(:)       !> flags for angular distribution arrays, for reaction MT(I), ang_flag(I), I=0,NXS(5)
  type (AngularDist), allocatable :: ang(:) !> angular distribution arrays, for reaction MT(I), ang_flag(I), I=0,NXS(5)

  !> DLW Block // Delayed Neutron Energy Distribution // DLWP Block
  type (EnergyDist), allocatable :: pneg(:)  !> prompt neutron energy distribution arrays, for reaction MT(I), pneg(I), I=1,NXS(5)
  type (EnergyDist), allocatable :: dneg(:)  !> delayed neutron energy distribution arrays, for delayed neutron group I, dneg(I), I=1,NXS(8)
  type (EnergyDist), allocatable :: ppeg(:)  !> prompt photon energy distribution arrays, for reaction MT(I), ppeg(I), I=1,NXS(8)

  !> Energy-Dependent Neutron Yields
  type(TabularDataForm), allocatable :: nyd(:)  !> neutron yield data, nyd(I), I=1,NXS(4) 

  !> UNR Block
  type(UNRtype) :: UNR

  !> UNIONIZED ENERGY GRID
  type(UNITED) :: UEG
  
  integer :: isab
  integer :: iso0K

  logical :: depletable = .false.

end type
type (AceFormat), allocatable, target :: ace(:)
type (AceFormat), allocatable, target :: ace_base(:)
integer :: num_iso = 0              !> total number of isotopes


! Hash-based Energy Search algorithm 
real(8) :: Emin, Emax = 0D0
integer, allocatable :: ugrid(:,:),&     !Lethargy-grid for hash-based search
                      & ugrid0K(:,:),&   !Lethargy-grid for hash-based search
                      & ugridsab(:,:)
real(8) :: udelta 
integer :: nugrid

! Hash-based E-search for UEG
integer, allocatable :: unigrid(:)
integer :: nuni = 8192
real(8) :: unidel


real(8), allocatable :: ueggrid(:)


! Doppler broadening
! on-the-fly Doppler broadening in resolved resonance region via Gauss Hermite
! method (Y.G. Jo, KNS 2017)
logical :: do_OTFDB
integer :: scat_kernel  ! scattering kernel : cons = 0, wcm = 1, dbrc = 2
real(8) :: DBRC_E_min,  DBRC_E_max  ! (MeV) For U-238, min = 0.4eV, max = 210eV
integer :: n_iso0K = 0      ! # of isotopes treated by exact scattering kernel 
type AceFormat0K
  character(200) :: library      !> name of library for each isotope
  character(200) :: xslib
  integer :: ZAID       ! ZAID number
  integer :: NXS(1:16)  ! number array in ace format
  integer :: JXS(1:32)  ! pointer array in ace format
  real(8) :: temp       ! temperature in [MeV]
  real(8) :: atn        ! ratio of atomic mass to neutron mass

  !Data block to read
  !> ESZ Block
  real(8), allocatable :: ERG(:)    ! energy grid
  real(8), allocatable :: XS0(:)    ! cross section at 0K

end type
type (AceFormat), pointer :: ace0Kptr
type (AceFormat0K), allocatable, target :: ace0K(:)



! S(alpha,beta) : scattering law table
type SAB_INEL_XS
    integer:: NE
    real(8), allocatable:: ERG(:)   !> energy
    real(8), allocatable:: XS(:)    !> cross section
end type

type SAB_EL_XS
    integer:: NE
    real(8), allocatable:: ERG(:)   !> energy
    real(8), allocatable:: XS(:)    !> cross section
end type

type SAB_INEL_E
    real(8), allocatable:: ERG(:,:)     !> (Ein,Eout)
    real(8), allocatable:: ANG(:,:,:)   !> (Ein,Eout,ang)
end type

type SAB_EL_ANG
    real(8), allocatable:: ANG(:,:)     !> (Ein,ang)
end type

type SAB_ACEFORMAT
    character(200):: library
    integer:: NXS(1:16)
    integer:: JXS(1:32)
    integer:: ZAID                !> ZAID number
    real(8):: temp                !> temperature in [MeV]
    real(8):: atn                 !> ratio of atomic mass to neutron mass

    type(SAB_INEL_XS):: ITIE
    type(SAB_EL_XS)::   ITCE
    type(SAB_INEL_E)::  ITXE
    type(SAB_EL_ANG)::  ITCA

    character(20) :: xslib

end type
type(Sab_AceFormat), allocatable, target:: sab(:)
integer:: sab_iso = 0   !> No. of isotopes considering thermal scattering S(a,b)

real(8), allocatable :: XSS(:)        !> temporary XSS array for currently reading isotope
real(8), allocatable :: sab_XSS(:,:)  !> temporary XSS array for currently reading isotope for S(alfa,beta)

type THERM_HEADER ! Identical to therm option in Serpent 2 ( ref. to wiki )
    character(20) :: tag                !> Tag for id (CE and SAB)
    character(20) :: lib_low, lib_high  !> Library for processing
    integer :: issab                    !> Determines whether they're sab
    integer :: iso_low, iso_high        !> ace (or sab) index for low and high E
    real(8) :: temp                     !> Temperature in K * Boltzmann
    real(8) :: f                        !> Fraction used for MAKXSF
end type
type(THERM_HEADER), allocatable, target :: therm(:)
integer :: therm_iso = 0


! Doppler broadening
type DOPPLER_BROADEN
    character(20):: library
    real(8):: Elow
    real(8):: Ehigh
    real(8), allocatable:: ERG(:)   ! energy grid
    real(8), allocatable:: XS0(:)   ! cross section at 0K
end type
type(DOPPLER_BROADEN), allocatable:: DB(:)
integer:: db_iso    ! # of isotopes for DBRC




!Fission ZAIDS for fission Q-values.
integer, parameter :: mfiss(22) =  & 
& (/     90232,     91233,     92233,     92234,     92235, &
&        92236,     92237,     92238,     92239,     92240, &
&        93237,     94238,     94239,     94240,     94241, &
&        94242,     94243,     95241,     95242,     95243, &
&        96242,     96244 /)

!Fission Q-values. [MeV]
real(8), parameter :: qfiss(23) =  & 
& (/ 171.91d+0, 175.57d+0, 180.84d+0, 179.45d+0, 180.88d+0, &
&    179.50d+0, 180.40d+0, 181.31d+0, 180.40d+0, 180.40d+0, &
&    183.67d+0, 186.65d+0, 189.44d+0, 186.36d+0, 188.99d+0, &
&    185.98d+0, 187.48d+0, 190.83d+0, 190.54d+0, 190.25d+0, &
&    190.49d+0, 190.49d+0, 180.00d+0 /)

!MT numbers for reactions in ENDF Library
integer, parameter :: &
& ENDF_TOT = 1,   ENDF_ELASTIC = 2,     ENDF_INELASTIC = 4, &
& ENDF_N2N = 16,  ENDF_N3N = 17, &
& ENDF_FISS = 19, ENDF_DISAPPEAR = 101, ENDF_NG = 102,  ENDF_NP = 103, &
& ENDF_NALFA = 107


!Gauss-Hermite Quadratures for on-the-fly Doppler broadening
real(8),parameter :: ghq(1:16) = (/-4.688738939305818364688, &
 -3.869447904860122698719, -3.176999161979956026814, &
 -2.546202157847481362159, -1.951787990916253977435, & 
 -1.380258539198880796372, -0.8229514491446558925825, &
 -0.2734810461381524521583, 0.2734810461381524521583, &
  0.8229514491446558925825, 1.380258539198880796372, &
  1.951787990916253977435,  2.546202157847481362159, &
  3.176999161979956026814,  3.869447904860122698719, &
  4.688738939305818364688/)
real(8),parameter :: wghq(1:16) = (/2.65480747401118224471d-10, &
 2.32098084486521065339d-7, 2.71186009253788151202d-5, &
 9.32284008624180529914d-4, 0.01288031153550997368346d0, &
 0.0838100413989858294154d0, 0.2806474585285336753695d0, &
 0.5079294790166137419135d0, 0.5079294790166137419135d0, &
 0.2806474585285336753695d0, 0.0838100413989858294154d0, &
 0.01288031153550997368346d0, 9.32284008624180529914d-4, &
 2.71186009253788151202d-5, 2.32098084486521065339d-7, &
 2.65480747401118224471d-10/)
real(8),parameter :: ghq2(1:16) = (/1.94840741569E-01, &
 5.84978765436E-01, 9.76500463590E-01, 1.37037641095E+00, &
 1.76765410946E+00, 2.16949918361E+00, 2.57724953773E+00, &
 2.99249082500E+00, 3.41716749282E+00, 3.85375548547E+00, &
 4.30554795335E+00, 4.77716450350E+00, 5.27555098652E+00, &
 5.81222594952E+00, 6.40949814927E+00, 7.12581390983E+00/)
real(8),parameter :: xghq2(1:16) = (/3.7962914575E-02, &
 3.4220015601E-01, 9.5355315539E-01, 1.8779315077E+00, &
 3.1246010507E+00, 4.7067267077E+00, 6.6422151797E+00, &
 8.9550013377E+00, 1.1677033674E+01, 1.4851431342E+01, &
 1.8537743179E+01, 2.2821300694E+01, 2.7831438211E+01, &
 3.3781970488E+01, 4.1081666525E+01, 5.0777223878E+01/)
real(8),parameter :: wghq2(1:16) = (/1.42451415249E-02, &
 9.49462195824E-02, 1.44243732244E-01, 1.13536229019E-01, &
 5.48474621706E-02, 1.72025699141E-02, 3.56200987793E-03, &
 4.85055175195E-04, 4.26280054876E-05, 2.33786448915E-06, &
 7.59830980027E-08, 1.35405428589E-09, 1.17309796257E-11, &
 4.04486402497E-14, 3.79255121844E-17, 3.71215853650E-21/)

 
 
    contains 
    function find_ACE_iso_idx (this, iso_id) result (idx) 
        type(AceFormat) :: this(:) 
        character(*) :: iso_id 
        integer :: i, idx
        
        do i = 1, size(this) 
            if (trim(this(i)%library) == iso_id) then 
                idx = i 
                return 
            endif
        enddo 
        print *, "ERROR :: no such isotope id : ", iso_id 
        stop 
        
    end function 

    function find_ACE_iso_idx_zaid (zaid, temp) result (idx) 
        use constants, only: K_B
        integer, intent(in) :: zaid 
        real(8), intent(in), optional :: temp
        integer :: i, idx
        integer :: hund
        !do i = 1, num_iso
        !    if ((ace(i)%zaid == zaid/10) .and. ((mod(zaid,10)/=0) .eqv. ace(i)%excited)) then 
        !        idx = i 
        !        return 
        !    endif
        !enddo
        
        if(mod(zaid,10)==0) then !stable
            do i = 1,num_iso
                if(ace(i)%zaid==zaid/10) then
                    idx = i
                    return
                endif
            enddo
        elseif(mod(zaid,10)==1) then !1st meta.
            do i = 1,num_iso
                if(.not. ace(i)%excited) cycle
                hund = 3-mod(zaid/1000,10)
                if(ace(i)%zaid-hund*100==zaid/10) then
                    if(present(temp)) then
                        if(abs(ace(i)%temp-temp)<=1E-3*K_B) then
                            idx = i
                            return
                        endif
                    else
                        idx = i
                        return
                    endif
                endif
            enddo
        endif

        idx = 0 
        !print *, "ERROR :: no such isotope ZAID number : ", zaid 
        !stop 
        
    end function 


end module 
