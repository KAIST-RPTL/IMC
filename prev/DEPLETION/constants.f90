module constants
    implicit none 
    
    !==============================================================================
    !Constants for I/O unit
    integer, parameter :: rd_geom  = 1, &  !I/O unit for reading problem
                          rd_ctrl  = 2, &  !I/O unit for reading control param.
                          rd_xs    = 3, &
                          rd_tally = 4, &
                          rd_mat   = 5, & 
                          rd_dep   = 6, &
                          rd_inven = 7, &
                          rd_temp  = 0, &
                          rd_xslib = 99, &
                          rd_mgt   = 98
    integer, parameter :: wt_coord = 8  !I/O unit for coordinate test 
    integer, parameter :: prt_spec = 9  !I/O unit for spectral analysis variable print 
    
    !integer, parameter :: rd_mg           = 3  !I/O unit for reading control param.
    !integer, parameter :: rd_event        = 4  !I/O unit for event for tmc
    !integer, parameter :: rd_burn         = 5  !I/O unit for burnup calculation
    integer, parameter :: prt_keff        = 10  !I/O unit for keff
    integer, parameter :: prt_ntpy        = 11  !I/O unit for ntpy 
    integer, parameter :: prt_flux        = 12  !I/O unit for flux
    integer, parameter :: prt_powr        = 13  !I/O unit for power
    integer, parameter :: prt_bumat       = 14  !I/O unit for burned materials
    integer, parameter :: prt_dynamic     = 15  !I/O unit for burned materials
    integer, parameter :: prt_prompt      = 16  !I/O unit for burned materials
    integer, parameter :: prt_delayed     = 17  !I/O unit for burned materials
    integer, parameter :: prt_prec        = 18  !I/O unit for burned materials
    integer, parameter :: prt_wgt         = 19  !I/O unit for burned materials
    integer, parameter :: prt_tet_vrc     = 20  !I/O unit for burned materials
    
    integer, parameter :: bumat_test      = 21    
    integer, parameter :: prt_adjoint     = 22
    integer, parameter :: prt_fuel_mv     = 23
    !integer, parameter :: prt_assflux      = 10 !I/O unit for flux
    !integer, parameter :: prt_asspowr      = 11 !I/O unit for power
    !integer, parameter :: prt_3Dflux      = 12 !I/O unit for flux
    !integer, parameter :: prt_3Dpowr      = 13 !I/O unit for power
    !integer, parameter :: prt_HB_2Dflux  = 14 !I/O unit for flux
    !integer, parameter :: prt_HB_2Dpowr  = 15 !I/O unit for power
    !integer, parameter :: prt_HB_assflux = 16 !I/O unit for flux
    !integer, parameter :: prt_HB_asspowr = 17 !I/O unit for power
    !integer, parameter :: prt_HB_3Dflux  = 18 !I/O unit for flux
    !integer, parameter :: prt_HB_3Dpowr  = 19 !I/O unit for power
    !integer, parameter :: monitor         = 20 !I/O unit for monitor to text
    !integer, parameter :: debug1          = 21 !I/O unit for debugging
    !integer, parameter :: debug2          = 22 !I/O unit for debugging
    !integer, parameter :: prt_pwr         = 23 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_pk          = 24 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_pk2         = 25 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_pk3         = 26 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_pop         = 27 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_tsize       = 28 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_mg          = 29 !I/O unit for transient Monte Carlo
    !integer, parameter :: prt_time        = 30 !I/O unit for time analysis
    !integer, parameter :: prt_rSF         = 31 !I/O unit for convergence analysis
    !integer, parameter :: prt_pCMFD       = 32 !I/O unit for keff from pCMFD calculation
    !integer, parameter :: prt_ratio       = 33 !I/O unit for ratio estimator for keff
    !integer, parameter :: prt_burn1       = 34 !I/O unit for burnup calculation
    !integer, parameter :: prt_burn2       = 35 !I/O unit for burnup calculation
    !integer, parameter :: prt_burn3       = 36 !I/O unit for burnup calculation
    !integer, parameter :: prt_burn4       = 37 !I/O unit for burnup calculation
    !integer, parameter :: prt_power       = 38 !I/O unit for fuel power distribution 
    !integer, parameter :: prt_burnup      = 39 !I/O unit for burnup distribution
    !integer, parameter :: prt_sca1        = 40 !I/O unit for fuel temperature distribution 
    !integer, parameter :: prt_sca2        = 41 !I/O unit for bulk temperature
    !integer, parameter :: prt_sca3        = 42 !I/O unit for bulk density
    !integer, parameter :: prt_cbc         = 43 !I/O unit for critical boron concentration 

    !==============================================================================
    ! Physical constants
    real(8), parameter:: m_u = 1.660540D-27 ! amu to Kg
    real(8), parameter:: m_n = 1.008664     ! neutron mass (amu)
    real(8), parameter:: mevj = 1.6022D-13  ! MeV to Joule
	
	
	
	
    !==============================================================================
    ! Precision kind
	integer, parameter :: sp = kind(1.0), &
						  dp = kind(1.d0)

	
	
    ! universe types 
    integer, parameter :: pure_univ = 0 , &
                          Pin_univ    = 1
    
    ! surface types 
    integer ::  px   = 1, &
                py   = 2, &
                pz   = 3, &
                sqcx = 4, &
                sqcy = 5, &
                sqcz = 6, &
                cylx = 7, &
                cyly = 8, &
                cylz = 9, &
                sph  = 10, & 
				!hexz = 11, &
                hexxc = 11, &
                hexyc = 12, &

                ! DERIVED FORMAT STARTS FROM 3-digits
                gcross = 101

    !==============================================================================
    
  ! Maximum number of collisions/crossings
  integer, parameter :: MAX_EVENTS = 1000000
  integer, parameter :: MAX_SAMPLE = 100000
  
  
  ! Maximum number of secondary particles created
  integer, parameter :: MAX_SECONDARY = 1000

  ! Maximum number of coordinates
  integer, parameter :: MAX_COORD = 10
  
  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_WORDS       = 500
  integer, parameter :: MAX_LINE_LEN    = 250
  integer, parameter :: MAX_WORD_LEN    = 150
  integer, parameter :: MAX_FILE_LEN    = 255

  ! Maximum number of external source spatial resamples to encounter before an
  ! error is thrown.
  integer, parameter :: EXTSRC_REJECT_THRESHOLD = 10000
  real(8), parameter :: EXTSRC_REJECT_FRACTION = 0.05

  
  real(8), parameter :: fes(33) = &  != Fission energy spectrum
    & (/  0.0d0,       .162524d0,  .266043d0,  .358425d0,  .445672d0, &
    &      .530293d0,  .613680d0,  .696783d0,  .780264d0,  .864702d0, &
    &      .950540d0, 1.038286d0, 1.128326d0, 1.221111d0, 1.317206d0, &
    &     1.417070d0, 1.521302d0, 1.630646d0, 1.745929d0, 1.868073d0, &
    &     1.998282d0, 2.138046d0, 2.289259d0, 2.454356d0, 2.636707d0, &
    &     2.840830d0, 3.073518d0, 3.344965d0, 3.672134d0, 4.086420d0, &
    &     4.656234d0, 5.588725d0, 9.0d0  /)

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  ! Values here are from the Committee on Data for Science and Technology
  ! (CODATA) 2014 recommendation (doi:10.1103/RevModPhys.88.035009).

  real(8), parameter ::                      &
       PI               = 3.1415926535898_8, & ! pi
       SQRT_PI          = 1.7724538509055_8, & ! square root of pi
       MASS_NEUTRON     = 1.00866491588_8,   & ! mass of a neutron in amu
       MASS_NEUTRON_EV  = 939.5654133e6_8,   & ! mass of a neutron in eV/c^2
       MASS_PROTON      = 1.007276466879_8,  & ! mass of a proton in amu
       AMU              = 1.660539040e-27_8, & ! 1 amu in kg
       C_LIGHT          = 2.99792458e8_8,    & ! speed of light in m/s
       N_AVOGADRO       = 0.60220434469282e24_8, & !  0.6022140857e24_8,    & ! Avogadro's number in 1/mol
       K_B              = 8.6173303e-11_8,    & ! Boltzmann constant in MeV/K
       INFINITY         = huge(0.0_8),       & ! positive infinity
       TOOLONG          = 1.0e6_8,           & ! too long distance limit
       ZERO             = 0.0_8,             &
       HALF             = 0.5_8,             &
       ONE              = 1.0_8,             &
       TWO              = 2.0_8,             &
       THREE            = 3.0_8,             &
       FOUR             = 4.0_8
  complex(8), parameter :: ONEI = (ZERO, ONE)

  ! ============================================================================
  ! GEOMETRY-RELATED CONSTANTS
  real(8), parameter :: TINY_BIT  = 1e-9_8
  real(8), parameter :: wgt_min   = 1e-2_8 !0.25!1e-2_8
  real(8), parameter :: wgt_split = 2.0d0 
  
  ! FLOATING POINT ARITHMETICS MARGIN 
  real(8), parameter :: EPSILON = 1.0e-10_8  
  

  ! Cell fill types
  integer, parameter ::  &
       FILL_MATERIAL = 1, & ! Cell with a specified material
       FILL_UNIVERSE = 2, & ! Cell filled by a separate universe
       FILL_LATTICE  = 3    ! Cell filled with a lattice

  ! Void material
  integer, parameter :: MATERIAL_VOID = -1

  ! Lattice types
  integer, parameter ::  &
       LATTICE_RECT = 1, & ! Rectangular lattice
       LATTICE_HEX  = 2    ! Hexagonal lattice

  ! Base univer
  integer, parameter :: BASE = 0

  ! Flag to say that the outside of a lattice is not defined
  integer, parameter :: NO_OUTER_UNIVERSE = -1

  ! Maximum number of lost particles
  integer, parameter :: MAX_LOST_PARTICLES = 10

  ! Maximum number of lost particles, relative to the total number of particles
  real(8), parameter :: REL_MAX_LOST_PARTICLES = 1e-6_8

  ! ============================================================================
  ! CROSS SECTION RELATED CONSTANTS
  real(8) :: barn = 1.0d-24
  
  ! Interpolation flag
  integer, parameter ::   &
       HISTOGRAM     = 1, & ! y is constant in x
       LINEAR_LINEAR = 2, & ! y is linear in x
       LINEAR_LOG    = 3, & ! y is linear in ln(x)
       LOG_LINEAR    = 4, & ! ln(y) is linear in x
       LOG_LOG       = 5    ! ln(y) is linear in ln(x)

  ! Particle type
  integer, parameter :: &
       NEUTRON  = 1, &
       PHOTON   = 2, &
       ELECTRON = 3

  ! Angular distribution type
  integer, parameter :: &
       ANGLE_ISOTROPIC = 1, & ! Isotropic angular distribution (CE)
       ANGLE_32_EQUI   = 2, & ! 32 equiprobable bins (CE)
       ANGLE_TABULAR   = 3, & ! Tabular angular distribution (CE or MG)
       ANGLE_LEGENDRE  = 4, & ! Legendre angular distribution (MG)
       ANGLE_HISTOGRAM = 5    ! Histogram angular distribution (MG)

  ! Number of mu bins to use when converting Legendres to tabular type
  integer, parameter :: DEFAULT_NMU = 33

  ! Secondary energy mode for S(a,b) inelastic scattering
  integer, parameter :: &
       SAB_SECONDARY_EQUAL  = 0, & ! Equally-likely outgoing energy bins
       SAB_SECONDARY_SKEWED = 1, & ! Skewed outgoing energy bins
       SAB_SECONDARY_CONT   = 2    ! Continuous, linear-linear interpolation

  ! Elastic mode for S(a,b) elastic scattering
  integer, parameter :: &
       SAB_ELASTIC_DISCRETE = 3, & ! Sample from discrete cosines
       SAB_ELASTIC_EXACT    = 4    ! Exact treatment for coherent elastic

  ! Reaction types
  integer, parameter :: &
       TOTAL_XS = 1,  ELASTIC = 2,   N_NONELASTIC = 3, N_LEVEL = 4, MISC  = 5,   &
       N_2ND   = 11,  N_2N    = 16,  N_3N   = 17,  N_FISSION = 18, N_F    = 19,  &
       N_NF    = 20,  N_2NF   = 21,  N_NA   = 22,  N_N3A   = 23,  N_2NA   = 24,  &
       N_3NA   = 25,  N_NP    = 28,  N_N2A  = 29,  N_2N2A  = 30,  N_ND    = 32,  &
       N_NT    = 33,  N_N3HE  = 34,  N_ND2A = 35,  N_NT2A  = 36,  N_4N    = 37,  &
       N_3NF   = 38,  N_2NP   = 41,  N_3NP  = 42,  N_N2P   = 44,  N_NPA   = 45,  &
       N_N1    = 51,  N_N40   = 90,  N_NC   = 91,  N_DISAPPEAR = 101, N_GAMMA = 102, &
       N_P     = 103, N_D     = 104, N_T    = 105, N_3HE   = 106, N_A     = 107, &
       N_2A    = 108, N_3A    = 109, N_2P   = 111, N_PA    = 112, N_T2A   = 113, &
       N_D2A   = 114, N_PD    = 115, N_PT   = 116, N_DA    = 117, N_5N    = 152, &
       N_6N    = 153, N_2NT   = 154, N_TA   = 155, N_4NP   = 156, N_3ND   = 157, &
       N_NDA   = 158, N_2NPA  = 159, N_7N   = 160, N_8N    = 161, N_5NP   = 162, &
       N_6NP   = 163, N_7NP   = 164, N_4NA  = 165, N_5NA   = 166, N_6NA   = 167, &
       N_7NA   = 168, N_4ND   = 169, N_5ND  = 170, N_6ND   = 171, N_3NT   = 172, &
       N_4NT   = 173, N_5NT   = 174, N_6NT  = 175, N_2N3HE = 176, N_3N3HE = 177, &
       N_4N3HE = 178, N_3N2P  = 179, N_3N3A = 180, N_3NPA  = 181, N_DT    = 182, &
       N_NPD   = 183, N_NPT   = 184, N_NDT  = 185, N_NP3HE = 186, N_ND3HE = 187, &
       N_NT3HE = 188, N_NTA   = 189, N_2N2P = 190, N_P3HE  = 191, N_D3HE  = 192, &
       N_3HEA  = 193, N_4N2P  = 194, N_4N2A = 195, N_4NPA  = 196, N_3P    = 197, &
       N_N3P   = 198, N_3N2PA = 199, N_5N2P = 200, N_P0    = 600, N_PC    = 649, &
       N_D0    = 650, N_DC    = 699, N_T0   = 700, N_TC    = 749, N_3HE0  = 750, &
       N_3HEC  = 799, N_A0    = 800, N_AC   = 849, N_2N0   = 875, N_2NC   = 891

  ! Depletion reactions
  integer, parameter :: DEPLETION_RX(6) = [N_GAMMA, N_P, N_A, N_2N, N_3N, N_4N]
  real(8), parameter :: eVtoJoule = 1.602175314d-19 !Joule/eV

  ! ACE table types
  integer, parameter :: &
       ACE_NEUTRON   = 1, & ! continuous-energy neutron
       ACE_THERMAL   = 2, & ! thermal S(a,b) scattering data
       ACE_DOSIMETRY = 3    ! dosimetry cross sections

  ! MGXS Table Types
  integer, parameter :: &
       MGXS_ISOTROPIC   = 1, & ! Isotropically Weighted Data
       MGXS_ANGLE       = 2    ! Data by Angular Bins

  ! Fission neutron emission (nu) type
  integer, parameter ::   &
       NU_NONE       = 0, & ! No nu values (non-fissionable)
       NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Secondary particle emission type
  integer, parameter :: &
       EMISSION_PROMPT = 1,  & ! Prompt emission of secondary particle
       EMISSION_DELAYED = 2, & ! Delayed emission of secondary particle
       EMISSION_TOTAL = 3      ! Yield represents total emission (prompt + delayed)

  ! Cross section filetypes
  integer, parameter :: &
       ASCII  = 1, & ! ASCII cross section file
       BINARY = 2    ! Binary cross section file




  ! ============================================================================
  ! MISCELLANEOUS CONSTANTS

  ! indicates that an array index hasn't been set
  integer, parameter :: NONE = 0
  integer, parameter :: C_NONE = -1

  ! Codes for read errors -- better hope these numbers are never used in an
  ! input file!
  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8



  !=============================================================================
  ! DELAYED NEUTRON PRECURSOR CONSTANTS

  ! Since cross section libraries come with different numbers of delayed groups
  ! (e.g. ENDF/B-VII.1 has 6 and JEFF 3.1.1 has 8 delayed groups) and we don't
  ! yet know what cross section library is being used when the tallies.xml file
  ! is read in, we want to have an upper bound on the size of the array we
  ! use to store the bins for delayed group tallies.
  integer, parameter :: MAX_DELAYED_GROUPS = 8
    


end module constants
