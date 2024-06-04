
include 'mkl_pardiso.f90'

module depletion_module 
    use variables
    use constants
    use ace_header, only: ace, find_ACE_iso_idx_zaid, num_iso
    use ace_xs,     only: getxs, getierg
    use material_header
    use mpi 
    
    implicit none 
    
    !==============================================================================
    !Nuclide data structure for depletion
    type nuclide_data
      logical :: data_exist !data is read from library
      logical :: reduced = .false. !check whether reduced
      logical :: fiss = .false.
      integer :: idx   !nuclide index
      !neutron interaction
      real(8) :: sng   !(n,g) leading to ground state
      real(8) :: sn2n  !(n,2n) leading to ground state
      real(8) :: sna   !(n,alpha) leading to ground state for fission products and activation products
      real(8) :: snp   !(n,proton) leading to ground state for fission products and activation products
      real(8) :: sn3n  !(n,3n) leading to ground state for actinides
      real(8) :: snf   !(n,f) fission for actinides
      real(8) :: sngx  !(n,g) leading to excited state
      real(8) :: sn2nx !(n,2n) leading to excited state
      integer :: fy_idx!0 : no fission yield data, otherwise: index for fission yield array

      !radioactive decay
      integer :: iu    !Time designation [1=seconds, 2=minutes, 3=hours, 4=days, 5=years, 6=stable, 7=10^3 years, 8=10^6 years, 9=10^9 years]
      real(8) :: lambda!total decay constant [s^-1]
      real(8), allocatable :: frac(:) !Fraction array of decay mode
      integer, allocatable :: daughter(:,:) !Daughter nuclide corresponding to decay mode 
      real(8) :: n_emit
      real(8) :: p_emit
      real(8) :: a_emit !neutron, proton, alpha emission fraction per decay
      integer :: fp=0
      integer :: sfp = 0
      integer :: react_num=0
      real(8) :: qrec
      real(8) :: removal !Removal rate per circulation

      real(8) :: amu !Isotope mass in amu
      !Yield calculation
      real(8),allocatable :: yield_E(:)
      !LOGICAL: appears during depletion
      logical :: isin = .false.
      !211116 : SFY
      integer :: sfy_idx
      integer :: sfiss = 0
    end type
    type(nuclide_data) :: nuclide(0:2,0:170,0:111) !nuclide data for Isomeric State (0=ground, 1=excited), neutron number, atomic number 
    ! NFY TESTING
    real(8), allocatable :: fratio(:,:)
    real(8), allocatable :: yieldE(:,:)
    real(8), allocatable :: yieldnE(:)

    ! Used for chain reduction
    real(8) :: HL_th = 1.d0
    real(8) :: lam_th
    !Nuclide data index
    integer :: zai_idx(1:4000) !ZAI index = Z*10000 + A*10 + I
    integer :: maxnnz          !Maximum allowed number of nonzero elements in burnup matrix
    integer :: nnuc            !Total number of nuclides
    integer :: inuc            !Current nuclide index


    !Fission yield data
    integer :: nfssn           !number of fissionable nuclide
    integer :: nfp             !number of fission product
    integer, allocatable :: fssn_zai(:)     !fissionable zai index
    integer, allocatable :: ace_fssn(:)     !fid with ace_iso
    integer, allocatable :: fp_zai(:)       !fission product zai index
    real(8), allocatable :: yield_data(:,:) !fission yield data (1:nfp, 1:nfssn)
    real(8), allocatable :: tmp_yield(:,:,:) !Energywise yield storage
    real(8), allocatable :: ify_yield(:,:,:) !Indep. Fission yield
    real(8), allocatable :: cfy_yield(:,:,:) !Cumul. Fission yield

    ! SFY DATA
    integer :: nsfssn, nsfp
    integer, allocatable :: sfssn_zai(:), ace_sfssn(:), sfp_zai(:)
    real(8), allocatable :: sfy_yield(:,:) !Spontaneous FY data (211116 Add)

    real(8) :: tot_mass
    real(8) :: tot_fmass
    real(8) :: tot_mass_init, tot_fmass_init

    !Burnup matrix
    real(8), allocatable :: bMat(:,:)  !2-D burnup matrix
    real(8), allocatable :: bMat0(:,:) !2-D burnup matrix (material independent)
    !real(8), allocatable :: bMat_tmp(:,:) 


    !Real Power to normalize flux
    real(8) :: RealPower  ![MW]
    real(8) :: ULnorm     !Unit less normalization factor
    real(8) :: Wnorm      !Power normalization factor [MeV/s to W]


    !Matrix Exponential Solver Option
    integer :: Matrix_Exponential_Solver !0=Chebyshev Rational Approximation Method (CRAM)


    !Chebyshev Rational Approximation
    integer :: cram_order
    logical :: cram_init
    integer :: job(1:8) !Array contains the conversion parameters
    complex(8), allocatable :: Acsr(:) !Non-zero elements of matrix A
    integer, allocatable :: iAcsr(:), jAcsr(:) !Non-zero indices of matrix A


    !Depletion time interval and step
    integer :: &
    & nstep_burnup, &             !Number of time step for burnup calculation
    & istep_burnup = 0            !Step index for burnup calculation
    real(8), allocatable :: &
    & burn_step(:)               !Burnup time for each step [sec]
    real(8) :: bstep_size          !Burnup time interval for each step [sec]
    logical :: auto_burnup      !Do burnup calculation until k<1.00300 or cbc < 30 ppm (1ppm -> 10 pcm)
    real(8) :: auto_bstep_size  !Burnup time interval for remaining burnup calculation [sec]
    real(8) :: total_HMmass     !Total initial heavy metal mass [kgHM]
    real(8) :: total_burnup     !Total burnup [MWday/kgHM]    integer :: Matrix_Exponential_Solver
    
    !Location of depletion library 
    character(100) :: dep_lib
    
    !Removal/Refueling operator
    real(8) :: eff_gas = 0.d0
    real(8) :: eff_noble = 0.d0
    logical :: refuel = .false.
    integer :: n_rf
    real(8),allocatable :: r_rf(:)
    integer,allocatable :: anum_rf(:), mnum_rf(:)
    !==============================================================================
    !In-line xenon equilibrium search for depletion calculation
    logical :: do_Xe_search =.false.                  !Inline Xenon and Iodine search
    integer :: Xe135_iso                      !Index for Xe135 in global list
    integer :: Xe_skp                         !Number of skipped iteration for Xenon tally
    integer :: Xe_sfifo                       !Length of fifo queue length for Xe_tally

    integer :: I135_iso

    integer, allocatable :: Xe_pointer(:)     !Pointer of Xe_tally = loc
    integer, allocatable :: Xe_source(:)      !Potential source of Xenon
    real(8), allocatable :: Xe_prod(:)        !Xe_prod(loc) = cumulative xenon production
    real(8), allocatable :: Xe_trans(:)       !Xe_trans(loc+1) = cumulative xenon transmutation
     
    real(8), allocatable :: I_prod(:)
    
    real(8), allocatable :: Xe_numden(:)      !Xenon number densities [#/cc*1.d-24]
    real(8), allocatable :: avg_Xe_numden(:)  !Averaged number densities of Xenon [#/cc*1.d-24]
    integer, allocatable :: Xe135_mt_iso(:)   !Index for Xenon135 in local list
    real(8), allocatable :: Xe_cum_yield(:)   !Cumulative xenon yield : summation over independent yields of Sn-135, Sb-135, Te-135, I-135, Xe-135, Xe-135m
    !Modification: In-135,136, Sn-135,137, Sb-135~137, Te-135,136, I-135,Xe-135(/m)
    !Doctoral thesis of HS Lee (UNIST MCS)

    real(8), allocatable :: I_numden(:)
    integer, allocatable :: I135_mt_iso(:)
    real(8), allocatable :: I_cum_yield(:)    !Cumulative Iodine-135 yield
    
    logical :: depmtx = .true.
    
    ! ISOMERIC BRANCHING RATIO (UNDER TESTING 10/18)
    integer, allocatable :: ZAIMT_ism(:) ! 1: ZAI, ! IGNORED MT: all 102 2:mT
    real(8), allocatable :: gnd_frac(:)   ! ground fraction
    integer :: num_brn ! Number of isom. branches
    integer,parameter :: ngrid = 1E6
    contains 

    subroutine getENDFdepletionlibrary
        implicit none
        logical :: file_exists
        character(80)  :: line0, line1, line2
        integer, parameter :: rd_decay = 20160818
        integer, parameter :: rd_1gcx  = 20160819
        integer, parameter :: rd_yield = 20160820
        integer, parameter :: rd_sfy   = 20211116
        integer, parameter :: isom_brn = 20211018
        integer :: nuclid
        integer :: anum, mnum, inum, nnum !atomic number, mass number, isomer state, neutron number of nuclide
        integer :: nlb
        real(8) :: thalf, yyn
        integer :: i, j, k
        
        integer :: MT, MAT
        integer :: nskip
        real(8) :: skip
        integer :: n_rad, rad
        integer :: n_react, react
        integer :: rt,typedet,ri,fid
        real(8) :: ZA,HL,ST,r_type,r_q,r_ratio,r_iso
        real(8) :: Etmp1,Etmp2,Etmp3
        integer,allocatable :: prod(:)
        real(8) :: Eavg !average E. only requires if fpy_option =2 
        
        integer :: eg,n_y_eg
        integer :: fp,ifp,flag
        integer :: n_fp_max, n_fp
        real(8) :: i_real, real_rad
        real(8), allocatable :: Ep(:)
        real(8), allocatable :: coeff(:)

        real(8), allocatable :: Etmp(:), tmp(:)
        real(8) :: yield
        integer :: g, n_E
        !21/11/16 SFY
        integer :: n_sfp_max, sfp, isfp, sfid
        real(8) :: sfy
        !21/10/18 CFY and ISOMERIC branching TESTING
        real(8) :: cfy, ify
        integer :: brn, tmpfp

        integer :: fssn_nnum, fssn_anum, fssn_inum
        integer :: fp_tmp, zai
        integer :: isofp, gndfp, prevnuclid
        integer, parameter :: yield_option = 1

        ! Removal
        integer,allocatable :: rem_gas(:), rem_noble(:)

        ! Reduction
        integer :: n_dau, tmp_n, n_grand
        integer,allocatable :: prod_red(:)
        real(8) :: lam
        real(8), allocatable :: frac_dau(:), frac_grand(:)
        integer, allocatable :: prod_dau(:,:), prod_grand(:,:)
        real(8) :: y_sum, f_tmp
        integer :: fp_dau, fp_prod
        real(8) :: frac_tmp
        integer, allocatable :: prod_tmp(:)

        integer :: mt_iso,iso
        !==============================================================================
        !==============================================================================
        !Initialization
        do inum = 0, 2
        do nnum = 0, 170
        do anum = 0, 111
          nuclide(inum,nnum,anum)%data_exist = .false.
          nuclide(inum,nnum,anum)%fy_idx = 0 
          nuclide(inum,nnum,anum)%idx = 0 
          nuclide(inum,nnum,anum)%removal = 0.d0
        end do
        end do
        end do
        allocate(rem_gas(2)); allocate(rem_noble(10))
        rem_gas = (/53,36/)
        rem_noble = (/34,41,42,43,44,45,46,47,51,52/)
        !==============================================================================
        !Read decay library
        open(rd_decay,file=trim(dep_lib)//'sss_endfb71.dec',action='read')
        10  continue
        ! Skipping beginning
        MT = 0
        do while (MT /= 457)
            read(rd_decay,'(A80)',end = 100) line0
            read(line0(73:75),'(i)') MT
        enddo
        ! Read 1st line: get ZAI
        read(line0(1:11),'(f)') ZA
        read(line0(34:44),'(i)') inum
        if(inum>2) goto 15
        read(line0(56:66),'(f)') real_rad
        n_rad = int(real_rad)
        nuclid = int(ZA)
        anum = nuclid/1000
        mnum = nuclid-anum*1000
        nnum = mnum-anum
        if (nuclide(inum,nnum,anum)%data_exist) go to 15
        nuclide(inum,nnum,anum)%data_exist = .true.
         
        nuclide(inum,nnum,anum)%sng = 0.d0
        nuclide(inum,nnum,anum)%sn2n = 0.d0
        nuclide(inum,nnum,anum)%sna = 0.d0
        nuclide(inum,nnum,anum)%snp = 0.d0
        nuclide(inum,nnum,anum)%sn3n = 0.d0
        nuclide(inum,nnum,anum)%snf = 0.d0
        nuclide(inum,nnum,anum)%sngx = 0.d0
        nuclide(inum,nnum,anum)%sn2nx = 0.d0
        nuclide(inum,nnum,anum)%fy_idx = 0.d0
        read(line0(12:22),'(f)') nuclide(inum,nnum,anum)%amu
        if (ANY(rem_gas==anum)) nuclide(inum,nnum,anum)%removal = eff_gas
        if (ANY(rem_noble==anum)) nuclide(inum,nnum,anum)%removal = eff_noble
        ! Initialize fractions; some requires summation
        nuclide(inum,nnum,anum)%qrec = 0.0
        ! Read 2nd line: get HL, determine skipping next line
        read(rd_decay,'(A80)',end = 100) line0
        read(line0(1:11),'(f)') HL !> HL [sec]
        read(line0(45:55),'(f)') ST !> 1 if stable, 0 if radioactive
        if (HL == 0.d0) then
            nuclide(inum,nnum,anum)%lambda = 0.0
            nuclide(inum,nnum,anum)%iu = 6
        else
            nuclide(inum,nnum,anum)%lambda = 6.931471806d-1/HL
            nuclide(inum,nnum,anum)%iu = 1
        endif
        read(line0(45:55),'(f)') skip
        nskip = ceiling(skip/6.0) ! Number of lines to skip
        if (nskip > 0) then
            do i = 1,nskip
                read(rd_decay,'(A80)',end=100) line0
                read(line0(1:11),'(f)') Etmp1
                read(line0(23:33),'(f)') Etmp2
                read(line0(45:55),'(f)') Etmp3
                nuclide(inum,nnum,anum)%qrec = &
                    nuclide(inum,nnum,anum)%qrec +&
                    (Etmp1+Etmp2+Etmp3)/1000000
            enddo
        endif
        read(rd_decay,'(A80)',end=100) line0
        read(line0(56:66),'(i)') n_react ! Number of reactions
        allocate(nuclide(inum,nnum,anum)%frac(n_react))
        allocate(nuclide(inum,nnum,anum)%daughter(n_react,3))
        nuclide(inum,nnum,anum)%daughter=0
        nuclide(inum,nnum,anum)%n_emit = 0.d0
        nuclide(inum,nnum,anum)%p_emit = 0.d0
        nuclide(inum,nnum,anum)%a_emit = 0.d0
        if (n_react==0) go to 15
        nuclide(inum,nnum,anum)%react_num = n_react
        do react = 1,n_react
            if(.not.allocated(prod)) allocate(prod(3))
            prod = (/inum,nnum,anum/)
            flag = 4
            read(rd_decay,'(A80)',end=100) line1
            read(line1(1:11),'(f)') r_type
            read(line1(12:22),'(f)') r_iso
            read(line1(23:33),'(f)') r_q
            read(line1(45:55),'(f)') r_ratio
            nuclide(inum,nnum,anum)%frac(react) = r_ratio
            rt = floor(r_type)
            do while(rt<r_type)
                r_type = r_type * 10.d0
                rt = floor(r_type)
            enddo
            ri = int(r_iso) 
            14 continue
            typedet = unitdigit(rt)
            select case(typedet)
            case(1) !beta- decay: n to p
                prod(1) = ri
                prod(2) = prod(2) - 1
                prod(3) = prod(3) + 1
            case(2) !electron capture: p to n
                prod(1) = ri
                prod(2) = prod(2) + 1
                prod(3) = prod(3) - 1
            case(3) !IT
                prod(1) = ri
            case(4) !alpha emission: emit a, a-2, n-2
                nuclide(inum,nnum,anum)%a_emit = &
                    nuclide(inum,nnum,anum)%a_emit + r_ratio
                prod(1) = ri
                prod(2) = prod(2) - 2
                prod(3) = prod(3) - 2
            case(5) !neutron emission: emit n, n-1
                nuclide(inum,nnum,anum)%n_emit = &
                    nuclide(inum,nnum,anum)%n_emit + r_ratio
                prod(1) = ri
                prod(2) = prod(2) - 1
            case(6)
                ! Spontaneous fission: prod = 0
                nuclide(inum,nnum,anum)%sfiss = react
                prod    = 0
            case(7) !Proton emission
                nuclide(inum,nnum,anum)%p_emit = &
                    nuclide(inum,nnum,anum)%p_emit + r_ratio
                prod(1) = ri
                prod(3) = prod(3) - 1
            case default
            end select
            if (rt>10) then !If further decay exists
                rt = rt/10
                go to 14
            endif
            if (prod(1)>=0 .and. prod(2)>=0 .and. prod(3)>=0) then
                nuclide(inum,nnum,anum)%daughter(react,1) = prod(1)
                nuclide(inum,nnum,anum)%daughter(react,2) = prod(2)
                nuclide(inum,nnum,anum)%daughter(react,3) = prod(3)
            endif
        enddo
        15  continue
        do while (MT /= 451)
            read(rd_decay,'(A80)',end = 100) line0 ! Skip until next isotope's description (MT=451)
            read(line0(73:75),'(i)') MT 
        enddo
        read(rd_decay,*,end=100)
        go to 10
        100 close(rd_decay)
        
        !==============================================================================
        !Read one-group transmutation cross section library
        open(rd_1gcx, file=trim(dep_lib)//"1gcx_library",action="read")
        20 continue
        read(rd_1gcx,'(A80)',end=200) line0 !Title
        do
          read(rd_1gcx,'(A80)',end=200) line1
          
          read(line1(1:4),'(i)') nlb; if(nlb==-1) go to 20 !library index (1=activation products, 2=actinides, 3=fission products)
          read(line1(7:12),'(i)') nuclid
            anum = nuclid/10000
            mnum = (nuclid - anum*10000)/10
            nnum = mnum - anum
            inum = nuclid - anum*10000 - mnum*10
          read(line1(14:22),'(f)') nuclide(inum,nnum,anum)%sng  !(n,g) leading to ground state
          nuclide(inum,nnum,anum)%sng=nuclide(inum,nnum,anum)%sng*1.d-24
          read(line1(24:32),'(f)') nuclide(inum,nnum,anum)%sn2n !(n,2n) leading to ground state
          nuclide(inum,nnum,anum)%sn2n=nuclide(inum,nnum,anum)%sn2n*1.d-24
          if(nlb/=3) then !fission products or activation products
            read(line1(34:42),'(f)') nuclide(inum,nnum,anum)%sna !(n,alpha) leading to ground state
            nuclide(inum,nnum,anum)%sna=nuclide(inum,nnum,anum)%sna*1.d-24
            read(line1(44:52),'(f)') nuclide(inum,nnum,anum)%snp !(n,proton) leading to ground state
            nuclide(inum,nnum,anum)%snp=nuclide(inum,nnum,anum)%snp*1.d-24
          else
            read(line1(34:42),'(f)') nuclide(inum,nnum,anum)%sn3n!(n,3n) leading to ground state
            nuclide(inum,nnum,anum)%sn3n=nuclide(inum,nnum,anum)%sn3n*1.d-24
            read(line1(44:52),'(f)') nuclide(inum,nnum,anum)%snf !(n,f)
            nuclide(inum,nnum,anum)%snf=nuclide(inum,nnum,anum)%snf*1.d-24
          end if 
          read(line1(54:62),'(f)') nuclide(inum,nnum,anum)%sngx  !(n,g) leading to excited state
          nuclide(inum,nnum,anum)%sngx=nuclide(inum,nnum,anum)%sngx*1.d-24
          read(line1(64:72),'(f)') nuclide(inum,nnum,anum)%sn2nx !(n,2n) leading to excited state
          nuclide(inum,nnum,anum)%sn2nx=nuclide(inum,nnum,anum)%sn2nx*1.d-24
          read(line1(74:79),'(f)') yyn               !yyn > 0 : fission yield card follows, yyn < 0 : no fission yield card
          if(yyn>0.d0) read(rd_1gcx,'(A80)',end=200) line2
        end do
        200 close(rd_1gcx)


        !==============================================================================
        open(rd_yield, file = trim(dep_lib)//'sss_endfb71.nfy', action = 'read')
        fid = 0; ifp = 0
        allocate(tmp_yield(1:2500,1:4,1:100)); tmp_yield = 0.d0
        allocate(ify_yield(1:2500,1:4,1:100)); ify_yield = 0.d0
        allocate(cfy_yield(1:2500,1:4,1:100)); cfy_yield = 0.d0
        allocate(fp_zai(1:2500)) ! ZAI of fp
        allocate(fssn_zai(1:100)) ! ZAI of fssn
        allocate(ace_fssn(1:num_iso)); ace_fssn = 0.d0
        allocate(yieldE(1:100,4)); yieldE = 0.d0
        allocate(yieldnE(1:100)); yieldnE = 0.d0
        if(Xe_search)then
        allocate(Xe_cum_yield(1:100)) ! Cumulative yield of Xe135
        allocate(I_cum_yield(1:100)) ! Cumulative yield of Iodine135
        endif
        n_fp_max = 0
        30  continue
        MT = 0
        do while (MT /= 454)
            read(rd_yield,'(A80)',end = 300) line0
            read(line0(73:75),'(i)') MT
        enddo
        fid = fid + 1
        read(line0(1:11),'(f)') ZA
        nuclid = int(ZA)
        anum = nuclid/1000
        mnum = nuclid-anum*1000
        nnum = mnum-anum
        inum = 0
        if(anum==95 .and. mnum==242) inum = 1
        nuclide(inum,nnum,anum)%fiss = .true.
        fssn_anum = anum; fssn_nnum = nnum; fssn_inum = inum
        nuclide(inum,nnum,anum)%fy_idx = fid
        fssn_zai(fid) = nuclid*10+inum
        if(find_ACE_iso_idx_zaid(nuclid*10)>0) then
            ace_fssn(find_ACE_iso_idx_zaid(nuclid*10)) = fid
            !if(icore==score) print *, find_ACE_iso_idx_zaid(nuclid*10), nuclid, fid
        endif
        if(nuclide(inum,nnum,anum)%amu==0.d0) then
            read(line0(12:22),'(f)') nuclide(inum,nnum,anum)%amu
        endif
        read(line0(23:33),'(i)') n_y_eg ! # of energy groups
        allocate(Ep(n_y_eg)) ! Energy points for yield
        !Assign boundary energy
        allocate(nuclide(inum,nnum,anum)%yield_E(n_y_eg))
                do eg = 1,n_y_eg
                    read(rd_yield,'(A80)',end = 300) line0
                    read(line0(1:11),'(f)') Ep(eg)
                    read(line0(56:66),'(i)') n_fp ! # of fission products
                    n_fp_max = max(n_fp_max,n_fp) ! Maximal # of FP for given fssn
                    ! If FPY out of range: constant 
                    do fp = 1,n_fp
                        flag = mod(fp*4-3,6)
                        if (flag<3) read(rd_yield,'(A80)',end=300) line1
                        read(line1(11*flag-10:11*flag),'(f)') ZA
                        read(line1(11*flag+1:11*flag+11),'(f)') i_real ! Flag for ZA & I
                        nuclid = int(ZA)
                        inum = int(i_real)
                        anum = nuclid/1000
                        mnum = nuclid-anum*1000
                        nnum = mnum-anum
                        ! 211109 MODIFICATION: If decay chain of meta. not exists, pass to ground
                        tmpfp = nuclide(inum,nnum,anum)%fp
                        if (nuclide(inum,nnum,anum)%fp == 0) then !Assign ifp to nuclide
                            if(nuclide(inum,nnum,anum)%data_exist) then
                                ifp = ifp + 1
                                nuclide(inum,nnum,anum)%fp = ifp
                                fp_zai(ifp) = nuclid*10+inum
                                tmpfp = nuclide(inum,nnum,anum)%fp
                            endif
                        endif
                        flag = mod(fp*4-1,6) ! Flag for yield
                        if (flag<3) read(rd_yield,'(A80)',end=300) line1
                        read(line1(11*flag-10:11*flag),'(f)') ify
                        if(ify>0.d0 .and. tmpfp>0) ify_yield(tmpfp,eg,fid) = &
                            ify_yield(tmpfp,eg,fid) + ify
                    enddo
                enddo
                nuclide(fssn_inum,fssn_nnum,fssn_anum)%yield_E = Ep*1.d-6
                do eg = 1,n_y_eg
                    yieldE(fid,eg) = Ep(eg)*1.d-6
                enddo
                yieldnE(fid) = n_y_eg
                deallocate(Ep)
                do while (MT /= 0)
                    read(line0(73:75),'(i)') MT
                    read(rd_yield,'(A80)',end=300) line0 ! Skip until end of individ. NFY
                enddo
                ! ====== 211013 ==========
                ! Cumulative NFY should be considered for fpcut operation
                ! fpcut: Neglect cumul. fp < fpcut
                ! May reduce computation burden
                ! TODO
                read(line0(23:33),'(i)') n_y_eg
                do eg = 1,n_y_eg
                    read(rd_yield,'(A80)',end = 300) line1
                    read(line1(56:66),'(i)') n_fp
                    do fp = 1,n_fp
                        flag = mod(fp*4-3,6)
                        if(flag<3) read(rd_yield,'(A80)',end=300) line1
                        read(line1(11*flag-10:11*flag),'(f)') ZA
                        read(line1(11*flag+1:11*flag+11),'(f)') i_real
                        nuclid = int(ZA); inum = 0; inum = int(i_real)
                        anum = nuclid/1000
                        mnum = nuclid-anum*1000
                        nnum = mnum-anum
                        tmpfp = nuclide(inum,nnum,anum)%fp
                        if(tmpfp==0 .and. nuclide(0,nnum,anum)%fp>0) tmpfp = nuclide(0,nnum,anum)%fp
                        flag = mod(fp*4-1,6)
                        if(flag<3) read(rd_yield,'(A80)',end=300) line1
                        read(line1(11*flag-10:11*flag),'(f)') cfy !CUMUL. FY.
                        if(cfy>0.d0 .and. tmpfp>0) cfy_yield(tmpfp,eg,fid) = &
                            cfy_yield(tmpfp,eg,fid) + cfy
                        ! === 211110 ===
                        ! Updates on CFY and IFY -> TMP_FY
                        if(tmpfp>0) then
                        if(cfy_yield(tmpfp,eg,fid)>fpcut) then
                            tmp_yield(tmpfp,eg,fid) = tmp_yield(tmpfp,eg,fid) + ify_yield(tmpfp,eg,fid)
                            nuclide(inum,nnum,anum)%fiss = .true.
                        endif
                        endif
                    enddo
                enddo
                ! =======================
                do while (MT /= 451)
                    read(rd_yield,'(A80)',end = 300) line0 ! Skip until next isotope's description (MT=451)
                    read(line0(73:75),'(i)') MT 
                enddo
                read(rd_yield,*,end=300)
                go to 30
                300 close (rd_yield)
               
                !if(Xe_search) then
                !Xe_source =(/491350,491360,501350,501370,511350,511360,511370,521350,521360,531350/)
                !do i = 1,size(Xe_source)
                !    do j = 1,fid
                !    zai = Xe_source(i)
                !    anum = zai/10000; mnum = (zai-anum*10000)/10; nnum = mnum-anum
                !    inum = zai-anum*10000-mnum*10; fp_tmp = nuclide(inum,nnum,anum)%fp
                !    ! Consider only thermal; Xe-135 impact severly on thermal
                !    if(fp_tmp>0) I_cum_yield(i) = I_cum_yield(i) + tmp_yield(fp_tmp,1,j)
                !    enddo
                !enddo
                !Xe_cum_yield = I_cum_yield
                !do inum = 0,1
                !    do j = 1,fid
                !    anum = 52; mnum = 135; nnum = mnum-anum
                !    fp_tmp = nuclide(inum,nnum,anum)%fp
                !    if(fp_tmp>0) Xe_cum_yield(i) = Xe_cum_yield(i) + tmp_yield(fp_tmp,1,j)
                !    enddo
                !enddo
                !endif
                !
                ! ======= 211116, SFY LIBRARY READ ========
                open(rd_sfy, file = trim(dep_lib)//'sss_endfb71.sfy', action = 'read')
                sfid = 0; isfp = 0
                allocate(sfy_yield(1:2500,1:100)); sfy_yield = 0.d0
                allocate(sfp_zai(1:2500)) ! ZAI of fp from SFY
                allocate(sfssn_zai(1:100)) ! ZAI of fssn from SFY
                allocate(ace_sfssn(1:num_iso)); ace_fssn = 0.d0
                n_sfp_max = 0
                50 continue
                MT = 0
                do while (MT/=454)
                    read(rd_sfy,'(A80)',end=500) line0
                    read(line0(73:75),'(i)') MT
                enddo
                sfid = sfid + 1
                read(line0(1:11),'(f)') ZA; nuclid = int(ZA)
                anum = nuclid/1000; mnum = nuclid-anum*1000
                nnum = mnum-anum; inum =0
                fssn_anum = anum; fssn_nnum = nnum
                nuclide(inum,nnum,anum)%sfy_idx = sfid
                sfssn_zai(sfid) = nuclid*10 + inum
                if(find_ACE_iso_idx_zaid(nuclid*10)>0) ace_sfssn(find_ACE_iso_idx_zaid(nuclid*10)) = sfid
                read(line0(23:33),'(i)') n_y_eg
                if(icore==score .and. n_y_eg>1) print *, 'EXCEPTION in SFY'
                read(rd_sfy,'(A80)',end=500) line0
                read(line0(56:66),'(i)') n_fp
                n_sfp_max = max(n_sfp_max,n_fp)
                do fp = 1,n_fp
                    flag = mod(fp*4-3,6)
                    if(flag<3) read(rd_sfy,'(A80)',end=500) line1
                    read(line1(11*flag-10:11*flag),'(f)') ZA
                    read(line1(11*flag+1:11*flag+11),'(f)') i_real
                    nuclid = int(ZA); inum = int(i_real)
                    anum = nuclid/1000; mnum = nuclid-anum*1000; nnum = mnum-anum
                    tmpfp = nuclide(inum,nnum,anum)%sfp
                    if(nuclide(inum,nnum,anum)%sfp ==0 .and. nuclide(inum,nnum,anum)%data_exist) then
                        isfp = isfp + 1
                        nuclide(inum,nnum,anum)%sfp = isfp
                        sfp_zai(isfp) = nuclid*10+inum
                        tmpfp = nuclide(inum,nnum,anum)%sfp
                    endif
                    flag = mod(fp*4-1,6)
                    if(flag<3) read(rd_sfy,'(A80)',end=500) line1
                    read(line1(11*flag-10:11*flag),'(f)') sfy
                    if(sfy>0.d0 .and. tmpfp>0) sfy_yield(tmpfp,sfid) = &
                        sfy_yield(tmpfp,sfid) + sfy
                enddo
                
                do while (MT/=451)
                    read(line0(73:75),'(i)') MT
                    read(rd_sfy,'(A80)',end=500) line0
                enddo
                read(rd_sfy,*,end=500)
                go to 50
                500 close(rd_sfy)



                !========== TESTING y_cumul CUTOFF ================
                ! If HL > CUTOFF, remove the nuclide
                !=============================================================================
                !DECAY CHAIN REDUCTION PROCESS
                if(.not. chain_reduction) go to 11
                do inum = 0,1
                do nnum = 0,170
                do anum = 0,111
                    ! Skip non-decay defined nuclides
                    if(.not. nuclide(inum,nnum,anum)%data_exist) cycle
                    ! Reduce until stable-enough (HL>HL_th)
                    lam_th = log(2.d0)/HL_th
                    n_dau = nuclide(inum,nnum,anum)%react_num
                    if(n_dau==0) cycle !Escape if no daughter nuclides
                    i = 1
                    REDU: do
                        prod_red = nuclide(inum,nnum,anum)%daughter(i,:)
                        f_tmp = nuclide(inum,nnum,anum)%frac(i)
                        lam = nuclide(prod_red(1),prod_red(2),prod_red(3))%lambda
                        if(lam>lam_th)then !If product is unstable enough
                            nuclide(prod_red(1),prod_red(2),prod_red(3))%reduced = .true.                    
                            ! Combine %REACT_NUM, %FRAC and %DAUGHTER
                            ! 1. Recalculate %REACT_NUM
                            n_grand = nuclide(prod_red(1),prod_red(2),prod_red(3))%react_num
                            tmp_n = n_dau
                            n_dau = n_dau + n_grand - 1 !itself
                            nuclide(inum,nnum,anum)%react_num = n_dau
                            
                            ! 2. Recalculate %FRAC
                            frac_grand =  nuclide(prod_red(1),prod_red(2),prod_red(3))%frac
                            frac_dau = nuclide(inum,nnum,anum)%frac
                            deallocate(nuclide(inum,nnum,anum)%frac)
                            allocate(nuclide(inum,nnum,anum)%frac(n_dau))
                            if(i>1) nuclide(inum,nnum,anum)%frac(1:i-1)= frac_dau(1:i-1)
                            if(i<tmp_n) nuclide(inum,nnum,anum)%frac(i:tmp_n-1) = frac_dau(i+1:tmp_n)
                            nuclide(inum,nnum,anum)%frac(tmp_n:n_dau) = frac_grand*frac_dau(i)

                            ! 3. Combine %DAUGHTER
                            ! Store previous daughters
                            prod_grand = nuclide(prod_red(1),prod_red(2),prod_red(3))%daughter(:,:)
                            prod_dau = nuclide(inum,nnum,anum)%daughter(:,:)
                            ! Combine two stored daughters
                            deallocate(nuclide(inum,nnum,anum)%daughter)
                            allocate(nuclide(inum,nnum,anum)%daughter(n_dau,3))
                            if(i>1) nuclide(inum,nnum,anum)%daughter(1:i-1,:) = prod_dau(1:i-1,:)
                            if(i<tmp_n) nuclide(inum,nnum,anum)%daughter(i:tmp_n-1,:) = prod_dau(i+1:tmp_n,:) 
                            nuclide(inum,nnum,anum)%daughter(tmp_n:n_dau,:) = prod_grand

                            ! 4. Assign FPY
                            !fp_dau = nuclide(prod_red(1),prod_red(2),prod_red(3))%fp
                            !if(fp_dau>0) then
                            !do k = 1,n_grand ! For every daughter nuclides
                            !    frac_tmp = frac_grand(k)
                            !    prod_tmp = prod_grand(k,:) ! Decay fraction and products
                            !    fp_tmp = nuclide(prod_tmp(1),prod_tmp(2),prod_tmp(3))%fp !FP value for daugther k
                            !    ! FPY of daughter += frac_tmp * yield
                            !    if(fp_tmp==0) then !daughter is not FP
                            !        fid = fid + 1
                            !        nuclide(prod_tmp(1),prod_tmp(2),prod_tmp(3))%fp = fid
                            !        fp_zai(fid) = 10010*prod_tmp(3)+10*prod_tmp(2)+prod_tmp(1)
                            !        tmp_yield(fid,:,:) = tmp_yield(fp_dau,:,:)*frac_tmp
                            !    else !Additional FPY
                            !        tmp_yield(fp_tmp,:,:) = tmp_yield(fp_tmp,:,:) + frac_tmp*tmp_yield(fp_dau,:,:)
                            !    endif
                            !enddo
                            !tmp_yield(fp_dau,:,:) = 0.d0
                            !endif
                        else !ith daughter nuclide is stable now
                            i = i + 1 ! Search for branches of next particles
                            if(i>n_dau) exit REDU
                        endif
                        ! Assign product alpha/proton/neutron emission to its mother
                        nuclide(inum,nnum,anum)%a_emit = nuclide(inum,nnum,anum)%a_emit&
                        +nuclide(prod_red(1),prod_red(2),prod_red(3))%a_emit*f_tmp
                        nuclide(inum,nnum,anum)%n_emit = nuclide(inum,nnum,anum)%n_emit&
                        +nuclide(prod_red(1),prod_red(2),prod_red(3))%n_emit*f_tmp
                        nuclide(inum,nnum,anum)%p_emit = nuclide(inum,nnum,anum)%p_emit&
                        +nuclide(prod_red(1),prod_red(2),prod_red(3))%p_emit*f_tmp
                    enddo REDU
                enddo
                enddo
                enddo
                do inum = 0,2
                do nnum = 0,170
                do anum = 0,111
                    !if(nuclide(inum,nnum,anum)%data_exist .and. nuclide(inum,nnum,anum)%reduced)&
                    !nuclide(inum,nnum,anum)%data_exist = .false.
                enddo
                enddo
                enddo
                ! If nuclide reduced but exists in 'inventory.inp':
                do i = 1,num_iso
                zai = ace(i)%zaid
                anum = zai/1000; mnum = zai-anum*1000
                if(mnum>300) then
                    inum = 1
                    mnum = mnum - 200
                endif
                nnum = mnum - anum
                !if(.not. nuclide(inum,nnum,anum)%data_exist .and.  nuclide(inum,nnum,anum)%reduced)&
                !nuclide(inum,nnum,anum)%data_exist = .true.
                enddo
                11 continue
                !==============================================================================
                ! CROP NFY DATA
                tmp_yield = tmp_yield(1:ifp,1:4,1:fid)
                deallocate(ify_yield); deallocate(cfy_yield)
                fp_zai = fp_zai(1:ifp)
                fssn_zai = fssn_zai(1:fid)
                yieldE = yieldE(1:fid,1:4)
                yieldnE = yieldnE(1:fid)
                nfp = ifp
                nfssn = fid
                allocate(fratio(nfssn,4))
                if(Xe_search) then
                I_cum_yield = I_cum_yield(1:nfssn)
                Xe_cum_yield = Xe_cum_yield(1:nfssn) 
                endif
                !=============================================================================
                ! CROP SFY DATA
                sfy_yield = sfy_yield(1:isfp,1:sfid)
                sfp_zai = sfp_zai(1:isfp)
                sfssn_zai = sfssn_zai(1:sfid)
                nsfp = isfp; nsfssn = sfid
                
                ! =============================
                ! REMOVE USELESS NUCLIDES
                !call reducenuc

                !============================================================================== 
                !Sort nuclide data
                nnuc = 0
                do anum=0,111
                do nnum=0,170
                do inum=0,2
                   if(nuclide(inum,nnum,anum)%data_exist) then
!                  if(nuclide(inum,nnum,anum)%data_exist .and. nuclide(inum,nnum,anum)%fiss) then
                    nnuc = nnuc + 1
                    nuclide(inum,nnum,anum)%idx = nnuc
                    zai_idx(nnuc) = anum*10000+(nnum+anum)*10+inum
                  end if
                end do
                end do
                end do
                 
                do i = 1,num_iso
                    zai = ace(i)%zaid
                    anum = zai/1000; mnum = zai-anum*1000
                    inum = 0
                    if(mnum>300) then
                        mnum = mnum-200
                        if(anum>88) mnum = mnum + 100
                        inum = 1
                    endif
                    nnum = mnum-anum
                    if(nuclide(inum,nnum,anum)%idx==0) then
                        nnuc = nnuc + 1
                        nuclide(inum,nnum,anum)%idx = nnuc
                        zai_idx(nnuc) = anum*10000+mnum*10+inum
                    endif
                 enddo
                maxnnz = int(1 * nnuc**2) ! Assuming less than 10% sparsity  
                allocate(Acsr(1:maxnnz))
                allocate(jAcsr(1:maxnnz))
                allocate(iAcsr(1:nnuc+1))
                if(icore==score) print *, "Number of isotopes:", nnuc-1 ! Checking number of isotopes
                if(icore==score) print *, "ACE FORMAT ISOTOPES:", num_iso 

                !==================================================================================
                !ISOMERIC BRANCHING RATIO (UNDER TESTING, 10/18)
                open(isom_brn,file=trim(dep_lib)//'isomeric_branching',action='read')
                read(isom_brn,'(A80)',end = 400) line0
                read(line0,'(i)') num_brn

                allocate(ZAIMT_ism(num_brn))
                allocate(gnd_frac(num_brn))

                do brn = 1,num_brn
                   read(isom_brn,'(A80)',end = 400) line0
                   read(line0(1:6),'(i)') ZAIMT_ism(brn)
                   read(line0(13:20),'(f)') gnd_frac(brn)
                   !TESTING
                   !if(icore==score) print *, 'ISOM',ZAIMT_ism(brn),gnd_frac(brn)
                enddo
400             close(isom_brn)
                !Burnup result files
                if(icore==score) then
                    !open(prt_ntpy, file="dep_monitor",action="write",status="replace") 
                    open(prt_bumat, file="bumat.out",action="write",status="replace") !position='append')
                    !Initial Density
                    if(istep_burnup==0) then ! Initial Number
                    write(prt_bumat, '(a45)')         '   =========================================='
                    write(prt_bumat, '(a17,i4)')     '      Burnup step', istep_burnup
                    write(prt_bumat, '(f14.2,a16)') burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
                    write(prt_bumat, '(a45)')         '   =========================================='
                    do i = 1,n_materials
                    ! Initial mass
                        do mt_iso = 1,materials(i)%n_iso
                            iso = materials(i)%ace_idx(mt_iso)
                            zai = ace(iso)%zaid
                            if(zai>0) then
                            anum = zai/1000; mnum = (zai-anum*1000)
                            nnum = mnum-anum; inum = zai-anum*1000-mnum
                            !$omp atomic
                            tot_mass_init = tot_mass_init + &
                            ace(iso)%atn* m_n * materials(i)%numden(mt_iso) &
                            * materials(i)%vol/N_AVOGADRO
                            if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
                            !$omp atomic
                            tot_fmass_init = tot_fmass_init + &
                            ace(iso)%atn * m_n * materials(i)%numden(mt_iso) &
                            * materials(i)%vol/N_AVOGADRO
                            endif
                            endif
                        enddo
                        if(.not. materials(i)%depletable) cycle
                        write(prt_bumat,*) ' '
                        write(prt_bumat,*) materials(i)%mat_name
                        write(prt_bumat,*) ' '
                        do mt_iso = 1,materials(i)%n_iso
                            write(prt_bumat,'(a15,e14.5)')&
                            ace(materials(i)%ace_idx(mt_iso))%library,&
                            materials(i)%numden(mt_iso)*barn
                        enddo
                        write(prt_bumat,*) 'Num isotope', materials(i)%n_iso
                        write(prt_bumat,*) ' '
                    enddo
                    write(prt_bumat,*) 'Total mass[g]:', tot_mass_init
                    write(prt_bumat,*) 'Fiss. mass[g]:', tot_fmass_init
                    endif
                end if
                if(icore==score) print *, 'BUMAT WRITTEN'

            end subroutine getENDFdepletionlibrary
            
            function E2G(ee,Etmp)
            integer :: E2G
            real(8),intent(in) :: ee
            real(8),intent(in) :: Etmp(:)
            integer :: ii, n_E

            n_E = size(Etmp)
            if(n_E==1) then
                E2G = 1
                return
            endif
            
            E2G = n_E+1
            do ii = 1,n_E
            if (ee<=Etmp(ii)) then
                E2G = ii; return
            endif
            enddo
            end function
            
            function unitdigit(n)
            integer :: unitdigit, n, tmp
            if(n<10) then
                unitdigit = n; return
            endif
            tmp = n/10
            unitdigit = n-tmp*10
            end function
            
            subroutine reducenuc
            !integer :: zaid
            ! STEP 0. Assign ZAID of ACE formatted nuclides
            !allocate(nuc(1:4000)); nuc = 0
            !do iso = 1,num_iso
            !    zai = ace(iso)%zaid/1000
            !    anum = zai/1000; mnum = zai-anum*1000
            !    inum = 0
            !    if(mnum>300) then
            !       inum = 1
            !       mnum = mnum - 200
            !       if(anum>88) mnum = mnum + 100
            !   endif
            !   zaid = anum*10000+mnum*10+inum
            !   nuc(iso) = zaid
            !enddo

            !! RECURRENCE read chain
            !flag = 0
            !do
            !flag = flag + 1
            !zaid = nuc(flag)
            !anum = zaid/10000; mnum = (zaid-anum*10000)/10
            !nnum = mnum - anum; inum = zaid-anum*10000-mnum*10
            !! PATH 1: TRANSMUTATION
            !! ...
            !! PATH 2: DECAY CHAIN
            !ndec = nuclide(inum,nnum,anum)%react_num
            !if(ndec>0) then
            !    do d = 1,ndec
            !        inum1 = nuclide(inum,nnum,anum)%daguther(d,1)
            !        nnum1 = nuclide(inum,nnum,anum)%daughter(d,2)
            !        anum1 = nuclide(inum,nnum,anum)%daughter(d,3)
            !    enddo
            end subroutine

            ! ===================================================================
            !         Depletion :: Make depletion matrix and solve 
            ! ===================================================================
            subroutine depletion 
            
            implicit none
            
            integer :: imat, jmem, jnuc, knuc, knuc1, knuc2
            integer :: mt_iso, iso
            integer :: anum, mnum, nnum, inum
            integer :: anum1, mnum1, nnum1, inum1
            integer :: i, j, k
            logical :: exist_in_MC
            real(8) :: real_flux, tot_flux=0.d0
            integer :: fy_midx, diff, tmp, ierr
            real(8) :: nxt_full_numden(1:nnuc)
            type (Material_CE), pointer :: mat
            integer :: iso_idx(nnuc)
            integer :: a1, m1,n1
            real(8), allocatable :: frac(:)
            integer, allocatable :: prod(:,:)
            real(8), allocatable :: f_decay(:)
            real(8), allocatable :: nucexist(:)
            
            real(8) :: t_circulation = 10.0 ! Circulation time

            real(8):: Ep(4)
            real(8), allocatable :: ftmp(:)
            integer :: g,n_E,nE
            integer :: zai
            integer :: ism !TEST 211018
            real(8) :: fis_ratio(3)
            real(8) :: ratio
            real(8) :: remsum
            ! TESTING for SM-149
            real(8) :: samarium !Estimate Samarium
            real(8) :: samabs
            real(8) :: ingrid
            character(len=10) :: fileid, matid
            character(len=50) :: filename

            ! Homogenized OGXS
            real(8) :: sigmaa
            real(8) :: nusigf
            integer :: ii

            if(do_burn==.false.) return

            avg_power = avg_power / dble(n_act)
            avg_fiss = avg_fiss / dble(n_act)
            tot_flux = 0.d0

            call MPI_BCAST(avg_power, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(avg_fiss, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)
            if(icore==score) print *, 'Avg power[MeV]',avg_power 
            if(icore==score) print *, 'KAPPA[MeV]', avg_power/avg_fiss
            ! Interpolation for yield with given avg_Ef
            fis_ratio = (/fis_thermal,fis_epi,fis_fast/)
            if(icore==score) print *, 'E DIST',fis_ratio
            allocate(yield_data(nfp,nfssn))
            yield_data = 0.d0
            do i = 1,nfssn !fissionables
                !if(sum(fratio(i,:))==0.d0) then
                !    fis_ratio = (/fis_thermal,fis_epi,fis_fast/)
                !    !write(*,'(i,f5.3,f5.3,f5.3)') fssn_zai(i),fis_ratio(1)/sum(fis_ratio(:)),fis_ratio(2)/sum(fis_ratio(:)),fis_ratio(3)/sum(fis_ratio(:))
                !else
                !    fis_ratio = fratio(i,1:3)
                !    if(icore==score) &
                !    write(*,'(i,f5.3,f5.3,f5.3)') fssn_zai(i),fratio(i,1)/sum(fratio(i,:)),fratio(i,2)/sum(fratio(i,:)),fratio(i,3)/sum(fratio(i,:))
                !endif
                zai = fssn_zai(i)
                anum = zai/10000; mnum = (zai-10000*anum)/10; nnum = mnum-anum
                inum = zai-anum*10000-mnum*10
                if(allocated(ftmp)) deallocate(ftmp)
                allocate(ftmp(min(3,n_E))); ftmp = 0.d0
                Ep = yieldE(i,1:4); nE = yieldnE(i)
                if(nE==0) then !ASSUME THERMAL
                    yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i)
                else
                    if(icore==score) print *, i, fssn_zai(i)
                    if(sum(fratio(i,1:nE))>0) then
                        fratio(i,1:nE) = fratio(i,1:nE)/sum(fratio(i,1:nE))
                    else
                        fratio(i,1) = 1.d0
                        fratio(i,2:4) = 0.d0
                    endif
                    do k = 1,nE
                        !if(icore==score) print *, k,Ep(k), fratio(i,k)
                        yield_data(1:nfp,i) = yield_data(1:nfp,i) + &
                            tmp_yield(1:nfp,k,i)*fratio(i,k)
                    enddo
                endif
                !ingrid = 1.d0
                !do k = 1,min(3,n_E)
                !   if(Ep(k)<thermal_ub) then
                !    ftmp(k) = ftmp(k) + fis_ratio(1)
                !   elseif(Ep(k)>fast_lb) then
                !    ftmp(k) = ftmp(k) + fis_ratio(3)
                !   else
                !    ftmp(k) = ftmp(k) + fis_ratio(2)
                !   endif
                !   if(k>1) then
                !    if(ftmp(k)==ftmp(k-1)) then
                !     ftmp(k) = ftmp(k)*ingrid/(ingrid+1.d0)
                !     ftmp(k-1) = ftmp(k-1)*ingrid/(ingrid+1.d0)
                !     ingrid = ingrid + 1.d0
                !    endif
                !    ingrid = 1.d0
                !   endif
                !enddo
                !if(ANY(Ep<thermal_ub)) ftmp(1) = fis_ratio(1)
                !if(ANY(Ep>fast_lb)) ftmp(3) = fis_ratio(3)
                !if(ANY(Ep>=thermal_ub .and. Ep<=fast_lb)) ftmp(2) = fis_ratio(2)
                !if(sum(ftmp)==0.d0)then
                !    yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i) !ASSUME THERMAL
                !else
                !if(icore==score) print *, 'FTMP',fssn_zai(i),Ep,ftmp
                !do k = 1,min(3,n_E)
                    !if(Ep(k)<thermal_ub)then
                    !    yield_data(1:nfp,i) = yield_data(1:nfp,i) + tmp_yield(1:nfp,k,i)*ftmp(1)
                    !elseif(Ep(k)>fast_lb)then
                    !    yield_data(1:nfp,i) = yield_data(1:nfp,i) + tmp_yield(1:nfp,k,i)*ftmp(3)
                    !else
                    !    yield_data(1:nfp,i) = yield_data(1:nfp,i) + tmp_yield(1:nfp,k,i)*ftmp(2)
                    !endif
                    !yield_data(1:nfp,i) = yield_data(1:nfp,i) + tmp_yield(1:nfp,k,i)*ftmp(k)
                !enddo
                !endif
            enddo

            allocate(nucexist(nfp)); nucexist = 0.d0
            do i = 1,nfp
                ! PRIOR to normalization; exclude non-existing nuclides in FPY/SFY
                !zai = fp_zai(i)
                !anum = zai/10000; mnum = (zai-anum*10000)/10; nnum = mnum-anum;
                !inum = zai-anum*10000-mnum*10;
                !if(nuclide(inum,nnum,anum)%idx>0) nucexist(i) = 1.d0
            enddo

            do i = 1,nfssn
                ! NORMALIZE NFY: sum(NFY) = 200
                !ratio = sum(yield_data(:,i)*nucexist)/2.d0
                !if(ratio>0.d0) yield_data(:,i) = yield_data(:,i)/ratio
                !if(icore==score) print *, 'NFY', fssn_zai(i), ratio
            enddo
            deallocate(nucexist)
            
            ! NORMALIZE SFY too
            allocate(nucexist(nsfp)); nucexist = 0.d0
            do i = 1,nsfp
                !zai = sfp_zai(i)
                !anum = zai/10000; mnum = (zai-anum*10000)/10; nnum = mnum-anum;
                !inum = zai-anum*10000-mnum*10;
                !if(nuclide(inum,nnum,anum)%idx>0) nucexist(i) = 1.d0
            enddo
            do i = 1,nsfssn
                !ratio = sum(sfy_yield(:,i)*nucexist)/2.d0
                !if(ratio>0.d0) sfy_yield(:,i) = sfy_yield(:,i)/ratio
                !if(icore==score) print *, 'SFY', sfssn_zai(i), ratio
            enddo
            deallocate(nucexist)

            !Initialize material independent burnup matrix
            if(istep_burnup==0) then
            allocate(bMat(1:nnuc,1:nnuc))   !2-D burnup matrix : row to column transition
            allocate(bMat0(1:nnuc,1:nnuc)) !2-D burnup matrix : row to column transition (material independent)
            
            bMat0 = 0.d0
            
            
            !Build material independent burnup matrix
            do jnuc = 1, nnuc
                anum = zai_idx(jnuc)/10000
                mnum = (zai_idx(jnuc) - anum*10000)/10
                nnum = mnum - anum
                inum = zai_idx(jnuc) - anum*10000 - mnum*10
                if(nnum<0) print *, zai_idx(jnuc)
                if(nuclide(inum,nnum,anum)%react_num==0) cycle
                allocate(f_decay(nuclide(inum,nnum,anum)%react_num))
                f_decay = nuclide(inum,nnum,anum)%frac
                allocate(prod(nuclide(inum,nnum,anum)%react_num,3))
                prod = nuclide(inum,nnum,anum)%daughter
                bMat0(jnuc,jnuc) = bMat0(jnuc,jnuc) - nuclide(inum,nnum,anum)%lambda
                do k = 1,nuclide(inum,nnum,anum)%react_num
                ! In case of non-fission
                if (prod(k,1)>0 .or. prod(k,2)>0 .or. prod(k,3)>0) then
                    knuc = nuclide(prod(k,1),prod(k,2),prod(k,3))%idx
                    if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + nuclide(inum,nnum,anum)%lambda*f_decay(k)
                elseif(nuclide(inum,nnum,anum)%sfiss==k) then !Spontaneous fission
                    if(icore==score) print *, zai_idx(jnuc), nuclide(inum,nnum,anum)%frac, prod(k,:)
                    if(nuclide(inum,nnum,anum)%sfy_idx>0) then
                        ! If SFY exists
                        do j=1,nsfp
                            anum1 = sfp_zai(j)/10000
                            mnum1 = (sfp_zai(j) - anum1*10000)/10
                            nnum1 = mnum1 - anum1
                            inum1 = sfp_zai(j) - anum1*10000 - mnum1*10
                            knuc = nuclide(inum1,nnum1,anum1)%idx
                            if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) &
                                + nuclide(inum,nnum,anum)%lambda*f_decay(k) &
                                *sfy_yield(j,nuclide(inum,nnum,anum)%sfy_idx)
                        end do
                    else
                        if(anum>=89) then !For actinides
                        ! When fission yield not exist while it fissions;
                        !1. check for its ground state
                        !2. check for its isotone
                        !3. use U-235
                        ! Check for ground state
                        fy_midx = 0
                        if(inum>=1 .and. nuclide(0,nnum,anum)%sfy_idx>0) then
                            fy_midx = nuclide(0,nnum,anum)%sfy_idx
                        else
                        ! Check for isotones (first found)
                        do i = 1,nfssn
                            a1 = sfssn_zai(i)/10000
                            m1 = (sfssn_zai(i)-a1*10000)/10
                            n1 = m1-a1
                            if (n1==nnum) then
                                fy_midx = i
                            endif
                        enddo
                        if(fy_midx==0) fy_midx = nuclide(0,146,92)%sfy_idx
                        ! Use U-238 otherwise
                        
                        do j=1,nsfp
                            anum1 = sfp_zai(j)/10000
                            mnum1 = (sfp_zai(j) - anum1*10000)/10
                            nnum1 = mnum1 - anum1
                            inum1 = sfp_zai(j) - anum1*10000 - mnum1*10
                            knuc = nuclide(inum1,nnum1,anum1)%idx
                            if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
                                nuclide(inum,nnum,anum)%lambda*f_decay(k) &
                                *sfy_yield(j,fy_midx)
                        end do
                        end if
                        end if
                    endif
                endif
                enddo
                deallocate(f_decay)
                deallocate(prod)
                ! == 211111 update: Ignore emit from decay? ===
                ! Count for emissions; alpha, neutron and proton
                if (nuclide(inum,nnum,anum)%a_emit>0.d0) then
                    knuc = nuclide(0,2,2)%idx
                    !if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
                    !    nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%a_emit
                endif
                if (nuclide(inum,nnum,anum)%n_emit>0.d0) then
                    knuc = nuclide(0,1,0)%idx
                    !if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
                    !    nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%n_emit
                endif
                if (nuclide(inum,nnum,anum)%p_emit>0.d0) then
                    knuc = nuclide(0,0,1)%idx
                    !if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
                    !    nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%p_emit
                endif
                
                ! Count for removal per circulation
                if (nuclide(inum,nnum,anum)%removal > 0.d0) then
                    bMat0(jnuc,jnuc) = bMat0(jnuc,jnuc) + &
                    log(1.d0-nuclide(inum,nnum,anum)%removal)/t_circulation
                end if
            end do
            end if
            
            if (icore==score) then 
                write(prt_bumat, '(a45)')         '   =========================================='
                write(prt_bumat, '(a17,i4)')     '      Burnup step', istep_burnup+1
                write(prt_bumat, '(f14.2,a16)') burn_step(istep_burnup+1)/86400.d0, ' CUMULATIVE DAYS'
                write(prt_bumat, '(a45)')         '   =========================================='
            endif 
            
            !Normalization constant to be real power
            ULnorm = RealPower/(avg_power*eVtoJoule)
            
            !Substitute burnup matrix element
            do imat = 1, n_materials
                if(.not. materials(imat)%depletable) cycle    !material imat is not burned
                !samarium = 0.d0
                mat => materials(imat)
                if (icore==score) then 
                    write(prt_bumat, *) '' 
                    write(prt_bumat, *) mat%mat_name 
                    write(prt_bumat, *) ''
                endif         
                !Call the material independent burnup matrix 
                bMat = bMat0*bstep_size
                
                !Calculate real flux (volume-averaged)
                real_flux = ULnorm*mat%flux
                tot_flux = tot_flux + real_flux*mat%vol
                !Build burnup matrix with cross section obtained from MC calculation
                !DO_ISO: do mt_iso = 1, mat%n_iso
                    !iso = mat%ace_idx(mt_iso)
                DO_ISO: do mt_iso = 1,num_iso
                    iso = mt_iso
                    anum = ace(iso)%zaid/1000
                    mnum = (ace(iso)%zaid - anum*1000)
                    inum = 0
                    if(mnum>300) then
                        inum = 1
                        mnum = mnum - 200
                        if(anum>88) mnum = mnum + 100
                        !if(icore==score) print *,'META',ace(iso)%zaid,anum,mnum,inum
                    endif
                    nnum = mnum-anum
                    jnuc = nuclide(inum,nnum,anum)%idx
                    if(jnuc==0) cycle
                    
                    ! OGXS modification: 211117
                    ! OGXS 1: (n,g) reaction (branching included)
                    zai = ace(iso)%zaid*10+inum
                    ism = 0
                    do i = 1,num_brn
                        if(zai==ZAIMT_ism(i)) then
                            ism = i
                            exit
                        endif
                    enddo
                    if(ism>0) then   
                        knuc = nuclide(0,nnum+1,anum)%idx
                        if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,1)*bstep_size*gnd_frac(ism) !To GROUND state
                        knuc = nuclide(1,nnum+1,anum)%idx
                        if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,1)*bstep_size*(1.d0-gnd_frac(ism)) !To ISOMERIC state
                    else
                        knuc = nuclide(0,nnum+1,anum)%idx
                        if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,1)*bstep_size
                    endif

                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,1)*bstep_size
                    
                    ! OGXS 2: (n,2n) reaction: if nnum < 1, cycle
                    if(nnum<1) goto 920
                    knuc = nuclide(0,nnum-1,anum)%idx
                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,2)*bstep_size

                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,2)*bstep_size

                    920 continue

                    ! OGXS 3: (n,3n) reaction: if nnum < 2, cycle
                    if(nnum<2) goto 930
                    knuc = nuclide(0,nnum-2,anum)%idx
                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,3)*bstep_size
                    
                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,3)*bstep_size

                    930 continue

                    ! OGXS 4: (n,4n) reaction: if nnum < 3, cycle
                    if(nnum<3) goto 940
                    knuc = nuclide(0,nnum-3,anum)%idx
                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,4)*bstep_size

                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,4)*bstep_size

                    940 continue

                    ! OGXS 5: (n,p) reaction: if anum <1, cycle (but not exists...)
                    if(anum<1) goto 950
                    knuc = nuclide(0,nnum+1,anum-1)%idx
                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,5)*bstep_size

                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,5)*bstep_size

                    950 continue

                    ! OGXS 6: (n,a) reaction: if anum<2 or nnum<1, cycle
                    if(anum<2 .or. nnum<1) goto 960
                    knuc = nuclide(0,nnum-1,anum-2)%idx
                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,6)*bstep_size

                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,5)*bstep_size

                    960 continue

                    ! OGXS 7: (n,f) reaction: if not fissionable, cycle
                    if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then !Fissionable
                        if(nuclide(inum,nnum,anum)%fy_idx>0) then   
                            fy_midx = nuclide(inum,nnum,anum)%fy_idx
                        else
                            if(anum<89) cycle
                        ! When fission yield not exist while it fissions;
                        !1. check for its ground state
                        !2. check for its isotone
                        !3. use U-235
                        
                        !Case1.  Check for ground state
                        fy_midx = 0
                        if(inum==1 .and. nuclide(0,nnum,anum)%fy_idx>0) then
                            fy_midx = nuclide(0,nnum,anum)%fy_idx
                        else
                        !Case2.  Check for isotones (first found)
                        do i = 1,nfssn
                            a1 = fssn_zai(i)/10000
                            m1 = (fssn_zai(i)-a1*10000)/10
                            if (m1==mnum) fy_midx = i
                        enddo

                        !Case3. Use U-235
                        if(fy_midx==0) fy_midx = nuclide(0,143,92)%fy_idx
                        endif
                        endif
                        ! NFP calc.
                        do j=1,nfp
                            anum1 = fp_zai(j)/10000
                            mnum1 = (fp_zai(j) - anum1*10000)/10
                            nnum1 = mnum1 - anum1
                            inum1 = fp_zai(j) - anum1*10000 - mnum1*10
                            knuc = nuclide(inum1,nnum1,anum1)%idx
                            if(knuc/=0) bMat(knuc,jnuc) = &
                                bMat(knuc,jnuc) + real_flux*mat%ogxs(mt_iso,7)*yield_data(j,fy_midx)*bstep_size
                        end do
                        
                    bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,7)*bstep_size !reduction of jnuc
                    
                    ! IF REFUELS
                    if(refuel) then
                    do i = 1,n_rf
                        knuc = nuclide(0,mnum_rf(i)-anum_rf(i),anum_rf(i))%idx
                        if(knuc>0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                        real_flux*mat%ogxs(mt_iso,7)*r_rf(i)*bstep_size
                    enddo
                    endif
                  end if
        end do DO_ISO
        ! WRITE BURNUP MATRIX (OPTIONAL)
        write(fileid,'(i2)') istep_burnup
        filename = 'BUMAT_mat_'//trim(adjustl(mat%mat_name(:)))//'_step_'//trim(adjustl(fileid))//'.m'
        open(bumat_test, file = trim(filename),action="write",status="replace")
        write(bumat_test,*) 'FLUX=',real_flux,';'
        write(bumat_test,*) 'ZAI1=zeros(1,',nnuc,');'
        do knuc = 1,nnuc
            write(bumat_test,*) 'ZAI1(',knuc,')=',zai_idx(knuc),';'
        enddo
        write(bumat_test,*) 'A1 = zeros(',nnuc,',',nnuc,');'
        do knuc1 = 1,nnuc
            do knuc2 = 1,nnuc
                if(bMat(knuc1,knuc2)/=0.d0) write(bumat_test,*) 'A1(',knuc1,',',knuc2,')=',bMat(knuc1,knuc2),';'
            enddo
        enddo

        !Solve the burnup matrix equation 
        select case(matrix_exponential_solver)
        case(0)
            call cram(bMat, materials(imat)%full_numden(:), nxt_full_numden(:)) 
            !allocate(bMat_tmp(1:nnuc,1:nnuc)) 
            !call r8mat_expm1 ( nnuc, bMat, bMat_tmp )
            !call r8mat_expm2 ( nnuc, bMat, bMat_tmp )
            !nxt_full_numden = matmul(bMat_tmp,materials(imat)%full_numden)
            !deallocate(bMat_tmp)
            do jnuc=1, nnuc
                if(nxt_full_numden(jnuc)<0.d0) nxt_full_numden(jnuc) = 0.d0
            end do
            
        case default 
            print *, "ERROR :: No such matrix_exponential_solver option", matrix_exponential_solver
            stop 
        end select
            
        
        ! WRITE ATOMIC DENSITY IN MATLAB .m FILE (OPTIONAL)
        write(bumat_test,*) 'Ncomp=zeros(',nnuc,',2);'
        do knuc = 1,nnuc
            if(mat%full_numden(knuc)>0) write(bumat_test,*) 'Ncomp(',knuc,',1)=',mat%full_numden(knuc)*barn,';'
            if(nxt_full_numden(knuc)>0) write(bumat_test,*) 'Ncomp(',knuc,',2)=',nxt_full_numden(knuc)*barn,';'
        enddo
            
        !Update number density
        mat%full_numden = nxt_full_numden 
        
        ! ======================================================================================
        ! Reset material isotope inventory
        knuc = 0 
        remsum = 0.d0
        do jnuc=1, nnuc 
            tmp = find_ACE_iso_idx_zaid(zai_idx(jnuc))
            if(mat%full_numden(jnuc)>0.d0 .and. tmp > 0) then 
                knuc = knuc + 1
                iso_idx(knuc) = jnuc
            elseif(mat%full_numden(jnuc)>0.d0) then
                ! MEASURING LOST TERM
                remsum = remsum + mat%full_numden(jnuc)
            endif
        end do

        mat%n_iso = knuc
        deallocate(mat%ace_idx); allocate(mat%ace_idx(1:knuc))
        deallocate(mat%numden); allocate(mat%numden(1:knuc))
        i = 0 
        do mt_iso=1, mat%n_iso
            ! find ace_idx
            tmp = find_ACE_iso_idx_zaid(zai_idx(iso_idx(mt_iso)))
            if (tmp /= 0 ) then 
                i = i + 1
                mat%ace_idx(mt_iso) = tmp
                mat%numden(mt_iso)  = mat%full_numden(iso_idx(mt_iso))
                
                if (icore==score) then 
                    !print '(i3,i10,a2, a15,e14.5)', &
                    !    i, zai_idx(iso_idx(mt_iso)), '  ',&
                    !    ace(mat%ace_idx(mt_iso))%library, mat%numden(mt_iso)
                    ! ADENS
                    write(prt_bumat, '(a15,e14.5,e14.5,e14.5)') ace(mat%ace_idx(mt_iso))%library, mat%numden(mt_iso)*barn, mat%numden(mt_iso)*barn*N_AVOGADRO*sum(mat%ogxs(mt_iso,1:7)), mat%numden(mt_iso)*barn*N_AVOGADRO*mat%ogxs(mt_iso,7)
                    ! A
                    !write(prt_bumat, '(a15,e14.5)') ace(mat%ace_idx(mt_iso))%library, mat%numden(mt_iso)*mat%vol*barn
                    zai = zai_idx(iso_idx(mt_iso))
                    anum = zai/10000; mnum = (zai-anum*10000)/10
                    nnum = mnum-anum; inum = zai-anum*10000-mnum*10
                    !if(zai==611490) print *, 'Pm-149',mat%numden(mt_iso)*barn
                    !if(zai==621490) print *, 'Sm-149',mat%numden(mt_iso)*barn
                    !!$omp atomic
                    !tot_mass = tot_mass + &
                    !nuclide(inum,nnum,anum)%amu*mat%numden(mt_iso)*mat%vol*m_u*1000.d0
                    !iso = mat%ace_idx(mt_iso)
                    !if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then 
                    !    !$omp atomic
                    !    tot_fmass = tot_fmass + &
                    !    ace(iso)%*mat%numden(mt_iso)*mat%vol*m_u
                    !endif
                endif 
            endif
            !if (zai_idx(iso_idx(mt_iso))/10 == 54135) Xe_pointer(imat) = mt_iso
        end do
        deallocate(mat%ogxs); allocate(mat%ogxs(1:num_iso,1:7)); mat%ogxs(1:num_iso,1:7) = 0.0d0
        ! ======================================================================================
        if (icore==score) then
            write(prt_bumat, *) 'Num isotope', i
            write(prt_bumat, *) 'LOSS TERM', remsum*barn
            write(prt_bumat, *) ''
        endif 
    end do

    ! INITIALIZE FRATIO
    fratio = 0.d0

    if(icore==score) print *, 'Total Flux', tot_flux
    do i = 1,n_materials
        mat => materials(i)
        do mt_iso = 1,mat%n_iso
            if(icore==score) then
            zai = zai_idx(iso_idx(mt_iso))
            if (zai>0) then
                anum = zai/10000; mnum = (zai-anum*10000)/10
                nnum = mnum-anum; inum = zai-anum*10000-mnum*10
                !$omp atomic
                tot_mass = tot_mass + &
                ace(mat%ace_idx(mt_iso))%atn*mat%numden(mt_iso)*mat%vol/N_AVOGADRO*m_n
                iso = mat%ace_idx(mt_iso)
                if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
                    tot_fmass = tot_fmass + &
                    ace(mat%ace_idx(mt_iso))%atn*mat%numden(mt_iso)*mat%vol/N_AVOGADRO*m_n
                endif
            endif
            endif
        enddo
    enddo
    if(icore==score) then
        write(prt_bumat,*) 'Total mass[g]:', tot_mass 
        write(prt_bumat,*) 'Fiss. mass[g]:', tot_fmass
    endif
    istep_burnup = istep_burnup + 1
    deallocate(yield_data) 
     
    end subroutine depletion 
    
    ! ===================================================================
    !         Tally burnup parameters 
    ! ===================================================================
    subroutine tally_burnup (imat, distance, wgt, erg)
        integer, intent(in) :: imat  ! material index (p%material) 
        real(8), intent(in) :: distance, wgt, erg
        integer :: iso, mt_iso
        integer :: i, ierg
        integer :: fy_midx, diff, tmp, Xe_ptr
        real(8) :: g
        real(8),allocatable :: Ep(:)
        integer :: nE, eg

        type (Material_CE), pointer :: mat
        real(8) :: micro_f, micro_d, micro_a, ipfac
        real(8) :: micro_2n, micro_3n, micro_4n, micro_p
        !real(8) :: dep_ogxs(7)
        integer :: anum, mnum, inum, nnum !atomic number, mass number, isomer state, neutron number of nuclide
        integer :: a1, m1
        real(8) :: fluxtmp
        integer :: fssn

        if(E_mode==0) return       !Multigroup  -> return
        if(do_burn==.false.) return     !No burnup calculation -> return
        if(curr_cyc <= n_inact) return 
        if(.not. materials(imat)%depletable) return ! not depletable -> return
        
        
        mat => materials(imat)
        
        !Inline Xenon equilibrium search
        if(do_Xe_search) then
            !Xenon production
            do mt_iso = 1, mat%n_iso
                iso = mat%ace_idx(mt_iso)
                if(.not. allocated(ace(iso)%sigf)) cycle
                call getierg(iso,ierg,erg)
                ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
                micro_d   = ace(iso)%sigd(ierg) + ipfac*(ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
                ! Fissionable Material
                if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
                    micro_f   = ace(iso)%sigf(ierg) + ipfac*(ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
                else
                    micro_f   = 0.d0
                endif
                micro_a = micro_d + micro_f
                
                anum = ace(iso)%zaid/1000
                mnum = (ace(iso)%zaid - anum*1000)
                nnum = mnum - anum
                inum = 0
                if(nuclide(inum,nnum,anum)%fy_idx>0) then
                    fy_midx = nuclide(inum,nnum,anum)%fy_idx
                else
                    !diff = 999999
                    !do i=1,nfssn
                    !    tmp = abs(fssn_zai(i) - zai_idx(nuclide(inum,nnum,anum)%idx))
                    !    if(diff > tmp) then
                    !    diff = tmp
                    !    fy_midx = i
                    !    end if
                    !end do
                    fy_midx = 0
                    if(inum==1 .and. nuclide(0,nnum,anum)%fy_idx>0) &
                        fy_midx = nuclide(0,nnum,anum)%fy_idx
                    do i = 1,nfssn
                        a1 = fssn_zai(i)/10000
                        m1 = (fssn_zai(i)-a1*10000)/10
                        if(m1==mnum) fy_midx = i
                    enddo
                    if(fy_midx==0) fy_midx = nuclide(0,143,92)%fy_idx
                endif
                Xe_prod(imat) = Xe_prod(imat) + wgt*distance*mat%numden(mt_iso)*micro_f*Xe_cum_yield(fy_midx)
                I_prod(imat) = I_prod(imat) + wgt*distance*mat%numden(mt_iso)*micro_f*I_cum_yield(fy_midx)
            end do
    
            !Xenon transmutation 
            iso = Xe135_iso
            call getierg(iso,ierg,erg)
            ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
            micro_a = ace(iso)%sigd(ierg) + ipfac*(ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
            Xe_trans(imat) = Xe_trans(imat) + wgt*distance*micro_a*1.0d-24
        end if
        
        !$omp atomic
        mat%flux = mat%flux + wgt*distance !Volume-integrated flux tally (not normalized)
         do iso = 1,num_iso
         ! TEMPORARY... need to reduce
            call getierg(iso,ierg,erg)
            ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
            micro_d = ace(iso)%sigd(ierg) + ipfac*(ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
            !micro_2n= getxs(N_2N,iso,erg,ierg)
            !micro_3n= getxs(N_3N,iso,erg,ierg)
            !micro_4n= getxs(N_4N,iso,erg,ierg)
            !micro_p = getxs(N_P ,iso,erg,ierg)
            !micro_a = getxs(N_A ,iso,erg,ierg)
            micro_f = 0.d0
            !FISSIONABLE
            if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
                micro_f = ace(iso)%sigf(ierg) + ipfac*(ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
                fssn = ace_fssn(iso)
                if(fssn>0) then
                    Ep = yieldE(fssn,1:4); nE = yieldnE(fssn)
                    if(nE==0) cycle !If all Ep are 0, pass (cannot happen)
                    ! ASSUME E monotonically increases in Ep
                    Ep = Ep(1:nE)
                    if(nE==1) then
                        !$OMP ATOMIC
                        fratio(fssn,1) = fratio(fssn,1) + micro_f*wgt*distance*barn
                    else
                        if(erg<Ep(1)) then !Solely thermal
                            !$OMP ATOMIC
                            fratio(fssn,1) = fratio(fssn,1) + micro_f*wgt*distance*barn
                        elseif(erg>=Ep(nE)) then !Solely fast
                            !%OMP ATOMIC
                            fratio(fssn,nE) = fratio(fssn,nE) + micro_f*wgt*distance*barn
                        else
                            do eg = 1,nE-1 !Already nE >= 2
                                if(erg>=Ep(eg) .and. erg<Ep(eg+1)) then
                                    g = (erg-Ep(eg))/(Ep(eg+1)-Ep(eg))
                                    !$OMP ATOMIC
                                    fratio(fssn,eg) = fratio(fssn,eg) + micro_f*wgt*distance*barn*(1.d0-g)
                                    !$OMP ATOMIC
                                    fratio(fssn,eg+1) = fratio(fssn,eg+1) + micro_f*wgt*distance*barn*g
                                endif
                            enddo
                        endif
                    endif
                endif
            endif
            !dep_ogxs = (/micro_d, micro_2n, micro_3n, micro_4n, micro_p, micro_a, micro_f/)
            !do i = 1,7
            !    if(dep_ogxs(i)==0.d0) cycle
            !    !$OMP ATOMIC
            !    mat%ogxs(iso,i) = mat%ogxs(iso,i) + dep_ogxs(i)*wgt*distance*barn
            !enddo
            !$OMP ATOMIC
            mat%ogxs(iso,1) = mat%ogxs(iso,1) + micro_d*wgt*distance*barn
            !$OMP ATOMIC
            mat%ogxs(iso,2) = mat%ogxs(iso,2) + getxs(N_2N,iso,erg,ierg)*wgt*distance*barn
            !$OMP ATOMIC
            mat%ogxs(iso,3) = mat%ogxs(iso,3) + getxs(N_3N,iso,erg,ierg)*wgt*distance*barn
            !$OMP ATOMIC
            mat%ogxs(iso,4) = mat%ogxs(iso,4) + getxs(N_4N,iso,erg,ierg)*wgt*distance*barn
            !$OMP ATOMIC
            mat%ogxs(iso,5) = mat%ogxs(iso,5) + getxs(N_P,iso,erg,ierg)*wgt*distance*barn
            !$OMP ATOMIC
            mat%ogxs(iso,6) = mat%ogxs(iso,6) + getxs(N_A,iso,erg,ierg)*wgt*distance*barn
            !$OMP ATOMIC
            mat%ogxs(iso,7) = mat%ogxs(iso,7) + micro_f*wgt*distance*barn
         enddo

    end subroutine tally_burnup
    
    
    ! ===================================================================
    !         Inline Xe-135  
    ! ===================================================================
    subroutine inline_xenon
        integer :: imat
        integer :: loc1, loc2
        real(8) :: inv_nactive, avg_Xe, avg_I, temp_I, temp_Xe, rcvbufXe(n_materials)
        type (Material_CE), pointer :: mat
        if (.not. do_Xe_search) return 
       
        call MPI_ALLREDUCE(Xe_prod, rcvbufXe, n_materials, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)    
        call MPI_ALLREDUCE(Xe_trans, rcvbufXe, n_materials, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)    
        call MPI_ALLREDUCE(I_prod, rcvbufXe, n_materials, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)

        inv_nactive = 1.0d0 / dble(n_act)
        avg_Xe = 0.d0; avg_I = 0.d0
        do imat = 1, n_materials
            mat => materials(imat)

            if(mat%fissionable == .false.) cycle
    
            !Obtain volume averaged Xe_tally
            Xe_prod(imat)   = Xe_prod(imat)/materials(imat)%vol   !Xenon production
            Xe_trans(imat) = Xe_trans(imat)/materials(imat)%vol !Xenon transmutation
            I_prod(imat) = I_prod(imat)/materials(imat)%vol
            temp_Xe = Xe_prod(imat)/(Xe_trans(imat) + nuclide(0,81,54)%lambda)*barn
            temp_I = I_prod(imat)/nuclide(0,82,53)%lambda*barn
            !avg_Xe = avg_Xe + temp_Xe; avg_I = avg_I + temp_I
            !Set Xenon equilibrium density
            materials(imat)%numden(Xe135_mt_iso(imat)) = temp_Xe
            materials(imat)%numden(I135_mt_iso(imat)) = temp_I
            if(curr_cyc > n_inact) then
                !Accumulate Xenon number density for active cycles
                Xe_numden(imat) = Xe_numden(imat) + temp_Xe
                I_numden(imat) = I_numden(imat) + temp_I
                !Update Xenon number density for depletion calculation
                if(curr_cyc == n_totcyc) then
                    Xe_numden(imat) = Xe_numden(imat)*inv_nactive
                    materials(imat)%numden(Xe135_mt_iso(imat)) =  Xe_numden(imat)
                    I_numden(imat) = I_numden(imat)*inv_nactive
                    materials(imat)%numden(I135_mt_iso(imat)) = I_numden(imat)
                end if
            end if
        end do
        !if(icore==score)then
        !print *, avg_Xe, avg_I
        !endif
        Xe_prod = 0.d0; Xe_trans = 0.d0; I_prod = 0.d0
        
        !if(icore==score) print *, "inline xenon equilibrium updated"
        
    end subroutine inline_xenon
    
    
    ! ===================================================================
    !         Initialize burnup tallies 
    ! ===================================================================
    subroutine Init_burnup 
        integer :: i, imat
        
        if (istep_burnup == NSTEP_BURNUP) do_burn = .false. 
        if (.not. do_burn) return
        
        if (istep_burnup > 0) do_Xe_search = (Xe_search)! .and. bstep_size/86400.d0 > 5.00001d0)
        cram_init = .false.
        do imat = 1, n_materials
            materials(imat)%flux = 0.0d0 
            !materials(imat)%pwr = 0.0d0 
            materials(imat)%ogxs(:,:) = 0.0d0
            !allocate(materials(imat)%full_numden(nnuc))
        enddo
        ! TESTING for NFY
        fratio = 0.d0

        bstep_size = burn_step(istep_burnup+1) - burn_step(istep_burnup) 
        avg_power = 0.d0; avg_fiss = 0.d0
        fis_thermal = 0.d0; fis_epi = 0.d0; fis_fast = 0.d0
        tot_mass = 0.d0; tot_fmass = 0.d0
        if(icore==score) print *, 'INIT DONE'
        !> Xe equilibrium
        if (.not. do_Xe_search) return 
        Xe_prod(:) = 0.0d0
        Xe_trans(:) = 0.0d0
        Xe_numden(:) = 0.0d0
        I_prod(:) = 0.0d0
        I_numden(:) = 0.0d0
        do imat = 1, n_materials 
            Xe_local: do i = 1, materials(imat)%n_iso 
                if (materials(imat)%ace_idx(i) == Xe135_iso) then
                    Xe135_mt_iso(imat) = i 
                    exit Xe_local
                endif 
            enddo Xe_local
            I_local: do i = 1,materials(imat)%n_iso
                if(materials(imat)%ace_idx(i)==I135_iso) then
                    I135_mt_iso(imat) = i
                    exit I_local
                endif
            enddo I_local
        enddo 
        
    end subroutine Init_burnup 
    
    
    ! ===================================================================
    !     Set material library for burnup and equilibrium Xe135 search 
    ! ===================================================================
    subroutine setmat
        integer :: i, imat, iso, mt_iso
        integer :: anum, mnum, inum, nnum 
        !atomic number, mass number, isomer state, neutron number of nuclide
        
        if (.not. do_burn) return 
        
        cram_init = .false.
        do imat = 1, n_materials
            materials(imat)%flux = 0.0d0 
            !materials(imat)%pwr = 0.0d0 
            !allocate(materials(imat)%ogxs(1:materials(imat)%n_iso, 1:4)) 
            allocate(materials(imat)%ogxs(1:num_iso,1:7))
            materials(imat)%ogxs(:,:) = 0.0d0
            !allocate(materials(imat)%eflux(1:ngrid))
            !materials(imat)%eflux(:) = 0.d0
            
            if (materials(imat)%depletable) then
                allocate(materials(imat)%full_numden(1:nnuc))
                materials(imat)%full_numden = 0.d0     !Initialize number density
            endif
            
            do mt_iso = 1, materials(imat)%n_iso
                iso = materials(imat)%ace_idx(mt_iso)
                
                anum = ace(iso)%zaid/1000
                mnum = ace(iso)%zaid - anum*1000
                nnum = mnum - anum
                if(materials(imat)%depletable==.true. .and. do_burn .and. nuclide(0,nnum,anum)%idx>0) then 
                    materials(imat)%full_numden(nuclide(0,nnum,anum)%idx) = materials(imat)%numden(mt_iso)
                endif 
            end do
            
        enddo 
        fratio = 0.d0 
        istep_burnup = 0 
        
        !> Xe 
        if (.not. Xe_search) return 
        
        !> Find Xe 135 index
        !Xe135_iso = 0 
        !I135_iso = 0
        !Xe_search: do i = 1, num_iso 
        !    if (ace(i)%zaid == 54135) then     
        !        Xe135_iso = i 
        !        exit Xe_search
        !    endif 
        !enddo Xe_search
        !I_search: do i = 1,num_iso
        !    if (ace(i)%zaid == 53135) then
        !        I135_iso = i
        !        exit I_search
        !    endif
        !enddo I_search

        allocate(Xe_prod(1:n_materials))
        allocate(Xe_trans(1:n_materials))
        allocate(Xe_numden(1:n_materials))
        allocate(Xe135_mt_iso(1:n_materials))
        
        allocate(I_prod(1:n_materials))
        allocate(I_numden(1:n_materials))
        allocate(I135_mt_iso(1:n_materials))
        
        Xe135_iso = find_ACE_iso_idx_zaid(541350) 
        I135_iso = find_ACE_iso_idx_zaid(531350)

        if (Xe135_iso == 0 .and. icore==score ) then 
            print *, "setmat() - WARNING :: NO Xe-135 LIBRARY IN THE INVENTORY"
        endif 

        if (I135_iso == 0 .and. icore==score) then
            print *, "NO I-135 library"
        endif
        
        do imat = 1, n_materials            
            do mt_iso = 1, materials(imat)%n_iso
                iso = materials(imat)%ace_idx(mt_iso)
                if(Xe135_iso == iso) Xe135_mt_iso(imat) = mt_iso
                if(I135_iso == iso) I135_mt_iso(imat) = mt_iso
            end do
        enddo 
        
        
        
    end subroutine setmat
    
    
    
    ! ===================================================================
    !         MPI reduce burnup tallies 
    ! ===================================================================
    subroutine MPI_reduce_burnup 
        integer :: iso, imat, i 
        real(8) :: rcvbuf, rcvbufarr(7)
        real(8), allocatable :: sndbufarrlong(:), rcvbufarrlong(:)
        real(8) :: val
        
        if (.not. do_burn) return 
        
        
        do imat = 1, n_materials
            if (.not. materials(imat)%depletable) cycle 
            
            call MPI_ALLREDUCE(materials(imat)%flux, rcvbuf, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, core, ierr)
            materials(imat)%flux = rcvbuf / (dble(n_act) * materials(imat)%vol)
            
            !call MPI_ALLREDUCE(materials(imat)%pwr, rcvbuf, 1, MPI_DOUBLE_PRECISION, &
            !                MPI_SUM, core, ierr)
            !materials(imat)%pwr = rcvbuf / dble(n_act)
            
            !allocate(rcvbufarrlong(1:materials(imat)%n_iso*4)) 
            !allocate(sndbufarrlong(1:materials(imat)%n_iso*4)) 
            allocate(rcvbufarrlong(1:num_iso*7))
            allocate(sndbufarrlong(1:num_iso*7))

            val = dble(n_act) * materials(imat)%flux * materials(imat)%vol
            !do iso = 1, materials(imat)%n_iso
            do iso = 1,num_iso
                !call MPI_ALLREDUCE(materials(imat)%ogxs(iso,:), rcvbufarr, 4, MPI_DOUBLE_PRECISION, &
                !                MPI_SUM, core, ierr)
                
                sndbufarrlong((iso-1)*7+1:(iso-1)*7+7) = materials(imat)%ogxs(iso,:)
                
                !materials(imat)%ogxs(iso,:) = rcvbufarr(:) / val
            enddo 
            !call MPI_ALLREDUCE(sndbufarrlong, rcvbufarrlong, materials(imat)%n_iso*4, &
            !                    MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)            
            call MPI_ALLREDUCE(sndbufarrlong, rcvbufarrlong, num_iso*7, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)                    
            rcvbufarrlong(:) = rcvbufarrlong(:) / val
            
            !do iso = 1, materials(imat)%n_iso
            do iso = 1,num_iso
                materials(imat)%ogxs(iso,:) = rcvbufarrlong((iso-1)*7+1:(iso-1)*7+7) 
            enddo
            deallocate(rcvbufarrlong)
            deallocate(sndbufarrlong) 

            ! EFLUX
            !allocate(rcvbufarrlong(1:ngrid)); allocate(sndbufarrlong(1:ngrid))

        enddo

    end subroutine MPI_reduce_burnup
    
    
    ! ===================================================================
    !         Chebyshev rational approximation method (CRAM) 
    ! ===================================================================
    subroutine cram(burnupMat, nuc0, nuc1) 
        use mkl_pardiso
        implicit none 
        include 'mkl_spblas.fi'   !! or mkl.fi (not mandatory but recommended)
        
        real(8), intent(in) :: burnupMat(:,:)
        real(8), intent(in) :: nuc0(:)
        real(8), intent(inout) :: nuc1(:)
    
        !Cram theta and alpha
        complex(8), parameter :: cram14_theta(1:7) = &
        & (/(-8.8977731864688888199d0,   1.6630982619902085304d1), &
        &   (-3.7032750494234480603d0,   1.3656371871483268171d1), &
        &   (-0.2087586382501301251d0,   1.0991260561901260913d1), &
        &   (3.9933697105785685194d0,  6.0048316422350373178d0), &
        &   (5.0893450605806245066d0,  3.5888240290270065102d0), &
        &   (5.6231425727459771248d0,  1.1940690463439669766d0), &
        &   (2.2697838292311127097d0,  8.4617379730402214019d0)/)
        complex(8), parameter :: cram14_alpha(0:7) = &
        & (/( 1.8321743782540412751d-14, 0.d0), &
        &   (-7.1542880635890672853d-5,  1.4361043349541300111d-4), &
        &   ( 9.4390253107361688779d-3, -1.7184791958483017511d-2), &
        &   (-3.7636003878226968717d-1,  3.3518347029450104214d-1), &
        &   (-2.3498232091082701191d1,  -5.8083591297142074004d0),  &
        &   ( 4.6933274488831293047d1,   4.5643649768827760791d1),  &
        &   (-2.7875161940145646468d1,  -1.0214733999056451434d2),  &
        &   ( 4.8071120988325088907d0,  -1.3209793837428723881d0)/)
        complex(8), parameter :: cram16_theta(1:8) = &
        & (/(-1.0843917078696988026d1,   1.9277446167181652284d1),  &
        &   (-5.2649713434426468895d0,   1.6220221473167927305d1),  &
        &   ( 5.9481522689511774808d0,   3.5874573620183222829d0),  &
        &   ( 3.5091036084149180974d0,   8.4361989858843750826d0),  &
        &   ( 6.4161776990994341923d0,   1.1941223933701386874d0),  &
        &   ( 1.4193758971856659786d0,   1.0925363484496722585d1),  &
        &   ( 4.9931747377179963991d0,   5.9968817136039422260d0),  &
        &   (-1.4139284624888862114d0,   1.3497725698892745389d1)/)
        complex(8), parameter :: cram16_alpha(0:8) = &
        & (/( 2.1248537104952237488d-16, 0.d0), &
        &   (-5.0901521865224915650d-7, -2.4220017652852287970d-5), &
        &   ( 2.1151742182466030907d-4,  4.3892969647380673918d-3), &
        &   ( 1.1339775178483930527d2,   1.0194721704215856450d2),  &
        &   ( 1.5059585270023467528d1,  -5.7514052776421819979d0),  &
        &   (-6.4500878025539646595d1,  -2.2459440762652096056d2),  &
        &   (-1.4793007113557999718d0,   1.7686588323782937906d0),  &
        &   (-6.2518392463207918892d1,  -1.1190391094283228480d1),  &
        &   ( 4.1023136835410021273d-2, -1.5743466173455468191d-1)/)
    
        integer :: i, j, k
        complex(8) :: A(1:nnuc,1:nnuc)
        complex(8) :: b(1:nnuc)
        complex(8) :: sol(1:nnuc)
        integer :: ii,jj,kk
    
        !MKL indicator
        integer :: info     !Indicator only for restoring the matrix A from CSR format
    
        !PARDISO control parameter
        integer :: phase
        integer :: idum(1)
        complex(8) :: zdum(1)
        integer :: error, error1
        TYPE(MKL_PARDISO_HANDLE) :: burn_pt(1:64)
		integer :: iparm(1:64) !solver control
		integer :: mtype       !matrix type 
		integer :: msglvl      !message level
    
		!if(istep_burnup > 0) then 
		!	cram_init = .true. 
		!else 
		!	cram_init = .false.
		!endif	
	
        !Set control parameters
        !if(cram_init == .false.) then
        job(1) = 0    !convert rectangular matrix A to CSR format
        job(2:3) = 1  !use 1-based indices
        job(4) = 2    !use whole matrix A as input
        job(5) = nnuc*nnuc !maximum allowed number of nonzere elements
        job(6) = 1    !generate Asparse, ia, and ja, as output
        
        iparm = 0
        iparm(1) = 1 !Not default values
        iparm(2) = 2 !Nested dissection algorithm from METIS package
        iparm(3) = 0 !Reserved (set to zero)
        iparm(4) = 0 !Factorization is always computed as required by phase
        iparm(5) = 0 !Perm is ignored
        iparm(6) = 0 !Write solution on x
        iparm(8) = 10000 !Max. numbers of iterative refinement steps
        iparm(10) = 13 !Pivoting perturbation (10**(-iparm(10)))
        iparm(11) = 1 !Scaling vectors (enabling scaling default for unsymmetric matrix)
        iparm(13) = 1 !Improved accuracy using nonsymmetric weighted matching
        
        mtype = 13   !(13 = complex and nonsymmetric)
        msglvl = 0   !Statistical information is printed to screen
        do i=1, 64
            burn_pt(i)%dummy = 0
        end do
		
		
		
        !end if
    
        select case(cram_order)
        case(14)
        nuc1 = nuc0*cram14_alpha(0)
        do i=1, 7
            A = burnupMat
            b = cram14_alpha(i)*nuc0
            do j=1, nnuc 
                A(j,j) = A(j,j) - cram14_theta(i)
            end do
    
            !Full matrix to CSR format
            if(cram_init==.false.) then
                call mkl_zdnscsr(job, nnuc, nnuc, A, nnuc, Acsr, jAcsr, iAcsr, info)
                cram_init = .true.
            else
                kk=0
                do ii=1, nnuc
                    do jj=1, iAcsr(ii+1)-iAcsr(ii)
                        kk = kk + 1
                        Acsr(kk) = A(ii,jAcsr(kk))
                    end do
                end do
            end if
        
            !Symbolic factorization (allocates all memory required for factorization)
            if(i==1) then
                phase = 11
                call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, zdum, zdum, error)
            end if
    
            !Factorization
            phase = 22
            call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, zdum, zdum, error)
    
            !Back substitution and iterative refinement
            phase = 33
            call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, b, sol, error)
            nuc1 = nuc1 + 2.d0*dble(sol)
        end do
    
        !Termination and release memory used for factorization
        phase = -1 
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, zdum, zdum, error)
    
        case(16)
        nuc1 = nuc0*cram16_alpha(0)
        do i=1, 8
            A = burnupMat
            b = cram16_alpha(i)*nuc0
            do j=1, nnuc 
            A(j,j) = A(j,j) - cram16_theta(i)
            end do
    
            !Full matrix to CSR format
            if(cram_init==.false.) then
            call mkl_zdnscsr(job, nnuc, nnuc, A, nnuc, Acsr, jAcsr, iAcsr, info)
            cram_init = .true.
            else
            kk=0
            do ii=1, nnuc
                do jj=1, iAcsr(ii+1)-iAcsr(ii)
                kk = kk + 1
                Acsr(kk) = A(ii,jAcsr(kk))
                end do
            end do
            end if
        
            !Symbolic factorization (allocates all memory required for factorization)
            if(i==1) then
            phase = 11
            call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, zdum, zdum, error)
            end if
    
            !Factorization
            phase = 22
            call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, zdum, zdum, error)
    
            !Back substitution and iterative refinement
            phase = 33
            call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, b, sol, error)
            nuc1 = nuc1 + 2.d0*dble(sol)
        end do
    
        !Termination and release memory used for factorization
        phase = -1 
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, idum, 1, iparm, msglvl, zdum, zdum, error)
    
    
        case default 
        if (icore==score) print *, "CRAM() :: WRONG CRAM ORDER", cram_order
        stop 
        end select
    
    
    end subroutine
    
!subroutine r8mat_expm1 ( n, a, e )
!
!  implicit none
!
!  integer ( kind = 4 ) n
!
!  real ( kind = 8 ) a(n,n)
!  real ( kind = 8 ) a2(n,n)
!  real ( kind = 8 ) a_norm
!  real ( kind = 8 ) c
!  real ( kind = 8 ) d(n,n)
!  real ( kind = 8 ) e(n,n)
!  integer ( kind = 4 ) ee
!  integer ( kind = 4 ) k
!  logical p
!  integer ( kind = 4 ) , parameter :: q = 6
!  a_norm = r8mat_norm_li ( n, n, a2 )
!
!  ee = int ( r8_log_2 ( a_norm ) ) + 1
!
!  s = max ( 0, ee + 1 )
!
!  a2(1:n,1:n) = a2(1:n,1:n) / 2.0D+00**s
!
!  x(1:n,1:n) = a2(1:n,1:n)
!
!  c = 0.5D+00
!
!  call r8mat_identity ( n, e )
!  e(1:n,1:n) = e(1:n,1:n) + c * a2(1:n,1:n)
!
!  call r8mat_identity ( n, d )
!  d(1:n,1:n) = d(1:n,1:n) - c * a2(1:n,1:n)
!
!  p = .true.
!
!  do k = 2, q
!
!    c = c * real ( q - k + 1, kind = 8 ) &
!      / real ( k * ( 2 * q - k + 1 ), kind = 8 )
!
!    x(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )
!
!    e(1:n,1:n) = e(1:n,1:n) + c * x(1:n,1:n)
!
!    if ( p ) then
!      d(1:n,1:n) = d(1:n,1:n) + c * x(1:n,1:n)
!    else
!      d(1:n,1:n) = d(1:n,1:n) - c * x(1:n,1:n)
!    end if
!
!    p = .not. p
!
!  end do
!!
!!  E -> inverse(D) * E
!!
!  call r8mat_minvm ( n, n, d, e, e )
!  
!  
!!  E -> E^(2*S)
!  do k = 1, s
!    if (icore==(ncore-1)) print *, 'check k'  , k, '/',s
!    e(1:n,1:n) = matmul ( e(1:n,1:n), e(1:n,1:n) )
!  end do
!  return
!end subroutine 
!
!subroutine r8mat_expm2 ( n, a, e )
!
!!*****************************************************************************80
!!
!!! R8MAT_EXPM2 uses the Taylor series for the matrix exponential.
!!
!!  Discussion:
!!
!!    Formally,
!!
!!      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...
!!
!!    This function sums the series until a tolerance is satisfied.
!!
!!  Licensing:
!!
!!    This code is distributed under the GNU LGPL license.
!!
!!  Modified:
!!
!!    26 November 2011
!!
!!  Author:
!!
!!    Cleve Moler, Charles Van Loan
!!
!!  Reference:
!!
!!    Cleve Moler, Charles VanLoan,
!!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!!    Twenty-Five Years Later,
!!    SIAM Review,
!!    Volume 45, Number 1, March 2003, pages 3-49.
!!
!!  Parameters:
!!
!!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!!
!!    Input, real ( kind = 8 ) A(N,N), the matrix.
!!
!!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!!
!  implicit none
!
!  integer ( kind = 4 ) n
!
!  real ( kind = 8 ) a(n,n)
!  real ( kind = 8 ) e(n,n)
!  real ( kind = 8 ) f(n,n)
!  real ( kind = 8 ) g(n,n)
!  integer ( kind = 4 ) k
!  logical r8mat_is_insignificant
!
!  e(1:n,1:n) = 0.0D+00
!
!  call r8mat_identity ( n, f )
!
!  k = 1
!
!  do
!
!    if ( r8mat_is_insignificant ( n, n, e, f ) ) then
!      exit
!    end if
!
!    e(1:n,1:n) = e(1:n,1:n) + f(1:n,1:n)
!
!    f(1:n,1:n) = matmul ( a(1:n,1:n), f(1:n,1:n) ) / real ( k, kind = 8 )
!    k = k + 1
!
!  end do
!
!  return
!end subroutine

end module
