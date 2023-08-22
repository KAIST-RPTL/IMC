
include 'mkl_pardiso.f90'

module depletion_module 
    use variables
    use constants
    use ace_header, only: ace, find_ACE_iso_idx_zaid, num_iso, Emin, Emax, nueg, ueggrid
    use ace_module, only: set_ace_iso
    use ace_xs,     only: getxs, getierg, getiueg, setueg
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
      integer :: conn  = -1
    end type
    type(nuclide_data) :: nuclide(0:2,0:170,0:111) !nuclide data for Isomeric State (0=ground, 1=excited), neutron number, atomic number 
    ! NFY TESTING
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
    integer, allocatable :: fiss_path(:,:)

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
    logical :: bumat_print = .false.

    integer :: RXMT(7)

    ! EFLUX
    integer :: depopt = 0
    integer :: ngrid
    real(8) :: gdelta
    integer, parameter :: numrx = 7
    real(8), allocatable :: XS(:,:)


    ! MPI REGION
    integer, allocatable :: mpigeom(:,:), tmpgeom(:)
    integer :: ngeom, totgeom

    integer :: NFYtype = 3

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
        integer :: i, j, k, ii
        
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

        real(8) :: erg
        integer :: ierg

        integer, allocatable :: decay(:), daugh(:)
        integer :: anum1, nnum1, mnum1, inum1, rnum
        integer :: tgt,   idx,   idx1,  conval
        logical :: ngbranch

        character(4) :: tail
        integer :: liblen

        integer :: nn, pn, dn, tn, an, a3n
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
        if (nuclide(inum,nnum,anum)%data_exist .and. nuclide(inum,nnum,anum)%iu>0) go to 15
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
                if(.not. nuclide(prod(1),prod(2),prod(3))%data_exist) then
                    nuclide(prod(1), prod(2), prod(3)) % data_exist = .true.
                    nuclide(prod(1), prod(2), prod(3)) % iu = 0
                endif
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
!        !Read one-group transmutation cross section library
!        open(rd_1gcx, file=trim(dep_lib)//"1gcx_library",action="read")
!        20 continue
!        read(rd_1gcx,'(A80)',end=200) line0 !Title
!        do
!          read(rd_1gcx,'(A80)',end=200) line1
!          
!          read(line1(1:4),'(i)') nlb; if(nlb==-1) go to 20 !library index (1=activation products, 2=actinides, 3=fission products)
!          read(line1(7:12),'(i)') nuclid
!            anum = nuclid/10000
!            mnum = (nuclid - anum*10000)/10
!            nnum = mnum - anum
!            inum = nuclid - anum*10000 - mnum*10
!          read(line1(14:22),'(f)') nuclide(inum,nnum,anum)%sng  !(n,g) leading to ground state
!          nuclide(inum,nnum,anum)%sng=nuclide(inum,nnum,anum)%sng*1.d-24
!          read(line1(24:32),'(f)') nuclide(inum,nnum,anum)%sn2n !(n,2n) leading to ground state
!          nuclide(inum,nnum,anum)%sn2n=nuclide(inum,nnum,anum)%sn2n*1.d-24
!          if(nlb/=3) then !fission products or activation products
!            read(line1(34:42),'(f)') nuclide(inum,nnum,anum)%sna !(n,alpha) leading to ground state
!            nuclide(inum,nnum,anum)%sna=nuclide(inum,nnum,anum)%sna*1.d-24
!            read(line1(44:52),'(f)') nuclide(inum,nnum,anum)%snp !(n,proton) leading to ground state
!            nuclide(inum,nnum,anum)%snp=nuclide(inum,nnum,anum)%snp*1.d-24
!          else
!            read(line1(34:42),'(f)') nuclide(inum,nnum,anum)%sn3n!(n,3n) leading to ground state
!            nuclide(inum,nnum,anum)%sn3n=nuclide(inum,nnum,anum)%sn3n*1.d-24
!            read(line1(44:52),'(f)') nuclide(inum,nnum,anum)%snf !(n,f)
!            nuclide(inum,nnum,anum)%snf=nuclide(inum,nnum,anum)%snf*1.d-24
!          end if 
!          read(line1(54:62),'(f)') nuclide(inum,nnum,anum)%sngx  !(n,g) leading to excited state
!          nuclide(inum,nnum,anum)%sngx=nuclide(inum,nnum,anum)%sngx*1.d-24
!          read(line1(64:72),'(f)') nuclide(inum,nnum,anum)%sn2nx !(n,2n) leading to excited state
!          nuclide(inum,nnum,anum)%sn2nx=nuclide(inum,nnum,anum)%sn2nx*1.d-24
!          read(line1(74:79),'(f)') yyn               !yyn > 0 : fission yield card follows, yyn < 0 : no fission yield card
!          if(yyn>0.d0) read(rd_1gcx,'(A80)',end=200) line2
!        end do
!        200 close(rd_1gcx)


        !==============================================================================
        open(rd_yield, file = trim(dep_lib)//'sss_endfb71.nfy', action = 'read')
        fid = 0; ifp = 0
        allocate(tmp_yield(1:2500,1:4,1:100)); tmp_yield = 0.d0
        allocate(ify_yield(1:2500,1:4,1:100)); ify_yield = 0.d0
        allocate(cfy_yield(1:2500,1:4,1:100)); cfy_yield = 0.d0
        allocate(fiss_path(1:2500,1:100));     fiss_path = 0
        allocate(fp_zai(1:2500)) ! ZAI of fp
        allocate(fssn_zai(1:100)) ! ZAI of fssn
        allocate(ace_fssn(1:2000)); ace_fssn = 0.d0
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
!        if(find_ACE_iso_idx_zaid(nuclid*10)>0) then
!            ace_fssn(find_ACE_iso_idx_zaid(nuclid*10+inum)) = fid
!        endif
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
                            else
                                tmpfp = nuclide(0, nnum, anum) % fp
                            endif
                        endif
                        flag = mod(fp*4-1,6) ! Flag for yield
                        if (flag<3) read(rd_yield,'(A80)',end=300) line1
                        read(line1(11*flag-10:11*flag),'(f)') ify

                        ! Alternate for Meta...

                        if(ify>0.d0) then
                            if(tmpfp>0) then
                                ify_yield(tmpfp,eg,fid) = &
                                ify_yield(tmpfp,eg,fid) + ify
                            endif
                        endif
                        !if(icore==score .and. anum==33 .and. mnum==74 .and. fssn_zai(fid)==922350) print *, 'NFY', nuclid, inum, ify, eg
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
                        tmpfp = nuclide(inum, nnum, anum) %fp
                        if(tmpfp>0) then
                            if(cfy_yield(tmpfp,eg,fid)>fpcut .or. ify_yield(tmpfp,eg,fid)>fpcut) then
                                tmp_yield(tmpfp,eg,fid) = tmp_yield(tmpfp,eg,fid) + ify_yield(tmpfp,eg,fid)
                                if(fiss_path(tmpfp,fid)==0) fiss_path(tmpfp,fid) = 1
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
                read(line0(23:33),'(i)') n_y_eg
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

                11 continue

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
                enddo
400             close(isom_brn)
                do brn = 1, num_brn
                !if(icore==score) print *, ZAIMT_ism(brn), gnd_frac(brn)
                enddo
                !==============================================================================
                ! CROP NFY DATA
                tmp_yield = tmp_yield(1:ifp,1:4,1:fid)
                fiss_path = fiss_path(1:ifp,1:fid)
                deallocate(ify_yield); deallocate(cfy_yield)
                fp_zai = fp_zai(1:ifp)
                fssn_zai = fssn_zai(1:fid)
                yieldE = yieldE(1:fid,1:4)
                yieldnE = yieldnE(1:fid)
                nfp = ifp
                nfssn = fid
                !allocate(fratio(nfssn,4))
                do i = 1,n_materials
                    allocate(materials(i)%fratio(nfssn,4))
                enddo
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
                RXMT = (/N_GAMMA, N_2N, N_3N, N_4N, N_P, N_A, N_FISSION/)
                allocate(decay(10000)); decay = 0; idx = 0

                do i = 1,num_iso
                    zai = ace(i)%zaid
                    liblen = len_trim(ace(i)%xslib)
                    tail= ace(i)%xslib(liblen-3:liblen)
                    anum= zai/1000
                    mnum= zai-anum*1000
                    inum= 0
                    if(mnum>300) then
                        mnum = mnum-200
                        if(anum>88) mnum = mnum + 100
                        inum = 1
                    elseif(mnum==0) then
                        cycle
                    endif
                    nnum= mnum - anum
                    if(nnum<0) print *, i, zai, tail
                    nuclide(inum,nnum,anum)%conn = 0

                    if(nuclide(0,2,2)%data_exist .and. nuclide(0,2,2)%conn<0) then ! ALPHA
                        nuclide(0,2,2)%conn = 0
                        idx = idx + 1
                        decay(idx) = 20040
                        call ADDACE(decay(idx),tail)
                    endif

                    if(nuclide(0,0,1)%data_exist .and. nuclide(0,0,1)%conn<0) then ! PROTON
                        nuclide(0,0,1)%conn = 0
                        idx = idx + 1
                        decay(idx) = 10010
                        call ADDACE(decay(idx),tail)
                    endif

                    if(nuclide(0,1,0)%data_exist .and. nuclide(0,1,0)%conn<0) then ! NEUTRON
                        nuclide(0,1,0)%conn = 0
                        idx = idx + 1
                        decay(idx) = 00010
                        call ADDACE(decay(idx),tail)
                    endif

                    if(nuclide(0,1,1)%data_exist .and. nuclide(0,1,0)%conn<0) then ! D
                        nuclide(0,1,0)%conn = 0
                        idx = idx + 1
                        decay(idx) = 00010
                        call ADDACE(decay(idx),tail)
                    endif

                    if(nuclide(0,2,1)%data_exist .and. nuclide(0,1,0)%conn<0) then ! T
                        nuclide(0,1,0)%conn = 0
                        idx = idx + 1
                        decay(idx) = 00010
                        call ADDACE(decay(idx),tail)
                    endif

                    if(nuclide(0,1,2)%data_exist .and. nuclide(0,1,0)%conn<0) then ! A3
                        nuclide(0,1,2)%conn = 0
                        idx = idx + 1
                        decay(idx) = 00010
                        call ADDACE(decay(idx),tail)
                    endif

                    do j = 1,ace(i)%NXS(4) ! FOR ALL AVAILABLE RX
!                        do k = 1,7
!                            if(ace(i)%MT(j)==RXMT(k)) then
!                                mt = RXMT(k)
                                mt = ace(i) % MT(j)
                                tgt= anum*1000+mnum
                                select case(mt)
                                case(N_GAMMA)
                                    ngbranch = .false.
                                    do ii = 1,num_brn
                                        if(tgt*10==ZAIMT_ism(ii)) then
                                            ngbranch = .true.
                                            exit
                                        endif
                                    enddo
                                    !if(icore==score) print *, 'GAMMA', mt, anum*1001 + nnum, ngbranch
                                    if(ngbranch .and. nuclide(1,nnum+1,anum)%data_exist .and. nuclide(1,nnum+1,anum)%conn < 0) then
                                        nuclide(1,nnum+1,anum)%conn = 1
                                        idx = idx + 1
                                        decay(idx) = tgt*10+11
                                        call ADDACE(decay(idx),tail)
                                        !if(icore==score) write(*,*) 'ADDED from ', ace(i)%xslib, decay(idx)
                                    endif
                                    if(nuclide(0,nnum+1,anum)%data_exist .and. nuclide(0,nnum+1,anum)%conn<0) then
                                        nuclide(0,nnum+1,anum)%conn = 1
                                        idx = idx + 1
                                        decay(idx) = tgt*10+10
                                        call ADDACE(decay(idx),tail)
                                        !if(icore==score) write(*,*) 'ADDED from ', ace(i)%xslib, ace(i)%MT(j)
                                    endif
                                case(N_FISSION)

                                case default
                                    ! IF direct, skip rest
                                    if(.not. do_rx_tally .or. ANY(RXMT==mt)) then
                                        call mtrxread(mt, nn, pn, dn, tn, an, a3n)  
                                        anum1 = anum - pn - dn - tn - 2*an - 2*a3n
                                        nnum1 = nnum + 1 - nn - dn - 2*tn - 2*an - a3n
                                        ! if(.not. (anum==anum1 .and. nnum==nnum1)) then ! If nuclide changes...
                                        if(nuclide(0, nnum1, anum1) % data_exist .and. nuclide(0, nnum1, anum1) % conn < 0) then
                                            !if(icore==score) print *, &
                                            !    'DAUGH', idx1, mt, anum*1001+ nnum, anum1 *1001 + nnum1
                                            nuclide(0, nnum1, anum1) % conn = 1
                                            idx = idx + 1
                                            decay(idx) = (nnum1+anum1)*10 + anum * 10000
                                            call ADDACE(decay(idx), tail)
                                        endif
                                    endif

                                end select
                            !endif
                        !enddo
                    enddo
                    if(nuclide(inum,nnum,anum)%lambda>0.d0)then !If Decays
                        rnum = nuclide(inum,nnum,anum)%react_num
                        do j = 1,rnum
                            inum1 = nuclide(inum,nnum,anum)%daughter(j,1)
                            nnum1 = nuclide(inum,nnum,anum)%daughter(j,2)
                            anum1 = nuclide(inum,nnum,anum)%daughter(j,3)
                            if(nuclide(inum1,nnum1,anum1)%data_exist .and. nuclide(inum1,nnum1,anum1)%conn<0) then
                                nuclide(inum1,nnum1,anum1)%conn = 99
                                idx = idx + 1
                                decay(idx) = anum1*10010+nnum1*10+inum1
                                call ADDACE(decay(idx),tail)
                                !if(icore==score) write(*,*) 'ADDED from ', ace(i)%xslib, 'decay' 

                            endif
                        enddo
                    endif
                enddo
                do j = 1,nfp
                    anum1 = fp_zai(j)/10000
                    mnum1 = (fp_zai(j)-anum1*10000)/10
                    nnum1 = mnum1 - anum1
                    inum1 = fp_zai(j)-anum1*10000-mnum1*10
                    if(nuclide(inum1,nnum1,anum1)%data_exist .and. nuclide(inum1,nnum1,anum1)%conn<0 .and. nuclide(inum1,nnum1,anum1)%fiss) then
                       nuclide(inum1,nnum1,anum1)%conn=1
                       idx = idx + 1
                       decay(idx) = fp_zai(j)
                       call ADDACE(decay(idx),tail)
                    endif
                enddo
                allocate(daugh(10000)); daugh = 0; conval= 1

                do while(idx>0)
                    !if(conval>5) exit
                    conval = conval + 1
                    idx1 = 0;
                    do i = 1,idx
                        anum = decay(i)/10000
                        mnum = (decay(i)-anum*10000)/10
                        nnum = mnum-anum
                        inum = decay(i)-anum*10000-mnum*10
                        if(nuclide(inum,nnum,anum)%lambda==0.d0) goto 83
                        rnum = nuclide(inum,nnum,anum)%react_num
                        do j = 1,rnum
                           inum1 = nuclide(inum,nnum,anum)%daughter(j,1)
                           nnum1 = nuclide(inum,nnum,anum)%daughter(j,2)
                           anum1 = nuclide(inum,nnum,anum)%daughter(j,3)
                            if(nuclide(inum1,nnum1,anum1)%data_exist .and. nuclide(inum1,nnum1,anum1)%conn<0) then
                               nuclide(inum1,nnum1,anum1)%conn = conval
                               idx1 = idx1 + 1
                               daugh(idx1) = anum1*10010+nnum1*10+inum1
                               call ADDACE(daugh(idx1),tail)
                               !if(icore==score) write(*,*) 'ADDED from ', decay(i), 'Decay'
                            endif
                        enddo

                        83 continue

                        ! Regarding their reactions
                        iso = find_ACE_iso_idx_zaid(decay(i)) 
                        if(iso==0) cycle
                        do j = 1,ace(iso)%NXS(4)
                            !do k = 1,7
                                !if(ace(iso)%MT(j)==RXMT(k)) then
                                !    mt = RXMT(k)
                                    mt = ace(iso) % MT(j)
                                    tgt = decay(i)/10
                                    select case(mt)
                                    case(N_GAMMA)
                                        ngbranch = .false.
                                        do ii = 1,num_brn
                                            if(tgt*10==ZAIMT_ism(ii)) then
                                                ngbranch = .true.
                                                exit
                                            endif
                                        enddo
                                        !if(icore==score) print *, 'GAMMA', mt, anum*1001 + nnum, ngbranch
                                        if(ngbranch .and. nuclide(1,nnum+1,anum)%data_exist .and. nuclide(1,nnum+1,anum)%conn < 0) then
                                            nuclide(1,nnum+1,anum)%conn = conval
                                            idx1 = idx1 + 1
                                            daugh(idx1) = tgt*10+11
                                            call ADDACE(daugh(idx1),tail)
                                            !if(icore==score) write(*,*) 'ADDED from ', ace(iso)%xslib, daugh(idx1)
                                        endif
                                        if(nuclide(0,nnum+1,anum)%data_exist .and. nuclide(0,nnum+1,anum)%conn<0) then
                                            nuclide(0,nnum+1,anum)%conn = conval
                                            idx1 = idx1 + 1
                                            daugh(idx1) = tgt*10+10
                                            call ADDACE(daugh(idx1),tail)
                                            !if(icore==score) write(*,*) 'ADDED from ', ace(iso)%xslib, ace(iso)%MT(j)
                                        endif
                                    case(N_FISSION)
                                    case default
                                        ! IF direct, skip rest
                                        if(.not. do_rx_tally .or. ANY(RXMT==mt)) then
                                            call mtrxread(mt, nn, pn, dn, tn, an, a3n)  
                                            anum1 = anum - pn - dn - tn - 2*an - 2*a3n
                                            nnum1 = nnum + 1 - nn - dn - 2*tn - 2*an - a3n
                                            if(.not. (anum==anum1 .and. nnum==nnum1)) then ! If nuclide changes...
                                                if(nuclide(0, nnum1, anum1) % data_exist .and. nuclide(0, nnum1, anum1) % conn < 0) then
                                                    !if(icore==score) print *, &
                                                    !    'DAUGH', idx1, mt, anum*1001+ nnum, anum1 *1001 + nnum1
                                                    nuclide(0, nnum1, anum1) % conn = conval
                                                    idx1 = idx1 + 1
                                                    daugh(idx1) = (nnum1+anum1)*10 + anum * 10000
                                                    call ADDACE(daugh(idx1), tail)
                                                endif
                                            endif
                                        endif
                                    end select
                                !endif
                            !enddo
                        enddo
!                        if(nuclide(inum,nnum,anum)%fy_idx>0) then
!                        do j = 1,nfp
!                            anum1 = fp_zai(j)/10000
!                            mnum1 = (fp_zai(j)-anum1*10000)/10
!                            nnum1 = mnum1 - anum1
!                            inum1 = fp_zai(j)-anum1*10000-mnum1*10
!                            if(nuclide(inum1,nnum1,anum1)%data_exist .and. nuclide(inum1,nnum1,anum1)%conn<0 .and. nuclide(inum1,nnum1,anum1)%fiss &
!                                .and. fiss_path(j,nuclide(inum,nnum,anum)%fy_idx)==1) then
!                               nuclide(inum,nnum,anum)%conn=1
!                               idx = idx + 1
!                               decay(idx) = fp_zai(j)
!                               call ADDACE(decay(idx),tail)
!                            endif
!                        enddo
!                        endif
                    enddo
                    decay(1:10000) = daugh(1:10000)
                    daugh = 0
                    idx = idx1
                enddo
                
                nnuc = 0
                do anum = 1,111
                do nnum = 0,170
                do inum = 0,2
                    if(nuclide(inum,nnum,anum)%conn>=0) then
                        nnuc = nnuc + 1
                        nuclide(inum,nnum,anum)%idx=nnuc
                        zai_idx(nnuc) = anum*10010+nnum*10+inum
                        !if(icore==score) print *, nnuc, zai_idx(nnuc)
                    endif
                enddo
                enddo
                enddo

                do i = 1,nfssn
                    if(find_ACE_iso_idx_zaid(fssn_zai(i))>0) &
                        ace_fssn(find_ACE_iso_idx_zaid(fssn_zai(i))) = i
                enddo
                ace_fssn = ace_fssn(1:num_iso)

                maxnnz = int(0.1 * nnuc**2) ! Assuming less than 10% sparsity  


                allocate(Acsr(1:maxnnz))
                allocate(jAcsr(1:maxnnz))
                allocate(iAcsr(1:nnuc+1))
                if(icore==score) print *, "Number of isotopes:", nnuc-1 ! Checking number of isotopes
                if(icore==score) print *, "ACE FORMAT ISOTOPES:", num_iso 

                !do i = 1, num_iso
                !    if(abs(ace(i)%TY(rx))==1 .and. .not. (mt==


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
                        write(prt_bumat,*) 'mat: ', materials(i)%mat_name
                        write(prt_bumat,*) ' '
                        do mt_iso = 1,materials(i)%n_iso
                            write(prt_bumat,'(a15,e14.5)')&
                            ace(materials(i)%ace_idx(mt_iso))%xslib,&
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
            
            subroutine setogxs
                integer :: i, j, ii, k
                integer :: ierg
                real(8) :: erg
                real(8), allocatable :: rcvbuf(:,:)
                ! EFLUX option
                if(.not. do_burn) return
                if(depopt == 0) return
                gdelta = log10(Emax/Emin)/dble(ngrid)
                allocate(XS(num_iso*numrx, ngrid)); XS = 0D0
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(erg, ierg, i, j)
                do i = 1,num_iso
                do j = 1,numrx
                    if(depopt==2 .and. j/=1 .and. j/=7) cycle 
                    do ii = 1,ace(i)%NXS(4)
                        if(ace(i)%MT(ii) == RXMT(j)) then
                            do k = 1,ngrid
                            
                            erg = Emin*10.d0**((dble(k)-0.5d0)*gdelta)
                            call getierg(i,ierg,erg)
                            XS((i-1)*numrx+j,k) = getxs(iso=i,mt_ENDF=RXMT(j),erg=erg,ierg=ierg)
                            !if(icore==score) write(*,'(I6,I4,F10.3,F10.3,I6)') ace(i)%ZAID, ace(i)%MT(ii), XS((i-1)*numrx+j,k), erg, ierg
                            enddo
                        endif
                    enddo
                enddo
                if(icore==score)print *, 'XS BUILT for ',ace(i)%xslib
                enddo
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                
                !allocate(rcvbuf(num_iso*numrx,ngrid))
                !call MPI_ALLREDUCE(XS,rcvbuf,num_iso*numrx*ngrid,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, ierr)
                !XS = rcvbuf
                !deallocate(rcvbuf)
                if(icore==score) print *, 'XS pre-calculation completed'
            end subroutine

            function buildflux(iso, n, eflux_tmp) result(flux)
            implicit none
            integer, intent(in) :: iso
            real(8), intent(in) :: eflux_tmp(:)
            real(8) :: eflux(0:nueg)
            real(8) :: flux(0:n-1)
            integer :: i, n, idx
            real(8) :: g, erg

            eflux(0:nueg) = eflux_tmp(1:nueg+1)

            flux(1:n-1) = 0d0
            flux(0)   = eflux(0)
            idx = 1
            FLX_LOOP: do i = 1, nueg
                if( ueggrid(i) >= ace(iso) % E(ace(iso)%NXS(3)) ) then
                    !flux(n-1) = flux(n-1) + sum(eflux(i:nueg))
                    flux(n-1) = flux(n-1) + eflux(i)
                    exit FLX_LOOP
!                elseif( ace(iso) % UEG % Egrid(i) == 0) then
!                    cycle
                elseif( ueggrid(i) < ace(iso) % E(1)) then
                    flux(0) = flux(0) + eflux(i)
                else ! In between
                    !idx = ace(iso) % UEG % Egrid(i)
                    EDO: do
                        if(ueggrid(i) < ace(iso) % E(idx+1)) exit EDO
                        idx = idx + 1
                    enddo EDO
                    flux(idx) = flux(idx) + eflux(i)
                endif
            enddo FLX_LOOP

            end function

            function buildogxs_iso(iso, rx, flux)
            implicit none
            real(8) :: buildogxs_iso
            integer, intent(in) :: iso, rx
            real(8), intent(in) :: flux(:)
            integer :: mt
            mt = ace(iso) % MT(rx)
            select case(mt)
            case(N_GAMMA)
                buildogxs_iso = dot_product(ace(iso) % sigd, flux)
            case(N_FISSION)
                buildogxs_iso = dot_product(ace(iso) % sigf, flux)
            case default
                if(rx/=0) buildogxs_iso = dot_product(ace(iso) % sig_MT(rx) % cx, flux)
                if(mt==N_3N .and. ace(iso)%zaid==42099) then
                    print *, 'MT', buildogxs_iso
                    !print *, 'CX', ace(iso)%sig_MT(rx)%cx
                    !print *, 'FLX', flux
                endif
            end select
            endfunction

            function buildogxs_e2(iso, rx, flux_tmp, flux2_tmp)
            implicit none
            real(8) :: buildogxs_e2
            integer, intent(in) :: iso, rx
            real(8), intent(in) :: flux_tmp(:), flux2_tmp(:)
            real(8) :: flux(0:ace(iso)%NXS(3)), flux2(0:ace(iso)%NXS(3))
            integer :: mt, i
            flux(0:ace(iso)%NXS(3)) = flux_tmp(1:ace(iso)%NXS(3)+1)
            flux2(0:ace(iso)%NXS(3))= flux2_tmp(1:ace(iso)%NXS(3)+1)
            mt = ace(iso) % MT(rx)
            buildogxs_e2 = 0d0
            select case(mt)
            case(N_GAMMA)
                buildogxs_e2 = flux(0) * ace(iso) % sigd(1)
                do i = 1, ace(iso) % NXS(3)-1
                    buildogxs_e2 = buildogxs_e2 + &
                        flux(i) * ace(iso) % sigd(i) + &
                        (flux2(i)-ace(iso)%E(i)*flux(i)) * (ace(iso)%sigd(i+1)-ace(iso)%sigd(i))/ (ace(iso)%E(i+1)-ace(iso)%E(i))
                    !buildogxs_e2 = buildogxs_e2 + flux(i) * (ace(iso)%sigd(i)+ace(iso)%sigd(i+1))*5d-1
                enddo
                buildogxs_e2 = buildogxs_e2 + flux(ace(iso)%NXS(3)) * ace(iso)%sigd(ace(iso)%NXS(3))


            case(N_FISSION)
                buildogxs_e2 = flux(0) * ace(iso) % sigf(1)
                do i = 1, ace(iso) % NXS(3)-1
                    buildogxs_e2 = buildogxs_e2 + &
                        flux(i) * ace(iso) % sigf(i) + &
                        (flux2(i)-ace(iso)%E(i)*flux(i)) * (ace(iso)%sigf(i+1)-ace(iso)%sigf(i))/ (ace(iso)%E(i+1)-ace(iso)%E(i))
                    !buildogxs_e2 = buildogxs_e2 + flux(i) * (ace(iso)%sigf(i)+ace(iso)%sigf(i+1))*5d-1
                enddo
                buildogxs_e2 = buildogxs_e2 + flux(ace(iso)%NXS(3)) * ace(iso)%sigf(ace(iso)%NXS(3))
            case default
                if(rx == 0) return
                !do i = ace(iso) % sig_MT(rx) % IE, ace(iso) % sig_MT(rx) % IE + ace(iso) % sig_MT(rx) % NE-2
                buildogxs_e2 = flux(0) * ace(iso) % sig_MT(rx) % cx(1)
                do i = 1, ace(iso) % NXS(3)-1
                    buildogxs_e2 = buildogxs_e2 + &
                        flux(i) * ace(iso) % sig_MT(rx) % cx(i) + &
                        (flux2(i)-ace(iso)%E(i)*flux(i)) * (ace(iso)%sig_MT(rx)%cx(i+1)-ace(iso)%sig_MT(rx)%cx(i))/ (ace(iso)%E(i+1)-ace(iso)%E(i))
                    !buildogxs_e2 = buildogxs_e2 + flux(i) * (ace(iso)%sig_MT(rx)%cx(i)+ace(iso)%sig_MT(rx)%cx(i+1))*5d-1
                enddo
                buildogxs_e2 = buildogxs_e2 + flux(ace(iso)%NXS(3)) * ace(iso)%sig_MT(rx)%cx(ace(iso)%NXS(3))
            end select
            endfunction

            function buildogxs_bias(iso, rx, flux, flux2)
            implicit none
            real(8) :: buildogxs_bias
            integer, intent(in) :: iso, rx
            real(8), intent(in) :: flux(:), flux2(:)
            integer :: mt, i
            mt = ace(iso) % MT(rx)
            buildogxs_bias = 0d0
            select case(mt)
            case(N_GAMMA)
                do i = 1, ace(iso) % NXS(3)-1
                    buildogxs_bias = buildogxs_bias + &
                        flux(i) * (ace(iso) % sigd(i)-ace(iso)%sigd(i+1))/2d0 + &
                        (flux2(i)-ace(iso)%E(i)*flux(i)) * (ace(iso)%sigd(i+1)-ace(iso)%sigd(i))/ (ace(iso)%E(i+1)-ace(iso)%E(i))
                enddo

            case(N_FISSION)
                do i = 1, ace(iso) % NXS(3)-1
                    buildogxs_bias = buildogxs_bias + &
                        flux(i) * (ace(iso) % sigf(i)-ace(iso)%sigf(i+1))/2d0 + &
                        (flux2(i)-ace(iso)%E(i)*flux(i)) * (ace(iso)%sigf(i+1)-ace(iso)%sigf(i))/ (ace(iso)%E(i+1)-ace(iso)%E(i))
                enddo
            case default
                if(rx == 0) return
                do i = ace(iso) % sig_MT(rx) % IE, ace(iso) % sig_MT(rx) % IE + ace(iso) % sig_MT(rx) % NE-2
                    buildogxs_bias = buildogxs_bias + &
                        flux(i) * (ace(iso) % sig_MT(rx) % cx(i)-ace(iso) % sig_MT(rx) % cx(i+1))/2d0 + &
                        (flux2(i)-ace(iso)%E(i)*flux(i)) * (ace(iso)%sig_MT(rx)%cx(i+1)-ace(iso)%sig_MT(rx)%cx(i))/ (ace(iso)%E(i+1)-ace(iso)%E(i))
                enddo
            end select
            endfunction

            function buildogxs(iso, mt, eflux, toteflux)
            implicit none
            !include 'mkl.fi'
            real(8) :: buildogxs
            integer, intent(in) :: iso      !> Isotope
            integer, intent(in) :: mt       !> MT of reaction
            real(8), intent(in) :: eflux(:) !> Material-wise Fine spectrum
            real(8), intent(in) :: toteflux
            
            select case(mt)
            case(N_GAMMA)
                buildogxs = dot_product(ace(iso) % UEG % sigd(:), eflux)/toteflux
                !buildogxs = dot(ace(iso) % UEG % sigd(:), eflux) / toteflux
            case(N_2N)
                buildogxs = dot_product(ace(iso) % UEG % sig2n(:), eflux)/toteflux
            case(N_3N)
                buildogxs = dot_product(ace(iso) % UEG % sig3n(:), eflux)/toteflux
            case(N_4N)
                buildogxs = dot_product(ace(iso) % UEG % sig4n(:), eflux)/toteflux
            case(N_P)
                buildogxs = dot_product(ace(iso) % UEG % sigp(:), eflux)/toteflux
            case(N_A)
                buildogxs = dot_product(ace(iso) % UEG % sigal(:), eflux)/toteflux
            case(N_FISSION)
                buildogxs = dot_product(ace(iso) % UEG % sigf(:), eflux)/toteflux
            case default
                    
            end select
            end function

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
            
            subroutine ADDACE(zai,tail)
                integer,      intent(in) :: zai
                character(4), intent(in) :: tail
                character(6) :: zaline
                character(10) :: line
                integer :: i, za, zalen
                
                za = zai / 10
                write(zaline, '(i6)') za
                zaline = adjustl(zaline)
                if(mod(zai,10)>0) then
                    zalen = len(zaline)
                    zaline(zalen-3:zalen-3) = '3'
                endif

                line = trim(trim(zaline)//tail)
                do i = 1, size(libname)
                    if(acerecord(i) .and. trim(line) == trim(libname(i))) then
                        num_iso = num_iso + 1
                        ace(num_iso) % xslib   = trim(line)
                        ace(num_iso) % library = trim(libpath(i))
                        acerecord(i) = .false.
                        call set_ace_iso(num_iso, trim(line))
                        return
                    endif
                enddo
            end subroutine
                
            ! ===================================================================
            !         Depletion :: Make depletion matrix and solve 
            ! ===================================================================
            subroutine depletion 
            
            implicit none
            
            integer :: imat, jmem, jnuc, knuc, knuc1, knuc2
            integer :: mt_iso, iso, niso
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
            real(8) :: ratio
            real(8) :: remsum
            ! TESTING for SM-149
            real(8) :: samarium !Estimate Samarium
            real(8) :: samabs
            real(8) :: ingrid
            character(len=10) :: fileid, matid
            character(len=50) :: filename, directory

            ! Homogenized OGXS
            real(8) :: sigmaa
            real(8) :: nusigf
            integer :: ii, jj, i_rx
            ! MTRXREAD
            integer :: mt, nn, pn, dn, tn, an, a3n, rx
            real(8) :: ogxs, erg, ogxs1
            real(8),allocatable :: flx(:), flx2(:)
            integer :: ierg,idx
            real(8) :: toteflux
            real(8),allocatable :: tmpogxs(:)

            integer :: eg
            integer :: aceval

            logical :: sorted
            real(8) :: tmpnumden
            integer :: tmpaceidx

            real(8) :: totfiss, numer, denom, g2
            integer :: cnt, rcv, addn

            logical :: do_exist

            if(do_burn==.false.) return
            avg_power = avg_power / dble(n_act)
            tot_flux = 0.d0

            call MPI_BCAST(avg_power, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)
            if(icore==score) print *, 'Avg power[MeV]',avg_power 

            do imat = 1, n_materials
                if(.not. materials(imat)%depletable) cycle
                call MPI_BCAST(materials(imat)%eflux, 1+nueg, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD,ierr)
                call MPI_BCAST(materials(imat)%flux , 1   , MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)
            enddo

            ! NORMALIZE SFY too
            !allocate(nucexist(nsfp)); nucexist = 0.d0
            !do i = 1,nsfp
            !    zai = sfp_zai(i)
            !    anum = zai/10000; mnum = (zai-anum*10000)/10; nnum = mnum-anum;
            !    inum = zai-anum*10000-mnum*10;
            !    if(nuclide(inum,nnum,anum)%idx>0) nucexist(i) = 1.d0
            !enddo
            !do i = 1,nsfssn
            !    ratio = sum(sfy_yield(:,i)*nucexist)/2.d0
            !    if(ratio>0.d0) sfy_yield(:,i) = sfy_yield(:,i)/ratio
            !enddo
            !deallocate(nucexist)

            !Initialize material independent burnup matrix
            if(istep_burnup==0) then
            allocate(bMat(1:nnuc,1:nnuc))   !2-D burnup matrix : row to column transition
            allocate(bMat0(1:nnuc,1:nnuc)) !2-D burnup matrix : row to column transition (material independent)
            bMat0 = 0.d0; bMat = 0d0
            
            !Build material independent burnup matrix
            do jnuc = 1, nnuc
                anum = zai_idx(jnuc)/10000
                mnum = (zai_idx(jnuc) - anum*10000)/10
                nnum = mnum - anum
                inum = zai_idx(jnuc) - anum*10000 - mnum*10
                if(nnum<0) then
                    print *, zai_idx(jnuc)
                    cycle
                endif
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
                        do i = 1,nsfssn
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
!                if (nuclide(inum,nnum,anum)%a_emit>0.d0) then
!                    knuc = nuclide(0,2,2)%idx
!                    if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
!                        nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%a_emit
!                endif
!                if (nuclide(inum,nnum,anum)%n_emit>0.d0) then
!                    knuc = nuclide(0,1,0)%idx
!                    if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
!                        nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%n_emit
!                endif
!                if (nuclide(inum,nnum,anum)%p_emit>0.d0) then
!                    knuc = nuclide(0,0,1)%idx
!                    if (knuc>0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
!                        nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%p_emit
!                endif
                
                ! Count for removal per circulation
!                if (nuclide(inum,nnum,anum)%removal > 0.d0) then
!                    bMat0(jnuc,jnuc) = bMat0(jnuc,jnuc) + &
!                    log(1.d0-nuclide(inum,nnum,anum)%removal)/t_circulation
!                end if
            end do
            end if
            !print *, 'BU' 
            if (icore==score) then 
                write(prt_bumat, '(a45)')         '   =========================================='
                write(prt_bumat, '(a17,i4)')     '      Burnup step', istep_burnup+1
                write(prt_bumat, '(f14.2,a16)') burn_step(istep_burnup+1)/86400.d0, ' CUMULATIVE DAYS'
                write(prt_bumat, '(a45)')         '   =========================================='
            endif 
            
            !Normalization constant to be real power
            ULnorm = RealPower/(avg_power*eVtoJoule)
            call MPI_BCAST(ULnorm, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)
            !Substitute burnup matrix element
            !do imat = 1, n_materials
            cnt = 0
            do ii = 1, ngeom
                imat = mpigeom(ii,icore)

                if(imat==0) cycle
                if(.not. materials(imat)%depletable) cycle    !material imat is not burned
                !samarium = 0.d0
                mat => materials(imat)
                !print *, icore, 'MATDEP', imat, '/', totgeom
!                    write(prt_bumat, *) '' 
!                    write(prt_bumat, *) mat%mat_name 
!                    write(prt_bumat, *) ''
                !Call the material independent burnup matrix 
                
                allocate(yield_data(nfp,nfssn))
                yield_data = 0.d0
                !TODO NFY interpolation Option
                select case(NFYtype)
                case(1) ! Nuc.wise Avg Energy 
                do i = 1,nfssn
                    numer = 0d0; denom = 0d0;
                    totfiss = 0d0
                    zai = fssn_zai(i)
                    anum = zai/10000; mnum = (zai-10000*anum)/10; nnum = mnum-anum
                    inum = zai-anum*10000-mnum*10
                    aceval = find_ACE_iso_idx_zaid(zai)
                    mt_iso = 0
                    do j = 1, mat % n_iso
                        if(mat % ace_idx(j) == aceval) then
                            mt_iso = j; exit
                        endif
                    enddo
                    if(mt_iso==0) cycle
                    Ep = yieldE(i,1:4); nE = yieldnE(i)
                    if(.not. allocated(ace(aceval) % UEG % sigf)) cycle

                    numer = numer + mat % eflux(0) * ace(aceval) % UEG % sigf(1) * ueggrid(1) * mat % numden(mt_iso)
                    denom = denom + mat % eflux(0) * ace(aceval) % UEG % sigf(1) * mat % numden(mt_iso)

                    do j = 1, nueg
                        numer = numer + mat % eflux(j) * ace(aceval) % UEG % sigf(j) * ueggrid(j) * mat % numden(mt_iso)
                        denom = denom + mat % eflux(j) * ace(aceval) % UEG % sigf(j) * mat % numden(mt_iso)
                    enddo
                    
                    erg = numer/denom
                    !print *, 'NFY', fssn_zai(i), erg
                    
                    if(nE<=1) then
                        yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i)
                    else
                        if(erg < Ep(1)) then
                            yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i)
                        elseif(erg > Ep(nE)) then
                            yield_data(1:nfp,i) = tmp_yield(1:nfp,nE,i)
                        else
                            do eg = 1, nE-1
                                if(erg>=Ep(eg) .and. erg<Ep(eg+1)) then
                                    g2 = (erg-Ep(eg))/(Ep(eg+1)-Ep(eg))
                                    !print *, fssn_zai(i), erg, g2, Ep(eg), Ep(eg+1)
                                    yield_data(1:nfp,i) = &
                                        tmp_yield(1:nfp,eg,i) * (1d0-g2) + &
                                        tmp_yield(1:nfp,eg+1,i) * g2
                                endif
                            enddo
                        endif
                    endif
                enddo

                case(2) ! Material-wise Avg. Energy
                ! USING UNIFIED ENERGY for NFY Interpolation
                numer = 0d0; denom = 0d0
                do i = 1,nfssn
                    totfiss = 0d0
                    zai = fssn_zai(i)
                    anum = zai/10000; mnum = (zai-10000*anum)/10; nnum = mnum-anum
                    inum = zai-anum*10000-mnum*10
                    aceval = find_ACE_iso_idx_zaid(zai)
                    mt_iso = 0
                    do j = 1, mat % n_iso
                        if(mat % ace_idx(j) == aceval) then
                            mt_iso = j; exit
                        endif
                    enddo
                    if(mt_iso==0) cycle
                    Ep = yieldE(i,1:4); nE = yieldnE(i)
                    if(.not. allocated(ace(aceval) % UEG % sigf)) cycle

                    numer = numer + mat % eflux(0) * ace(aceval) % UEG % sigf(1) * ueggrid(0) * mat % numden(mt_iso)
                    denom = denom + mat % eflux(0) * ace(aceval) % UEG % sigf(1) * mat % numden(mt_iso)
                    do j = 1, nueg
                        numer = numer + mat % eflux(j) * ace(aceval) % UEG % sigf(j) * ueggrid(j) * mat % numden(mt_iso)
                        denom = denom + mat % eflux(j) * ace(aceval) % UEG % sigf(j) * mat % numden(mt_iso)
                    enddo
                enddo
                !print *, 'NFY', numer/denom
                erg = numer/denom
                !erg = 0.85355
                do i = 1, nfssn
                    Ep = yieldE(i,1:4); nE = yieldnE(i)
                    if(nE<=1) then
                        yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i)
                    else
                        if(erg < Ep(1)) then
                            yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i)
                        elseif(erg > Ep(nE)) then
                            yield_data(1:nfp,i) = tmp_yield(1:nfp,nE,i)
                        else
                            do eg = 1, nE-1
                                if(erg>=Ep(eg) .and. erg<Ep(eg+1)) then
                                    g2 = (erg-Ep(eg))/(Ep(eg+1)-Ep(eg))
                                    yield_data(1:nfp,i) = &
                                        tmp_yield(1:nfp,eg,i) * (1d0-g2) + &
                                        tmp_yield(1:nfp,eg+1,i) * g2
                                    !print *, 'POS', erg, Ep(eg), Ep(eg+1), g2
                                endif
                            enddo
                        endif
                    endif
                enddo
                case(3)
                do i = 1,nfssn
                    Ep = yieldE(i,1:4); nE = yieldnE(i)
                    zai = fssn_zai(i)
                    if(nE==0) then
                        yield_data(1:nfp,i) = tmp_yield(1:nfp,1,i)
                    else
                        if(nE==1) then
                            aceval = find_ACE_iso_idx_zaid(zai)
                            mat%fratio(i,1) = 1.d0
                        else
                            aceval = find_ACE_iso_idx_zaid(zai)
                            if(aceval==0) cycle
                            if(.not. allocated(ace(aceval) % UEG % sigf)) cycle
                            do j = 1,nueg
                                erg = ueggrid(j)
                                if(erg<Ep(1)) then
                                    mat%fratio(i,1) = mat%fratio(i,1) + mat%eflux(j) * ace(aceval) % UEG % sigf(j)
                                elseif(erg>=Ep(nE)) then
                                    mat%fratio(i,nE) = mat%fratio(i,nE) + mat%eflux(j) * ace(aceval) % UEG % sigf(j)
                                else
                                    do eg = 1,nE-1
                                        if(erg>=Ep(eg) .and. erg<Ep(eg+1)) then
                                            g2 = (erg-Ep(eg))/(Ep(eg+1)-Ep(eg))
                                            mat%fratio(i,eg) = mat%fratio(i,eg) + mat%eflux(j) * ace(aceval) % UEG % sigf(j)*(1D0-g2)
                                            mat%fratio(i,eg+1) = mat%fratio(i,eg+1) + mat%eflux(j) * ace(aceval) % UEG % sigf(j)*g2
                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                        if(sum(mat%fratio(i,1:nE))>0) then
                            mat%fratio(i,1:nE) = mat%fratio(i,1:nE)/sum(mat%fratio(i,1:nE))
                        else
                            mat%fratio(i,1) = 1.d0
                            mat%fratio(i,2:4) = 0.d0
                        endif
                        do k = 1,nE
                            yield_data(1:nfp,i) = yield_data(1:nfp,i) + &
                                tmp_yield(1:nfp,k,i) * mat%fratio(i,k)
                        enddo
                    endif
                enddo
                end select

                allocate(nucexist(nfp)); nucexist = 0.d0
                do i = 1,nfp
                    ! PRIOR to normalization; exclude non-existing nuclides in FPY/SFY
                    zai = fp_zai(i)
                    anum = zai/10000; mnum = (zai-anum*10000)/10; nnum = mnum-anum;
                    inum = zai-anum*10000-mnum*10;
                    if(nuclide(inum,nnum,anum)%idx>0) nucexist(i) = 1.d0
                enddo

                !do i = 1,nfssn
                !    ! NORMALIZE NFY: sum(NFY) = 200
                !    ratio = sum(yield_data(:,i)*nucexist)/2.d0
                !    !if(ratio>0.d0) yield_data(:,i) = yield_data(:,i)/ratio
                !    !if(icore==score) print *, 'NFY', fssn_zai(i), ratio
                !enddo
                deallocate(nucexist)
                !Calculate real flux (volume-averaged)
                real_flux = ULnorm*mat%flux
                print *, 'REAL FLUX', trim(mat%mat_name), real_flux, ULnorm, mat%flux
                !$OMP ATOMIC
                tot_flux = tot_flux + real_flux*mat%vol
                toteflux = sum(mat%eflux(0:nueg))
                !Build burnup matrix with cross section obtained from MC calculation
                bMat = bMat0*bstep_size
                !if(icore==score) print *, 'bMat0', bMat(nnuc,:), bMat0(nnuc,:)
                !!$omp parallel do default(private) &
                !!$omp shared(bMat, real_flux, toteflux, bstep_size, yield_data, ZAIMT_ism, gnd_frac, ace, nuclide, fssn_zai, fp_zai, nfssn, RXMT, num_iso)
                !$OMP PARALLEL DO &
                !$OMP PRIVATE(iso, anum, mnum, nnum, inum, jnuc, flx, flx2, mt, ogxs, ogxs1, fy_midx, anum1, nnum1, mnum1, inum1, knuc, a1, m1, n1, ism, zai, pn, dn, tn, an, a3n, nn, addn)  
                DO_ISO: do mt_iso = 1,num_iso
                    !print *, 'ISO', mt_iso
                    iso = mt_iso
                    anum = ace(iso)%zaid/1000
                    mnum = (ace(iso)%zaid - anum*1000)
                    inum = 0
                    if(mnum>300) then
                        inum = 1
                        mnum = mnum - 200
                        if(anum>88) mnum = mnum + 100
                    elseif(mnum==0) then
                        cycle
                    endif
                    nnum = mnum-anum

                    jnuc = nuclide(inum,nnum,anum)%idx
                    if(jnuc==0) cycle

                    ! BUILD ISO-WISE FLUX
                    if(do_ueg) then
                        flx  = buildflux(iso,ace(iso)%NXS(3)+1, mat %eflux(0:nueg))
                        flx2 = buildflux(iso,ace(iso)%NXS(3)+1, mat%e2flux(0:nueg))
                        flx  = flx / toteflux; flx = flx * real_flux
                        flx2 = flx2/ toteflux; flx2= flx2* real_flux
                    endif
                    
                    do rx = 1,ace(iso)%NXS(4)
                        mt = ace(iso)%MT(rx)
                        !if(abs(ace(iso)%TY(rx))==1) print *, 'EXCLUDED', iso, mt
                        !if(icore==score .and. ace(iso)%zaid==57138) print *, 'MTS', rx, mt
                        if(abs(ace(iso)%TY(rx))==1 .and. .not.(mt==N_NA .or. mt==N_NF .or. mt==N_NA .or. mt==N_N3A .or. mt==N_NP .or. mt==N_N2A .or. mt==N_ND .or. mt==N_NT .or. mt==N_N3HE .or. mt==N_ND2A .or. mt==N_NT2A .or. mt==N_N2P .or. mt==N_NPA .or. mt==N_NDA .or. mt==N_NPD .or. mt==N_NPT .or. mt==N_NDT .or. mt==N_NP3HE .or. mt==N_ND3HE .or. mt==N_NT3HE .or. mt==N_NTA .or. mt==N_N3P)) cycle ! Maybe inelastic?
                        if(mt == 4 .or. mt > 200) cycle ! Inelastic and Damage
                        !if(abs(ace(iso)%TY(rx))==1 .and. icore==score) print *, 'WOW',ace(iso)%zaid, mt
                        ! TALLY OGXS
                        !elseif(
                        if(do_ueg) ogxs = buildogxs_e2(iso, rx, flx, flx2) * barn
                        if(ace(iso)%zaid==92230) print *, 'U230', rx, mt, ogxs

                        if(do_rx_tally) then
                            if(ANY(RXMT==mt)) then
                                ogxs1 = ogxs
                                do i_rx = 1, 7
                                    if(RXMT(i_rx)==mt) then
                                        ogxs = mat % ogxs(iso, i_rx) * real_flux
                                        if(ace(iso)%zaid==92235 .or. ace(iso)%zaid==92238) write(*,'(A,A,I6, I3,E15.5,E15.5,E15.5)') 'BIAS ', trim(materials(imat)%mat_name), ace(iso)%zaid, mt, ogxs, ogxs1, (ogxs1-ogxs)/ogxs*1E2
                                        exit
                                    endif
                                enddo
                            endif
                        endif
!                        if(ace(iso)%zaid==57138) then
!                            print *, trim(materials(imat)%mat_name), ace(iso)%zaid, mt, ogxs
!                            call mtrxread(mt, nn, pn, dn, tn, an, a3n)
!                            print *, 'MT', 'N', 'P', 'D', 'T', 'A', 'A3'
!                            print *, nn, pn, dn, tn, an, a3n
!                        endif
                        !if(mt==N_NF .or. mt==N_2NF .or. mt==N_3NF) cycle
                        ! FIND DESTINATION
                        if(mt==18 .or. ace(iso)%TY(rx)==19) then ! In case of Fission
                            if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then 
                                addn = 0
                                if(ace(iso)%MT(rx)==N_NF) then ! NNF
                                    addn = 1
                                elseif(ace(iso)%MT(rx)==N_2NF) then
                                    addn = 2
                                elseif(ace(iso)%MT(rx)==N_3NF) then
                                    addn = 3
                                endif
                                if(nuclide(inum,nnum-addn,anum)%fy_idx>0) then
                                    fy_midx = nuclide(inum,nnum-addn,anum)%fy_idx
                                else
                                    fy_midx = 0
                                    if(inum>1 .and. nuclide(0,nnum-addn,anum)%fy_idx>0) then
                                        fy_midx = nuclide(0,nnum-addn,anum)%fy_idx
                                    elseif(inum==0 .and. nuclide(1, nnum-addn, anum) % fy_idx > 0) then 
                                        fy_midx = nuclide(1, nnum-addn, anum) % fy_idx
                                    else
                                        do i = 1,nfssn
                                            a1 = fssn_zai(i)/10000
                                            m1 = (fssn_zai(i)-a1*10000)/10
                                            n1 = m1-a1
                                            if(n1==nnum-addn) fy_midx = i
                                        enddo
                                    if(fy_midx==0) fy_midx = nuclide(0,143,92)%fy_idx
                                    endif
                                endif
                                !if(ace(iso)%zaid==92230) print *, 'FY', fssn_zai(fy_midx)
                                do i = 1,nfp
                                    anum1 = fp_zai(i)/10000
                                    mnum1 = (fp_zai(i)-anum1*10000)/10
                                    nnum1 = mnum1 - anum1
                                    inum1 = fp_zai(i)-anum1*10000-mnum1*10
                                    knuc = nuclide(inum1,nnum1,anum1)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                        + ogxs * yield_data(i,fy_midx) * bstep_size

                                    !if(anum1==42 .and. mnum1==97) print *, 'MO97', fssn_zai(fy_midx), yield_data(i,fy_midx), tmp_yield(i,:,fy_midx) 
                                    !if(knuc/=0) print *, 'FP', knuc, jnuc, bMat(knuc,jnuc)
                                    !if(anum == 94 .and. nnum == 239-anum .and. icore==score) &
                                        !print *, 'FY', anum1, mnum1, yield_data(i,fy_midx)*ogxs,tmp_yield(i,:,fy_midx) 
                                enddo
                                !if(icore==score) print *, 'TSTING', jnuc,  bMat(nnuc,jnuc)
                                !print *, 'SUMFY', fy_midx, fssn_zai(fy_midx), sum(yield_data(1:nfp,fy_midx))

                                bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - ogxs * bstep_size
                            endif
                        else ! NON-fission reaction
                            if(mt==102) then ! (n,gamma), isomeric branching
                                ism = 0
                                zai = anum*10000 + mnum*10 + inum
                                do i = 1,num_brn
                                    if(zai==ZAIMT_ism(i)) then
                                        ism = i
                                        exit
                                    endif
                                enddo
                                if(ism>0) then
                                    knuc = nuclide(0,nnum+1,anum)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) &
                                        = bMat(knuc,jnuc) + ogxs * bstep_size * gnd_frac(ism)
                                    !if(knuc/=0) print *, 'NGI', knuc, jnuc, bMat(knuc,jnuc)
                                    knuc = nuclide(1,nnum+1,anum)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) &
                                        = bMat(knuc,jnuc) + ogxs * bstep_size * (1.d0-gnd_frac(ism))
                                    !if(knuc/=0) print *, 'NG', knuc, jnuc, bMat(knuc,jnuc)
                                else
                                    knuc = nuclide(0,nnum+1,anum)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) &
                                        = bMat(knuc,jnuc) + ogxs * bstep_size
                                    !if(knuc/=0) print *, 'NG', knuc, jnuc, bMat(knuc,jnuc)
                                endif
                                bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - ogxs * bstep_size
                            else
                                call mtrxread(mt,nn,pn,dn,tn,an,a3n)
                                anum1 = anum - pn - dn - tn - 2*an - 2*a3n
                                nnum1 = nnum + 1 - nn - dn - 2*tn - 2*an - a3n
                                if(anum1 == anum .and. nnum1 == nnum) cycle
                                knuc = 0
                                if(anum1>0 .and. nnum1>0) knuc = nuclide(0,nnum1,anum1)%idx
                                !do i = 1, nnuc
                                !    if(icore==score .and. bMat(nnuc,i)/=0) print *, 'RXNN1', rx, ace(iso)%MT(rx), i, bMat(nnuc,i)
                                !enddo
                                !if(knuc/=0) print *, 'MTR', knuc, jnuc, mt, bMat(knuc,jnuc)
                                if(knuc/=0 .and. (nn+pn+dn+tn+an+a3n)>0) then
                                    bMat(knuc,jnuc) = bMat(knuc,jnuc) + ogxs * bstep_size
                                    !do i = 1, nnuc
                                    !    if(icore==score .and. bMat(nnuc,i)/=0) print *, 'RXNN2', rx, ace(iso)%MT(rx), i, bMat(nnuc,i)
                                    !enddo
                                    !if(knuc/=0) print *, 'RX', knuc, jnuc, anum1, nnum1,  bMat(knuc,jnuc)
                                    !if(knuc/=0) print *, 'N', nn, pn, dn, tn, an, a3n
                                    if(nn>0) then
                                    knuc = nuclide(0,1,0)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                                        ogxs * bstep_size * nn
                                    !if(knuc/=0) print *, 'NN', knuc, jnuc, bMat(knuc,jnuc)
                                    endif
                                    if(pn>0) then
                                    knuc = nuclide(0,0,1)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                                        ogxs * bstep_size * pn
                                    !if(knuc/=0) print *, 'PN', knuc, jnuc, bMat(knuc,jnuc)
                                    endif
                                    if(dn>0) then
                                    knuc = nuclide(0,1,1)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                                        ogxs * bstep_size * dn
                                    !if(knuc/=0) print *, 'DN', knuc, jnuc, bMat(knuc,jnuc)
                                    endif
                                    if(tn>0) then
                                    knuc = nuclide(0,2,1)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                                        ogxs * bstep_size * tn
                                    !if(knuc/=0) print *, 'TN', knuc, jnuc, bMat(knuc,jnuc)
                                    endif
                                    if(an>0) then
                                    knuc = nuclide(0,2,2)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                                        ogxs * bstep_size * an
                                    endif
                                    if(a3n>0) then
                                    knuc = nuclide(0,1,2)%idx
                                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) + &
                                        ogxs * bstep_size * a3n
                                    !if(knuc/=0) print *, 'A3', knuc, jnuc, bMat(knuc,jnuc)
                                    endif
                                endif
                                bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - ogxs * bstep_size
                            endif
                        endif
                        !if(icore==score) print *, 'TESTINGN', jnuc, mt,  bMat(nnuc,jnuc)
                    enddo
                    !if(icore==score) print *, 'TESTING', jnuc, bMat(nnuc,jnuc)
                end do DO_ISO
                !$omp end parallel do
                !print *, 'END DOISO', imat, materials(imat) % mat_name, icore

        !if(icore==score) print *, 'bMat', bMat(nnuc,:)
        deallocate(yield_data)

        ! WRITE BURNUP MATRIX (OPTIONAL)
        if(bumat_print)then
        
        write(fileid,'(i3)') istep_burnup
        directory = './BUMAT/'//adjustl(trim(fileid))
        !directory = './BUMAT/UEG'
        call execute_command_line('mkdir -p '//adjustl(trim(directory)))
        filename = 'mat_'//trim(adjustl(mat%mat_name(:)))//'_step_'//trim(adjustl(fileid))//'.m'
!        idx = 0
!        do
!            idx = idx + 1
!            write(fileid,'(i3)') idx
!            filename = trim(adjustl(mat%mat_name(:)))//'_trial'//trim(adjustl(fileid))//'.m'
!            inquire(file=trim(directory)//'/'//trim(filename),exist=do_exist)
!            if(.not. do_exist) then
!                exit
!            endif
!        enddo

        open(bumat_test, file = trim(directory)//'/'//trim(filename),action="write",status="replace")
        write(bumat_test,*) 'FLUX=',real_flux,';'
        write(bumat_test,*) 'ZAI1=zeros(',nnuc,',1);'
        do knuc = 1,nnuc
            write(bumat_test,*) 'ZAI1(',knuc,')=',zai_idx(knuc),';'
        enddo
        write(bumat_test,*) 'A1 = zeros(',nnuc,',',nnuc,');'
        do knuc1 = 1,nnuc
            do knuc2 = 1,nnuc
                if(bMat(knuc2,knuc1)/=0.d0) write(bumat_test,*) 'A1(',knuc2,',',knuc1,')=',bMat(knuc2,knuc1)/bstep_size,';'
            enddo
        enddo
        endif
        !print *, 'B4 CRAM', imat, materials(imat) % mat_name, icore
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
        !print *, 'DONE', icore, imat, totgeom
!        if(totgeom < 10) then
!            print *, 'DEP solved for mat:',imat,'/',totgeom, ':', icore
!        else
!            !$OMP ATOMIC
!            cnt = cnt + 1
!
!            call MPI_ALLREDUCE(cnt,rcv,1,MPI_INTEGER,MPI_SUM,core,ierr)
!            cnt = rcv
!
!            if(icore==score) then
!
!            if(cnt == totgeom/4) then
!                print *, 'DEPLETION CALC ( 25%/100%)'
!            elseif(cnt == totgeom/4*2) then
!                print *, 'DEPLETION CALC ( 50%/100%)'
!            elseif(cnt == totgeom/4*3) then
!                print *, 'DEPLETION CALC ( 75%/100%)'
!            elseif(cnt == totgeom) then
!                print *, 'DEPLETION CALC (100%/100%)'
!            endif
!            
!            endif
!        endif

        !print *, 'COUNT?', icore, imat, totgeom
        
        ! WRITE ATOMIC DENSITY IN MATLAB .m FILE (OPTIONAL)
        if(bumat_print)then
        write(bumat_test,*) 'Ncomp=zeros(',nnuc,',2);'
        do knuc = 1,nnuc
            if(mat%full_numden(knuc)>0 .or. nxt_full_numden(knuc)>0) then
                write(bumat_test,*) 'Ncomp(',knuc,',1)=',mat%full_numden(knuc)*barn,';'
                write(bumat_test,*) 'Ncomp(',knuc,',2)=',nxt_full_numden(knuc)*barn,';'
            endif
            !write(*,*) 'X', imat, zai_idx(knuc), mat%full_numden(knuc)*barn, nxt_full_numden(knuc)*barn, real_flux
        enddo
        endif
            
        !Update number density
        mat%full_numden = nxt_full_numden 
        
        !print *, 'FLAG1', imat, icore
        ! ======================================================================================
        ! Reset material isotope inventory
        knuc = 0 
        remsum = 0.d0
        iso_idx = 0
        do jnuc=1, nnuc
            !write(*,*) imat, icore, zai_idx(jnuc), mat%full_numden(jnuc)
            tmp = find_ACE_iso_idx_zaid(zai_idx(jnuc))
            print *, "MAT", imat, jnuc, zai_idx(jnuc), tmp
            !if(mat%full_numden(jnuc)>0.d0 .and. tmp > 0) then 
            if(tmp>0 .and. mat%full_numden(jnuc)>1d0) then
                knuc = knuc + 1
                iso_idx(knuc) = jnuc
            elseif(mat%full_numden(jnuc)>0.d0) then
                ! MEASURING LOST TERM
                remsum = remsum + mat%full_numden(jnuc)
            endif
        end do

        !print *, 'FLAG2', imat, icore

        mat%n_iso = knuc
        deallocate(mat%ace_idx);allocate(mat%ace_idx(1:knuc));     mat%ace_idx= 0D0
        deallocate(mat%numden); allocate(mat%numden(1:knuc));      mat%numden = 0D0
        deallocate(mat%ogxs);   allocate(mat%ogxs(1:num_iso,1:7)); mat%ogxs   = 0d0

        i = 0 
        do mt_iso=1, mat % n_iso
            ! find ace_idx
            tmp = find_ACE_iso_idx_zaid(zai_idx(iso_idx(mt_iso)))
            if (tmp /= 0 .and. mat%full_numden(iso_idx(mt_iso))>1d0) then 
                i = i + 1
                mat%ace_idx(mt_iso) = tmp
                mat%numden(mt_iso)  = mat%full_numden(iso_idx(mt_iso))
                !print *, mat%ace_idx(mt_iso), ace(mat%ace_idx(mt_iso))%zaid, zai_idx(iso_idx(mt_iso)), iso_idx(mt_iso), mt_iso
            endif
        end do

!        do mt_iso = 1, mat%n_iso
!            write(*,*) 'Before',imat, mt_iso, mat%numden(mt_iso), mat%ace_idx(mt_iso)
!        enddo
        ! SORT FOR EFFICIENCY
!        sorted = .false.
!        if(mat%n_iso>1) then
!            do while ( .not. sorted )
!                sorted = .true.
!                do mt_iso = 1, mat%n_iso-1
!                    if(mat%numden(mt_iso) < mat%numden(mt_iso+1)) then
!                        if(sorted) sorted = .false.
!                        tmpnumden = mat%numden(mt_iso+1)
!                        tmpaceidx = mat%ace_idx(mt_iso+1)
!
!                        mat%numden(mt_iso+1)  = mat%numden(mt_iso)
!                        mat%ace_idx(mt_iso+1) = mat%ace_idx(mt_iso)
!                        mat%numden(mt_iso)    = tmpnumden
!                        mat%ace_idx(mt_iso)   = tmpaceidx
!                    endif
!                enddo
!            enddo
!        endif
                    
    !print *, 'FLAG3', imat, icore
            

    enddo
    !print *, 'ARRIVED', icore
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(icore==score) print *, 'Total Flux', tot_flux



    ! data sharing
    do ii = 0, ncore-1
    do jj = 1, ngeom
        imat = mpigeom(jj,ii)
        if( imat == 0 ) cycle
        mat => materials(imat)
        call MPI_BCAST(mat%n_iso,1,MPI_INTEGER,ii,MPI_COMM_WORLD,ierr)
        niso = mat%n_iso
    
        if ( icore /= ii ) then
           deallocate(mat%ace_idx); allocate(mat%ace_idx(niso)); mat%ace_idx(:) = 0
           deallocate(mat%numden);  allocate(mat%numden(niso));  mat%numden(:)  = 0
       end if
    
       call MPI_BCAST(mat%ace_idx,niso,MPI_INTEGER,ii,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(mat%numden,niso,MPI_REAL8,ii,MPI_COMM_WORLD,ierr)
    end do
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if(icore==score) then
        do imat = 1,n_materials
        if(.not. materials(imat)%depletable) cycle
            mat => materials(imat)
            write(prt_bumat,*) ''
            write(prt_bumat,*) 'mat: ', mat%mat_name
            write(prt_bumat,*) ''
            do mt_iso = 1, mat%n_iso
               write(prt_bumat, '(a15,e14.5)') ace(mat%ace_idx(mt_iso))%xslib, mat%numden(mt_iso)*barn 
            enddo
            write(prt_bumat, *) 'Num isotope', i
            write(prt_bumat, *) ''
        enddo
    endif

    do i = 1,n_materials
        mat => materials(i)
        print *, 'ISO', iso_idx
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
     
    end subroutine depletion 
    
    subroutine mtrxread(mt,nn,pn,dn,tn,an,a3n)
    integer, intent(in) :: mt
    integer, intent(out):: nn,pn,dn,tn,an,a3n
    nn = 0; pn = 0; dn = 0; tn = 0; an = 0; a3n = 0
    if(mt==N_P) then !NP RX
        pn = 1
    elseif(mt==N_D) then !ND RX
        dn = 1
    elseif(mt==N_T) then !NT RX
        tn = 1
    elseif(mt==N_3HE) then !N3HE RX
        a3n = 1
    elseif(mt==N_A) then !NA RX
        an = 1
    elseif(mt==N_2A) then !N 2A RX
        an = 2
    elseif(mt==N_3A) then !N 3A RX
        an = 3
    elseif(mt==N_2P) then
        pn = 2
    elseif(mt==N_PA) then
        pn = 1; an = 1
    elseif(mt==N_T2A) then
        tn = 1; an = 2
    elseif(mt==N_D2A) then
        dn = 1; an = 2
    elseif(mt==N_PD) then
        pn = 1; dn = 1
    elseif(mt==N_PT) then
        pn = 1; tn = 1
    elseif(mt==N_DA) then
        dn = 1; an = 1
    elseif(mt==N_5N) then
        nn = 5
    elseif(mt==N_2ND) then
        nn = 2; dn = 1
    elseif(mt==N_2NP) then
        nn = 2; pn = 1
    elseif(mt==N_2N) then
        nn = 2
    elseif(mt==N_3N) then
        nn = 3
    elseif(mt==N_NA) then
        nn = 1; an = 1
    elseif(mt==N_N3A) then
        nn = 1; a3n = 1
    elseif(mt==N_2NA) then
        nn = 2; an = 1
    elseif(mt==N_3NA) then
        nn = 3; an = 1
    elseif(mt==N_NP) then
        nn = 1; pn = 1
    elseif(mt==N_N2A) then
        nn = 1; an = 2
    elseif(mt==N_2N2A) then
        nn = 2; an = 2
    elseif(mt==N_ND) then
        nn = 1; dn = 1
    elseif(mt==N_NT) then
        nn = 1; tn = 1
    elseif(mt==N_N3HE) then
        nn = 1; a3n = 1
    elseif(mt==N_ND2A) then
        nn = 1; dn = 1; an = 2
    elseif(mt==N_NT2A) then
        nn = 1; tn = 1; an = 2
    elseif(mt==N_4N) then
        nn = 4
    endif
    end subroutine
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
        real(8) :: Ep(4)
        integer :: nE, eg

        type (Material_CE), pointer :: mat
        real(8) :: micro_f, micro_d, micro_a, ipfac
        real(8) :: micro_2n, micro_3n, micro_4n, micro_p
        !real(8) :: dep_ogxs(7)
        integer :: anum, mnum, inum, nnum !atomic number, mass number, isomer state, neutron number of nuclide
        integer :: a1, m1
        real(8) :: fluxtmp, val1, val2
        integer :: fssn

        real(8), pointer :: ogxs(:,:)

        if(E_mode==0) return       !Multigroup  -> return
        if(do_burn==.false.) return     !No burnup calculation -> return
        if(curr_cyc <= n_inact) return 
        if(.not. materials(imat)%depletable) return ! not depletable -> return
        
        if(erg > ueggrid(nueg) .or. erg < ueggrid(1)) return
        
        mat => materials(imat)
        
        !$omp atomic
        mat%flux = mat%flux + wgt*distance !Volume-integrated flux tally (not normalized)
        if(do_ueg) then
            ! EFLUX TALLY
            call getiueg(erg, ierg)
            
            
            if(ierg>=0) then
   
                !$OMP ATOMIC
                mat%eflux(ierg) = mat%eflux(ierg) + wgt * distance
        
                !$OMP ATOMIC
                mat%e2flux(ierg)= mat%e2flux(ierg)+ wgt * distance * max(ueggrid(1), min(ueggrid(nueg), erg))

            endif
        endif
        
        if(do_rx_tally)then
            do iso = 1,num_iso
               call getierg(iso,ierg,erg)
               ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
               micro_d = ace(iso)%sigd(ierg) + ipfac*(ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
               micro_f = 0.d0
               !FISSIONABLE
               if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
                   micro_f = ace(iso)%sigf(ierg) + ipfac*(ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
               endif
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
        endif
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
            materials(imat)%eflux = 0D0
            materials(imat)%e2flux= 0d0
            !materials(imat)%pwr = 0.0d0 
            materials(imat)%ogxs(:,:) = 0.0d0
            !allocate(materials(imat)%full_numden(nnuc))
            materials(imat)%fratio = 0.d0
        enddo
        ! TESTING for NFY

        bstep_size = burn_step(istep_burnup+1) - burn_step(istep_burnup) 
        avg_power = 0.d0
        tot_mass = 0.d0; tot_fmass = 0.d0
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
            allocate(materials(imat)%eflux(0:nueg))
            allocate(materials(imat)%e2flux(0:nueg))
            materials(imat)%eflux(:) = 0.d0
            materials(imat)%e2flux(:)= 0d0
            
            if (materials(imat)%depletable) then
                allocate(materials(imat)%full_numden(1:nnuc))
                materials(imat)%full_numden = 0.d0     !Initialize number density
            endif
            
            do mt_iso = 1, materials(imat)%n_iso
                iso = materials(imat)%ace_idx(mt_iso)
                inum = 0
                anum = ace(iso)%zaid/1000
                mnum = ace(iso)%zaid - anum*1000
                if(mnum>300) then
                    inum = 1
                    mnum = mnum - 200
                    if(anum>88) mnum = mnum + 100
                elseif(mnum==0) then
                    cycle
                endif
                nnum = mnum - anum

                if(materials(imat)%depletable==.true. .and. do_burn .and. nuclide(inum,nnum,anum)%idx>0) then 
                    materials(imat)%full_numden(nuclide(inum,nnum,anum)%idx) = materials(imat)%numden(mt_iso)
                endif 
            end do
            materials(imat)%fratio = 0.d0
        enddo 
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

            val = dble(n_act) * materials(imat)%flux * materials(imat)%vol

            ! EFLUX
            allocate(rcvbufarrlong(0:nueg)); allocate(sndbufarrlong(0:nueg))
            sndbufarrlong(0:nueg) = materials(imat)%eflux(0:nueg)
            call MPI_ALLREDUCE(sndbufarrlong, rcvbufarrlong, 1+nueg, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)
            rcvbufarrlong(0:nueg) = rcvbufarrlong(0:nueg)/ val
            materials(imat)%eflux = rcvbufarrlong
            deallocate(rcvbufarrlong); deallocate(sndbufarrlong)
            
            allocate(rcvbufarrlong(0:nueg)); allocate(sndbufarrlong(0:nueg))
            sndbufarrlong(0:nueg) = materials(imat)%e2flux(0:nueg)
            call MPI_ALLREDUCE(sndbufarrlong, rcvbufarrlong, 1+nueg, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)
            rcvbufarrlong(0:nueg) = rcvbufarrlong(0:nueg)/ val
            materials(imat)%e2flux = rcvbufarrlong
            deallocate(rcvbufarrlong); deallocate(sndbufarrlong)

            if(do_rx_tally) then
                ! DIRECT RX rate TALLY
                allocate(rcvbufarrlong(1:num_iso*7))
                allocate(sndbufarrlong(1:num_iso*7))
    
                val = dble(n_act) * materials(imat)%flux * materials(imat)%vol
                do iso = 1,num_iso
                    sndbufarrlong((iso-1)*7+1:(iso-1)*7+7) = materials(imat)%ogxs(iso,:)
                enddo 
                call MPI_ALLREDUCE(sndbufarrlong, rcvbufarrlong, num_iso*7, MPI_DOUBLE_PRECISION, MPI_SUM, core, ierr)                    
                rcvbufarrlong(:) = rcvbufarrlong(:) / val
                do iso = 1,num_iso
                    materials(imat)%ogxs(iso,:) = rcvbufarrlong((iso-1)*7+1:(iso-1)*7+7) 
                enddo
                deallocate(rcvbufarrlong)
                deallocate(sndbufarrlong) 
            endif
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
