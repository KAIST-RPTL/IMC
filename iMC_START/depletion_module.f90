
include 'mkl_pardiso.f90'

module depletion_module 
    use FMFD_HEADER, only: DTMCBU, buflux, bukapa
    use variables
    use constants
    use ace_header, only: ace, find_ACE_iso_idx_zaid, num_iso
    use ace_xs,     only: getxs, getierg
    use material_header
    use mpi 
    
    implicit none 
    character(80):: bufile
    
    !==============================================================================
    !Nuclide data structure for depletion
    type nuclide_data
      logical :: data_exist !data is read from library
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
      real(8) :: fb    !Fraction of beta- decay to ground state
      real(8) :: fbx   !Fraction of beta- decay to excited state
      real(8) :: fpec  !Fraction of all beta+ decay or electron capture
      real(8) :: fpecx !Fraction of beta+ decay or electron capture to excited state
      real(8) :: fa    !Fraction of all alpha decay
      real(8) :: fit   !Fraction of gamma-ray release from excited state to ground state
      real(8) :: fsf   !Fraction of all decay events of spontaneous fission
      real(8) :: fn    !Fraction of all decay events that are (beta- & neutron) decays
      real(8) :: qrec  !Total recoverable energy for all decay events in MeV
      real(8) :: abund !Naturally occuring isotopic abundance of NUCLID in atom percent
      real(8) :: arcg  !Radioactivity concentration guide (RCG) for continuous inhaliation
      real(8) :: wrcg  !Radioactivity concentration guide (RCG) for continuous ingestion
    end type
    !nuclide data for Isomeric State 
    ! (0=ground, 1=excited), neutron number, atomic number 
    type(nuclide_data) :: nuclide(0:1,0:160,1:100) 


    !Nuclide data index
    integer :: zai_idx(1:4000) !ZAI index = Z*10000 + A*10 + I
    integer :: maxnnz          !Maximum allowed number of nonzero elements in burnup matrix
    integer :: nnuc            !Total number of nuclides
    integer :: inuc            !Current nuclide index
    real(8) :: tmass           !Total mass
    real(8) :: fmass           !Fuel mass

    ! Material index for parallel calculation
    integer, allocatable:: mp1(:), mp2(:)

    !Fission yield data
    integer :: nfssn           !number of fissionable nuclide
    integer :: nfp             !number of fission product
    integer, allocatable :: fssn_zai(:)     !fissionable zai index
    integer, allocatable :: fp_zai(:)       !fission product zai index
    real(8), allocatable :: yield_data(:,:) !fission yield data (1:nfp, 1:nfssn)


    !Burnup matrix
    real(8), allocatable :: bMat(:,:)  !2-D burnup matrix
    real(8), allocatable :: bMat0(:,:) !2-D burnup matrix (material independent)
    !real(8), allocatable :: bMat_tmp(:,:) 


    !Real Power to normalize flux
    real(8) :: RealPower  !MWth
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
    
    !==============================================================================
    !In-line xenon equilibrium search for depletion calculation
    logical :: Xe_search=.true.               !Inline Xenon and Iodine search
    logical :: do_Xe_search                   !Inline Xenon and Iodine search
    integer :: Xe135_iso                      !Index for Xe135 in global list
    integer :: Xe_skp                         !Number of skipped iteration for Xenon tally
    integer :: Xe_sfifo                       !Length of fifo queue length for Xe_tally
    integer, allocatable :: Xe_pointer(:)     !Pointer of Xe_tally = loc

    real(8), allocatable :: Xe_prod(:)        !Xe_prod(loc) = cumulative xenon production
    real(8), allocatable :: Xe_trans(:)       !Xe_trans(loc+1) = cumulative xenon transmutation
    
    
    real(8), allocatable :: Xe_numden(:)      !Xenon number densities [#/cc*1.d-24]
    real(8), allocatable :: avg_Xe_numden(:)  !Averaged number densities of Xenon [#/cc*1.d-24]
    integer, allocatable :: Xe135_mt_iso(:)   !Index for Xenon135 in local list
    real(8), allocatable :: Xe_cum_yield(:)   !Cumulative xenon yield : summation over independent yields of Sn-135, Sb-135, Te-135, I-135, Xe-135, Xe-135m

    ! predictor corrector
    integer:: preco = 0     ! predictor-corrector (PC)
    integer:: porc = 0      ! predictor or corrector ?
    real(8):: stepsize0
    

    contains 
    
! ===================================================================
!     getdepletionlibrary :: Read depletion library and store data
! ===================================================================
subroutine getdepletionlibrary
    implicit none
    logical :: file_exists
    character(80)  :: line0, line1, line2
    integer, parameter :: rd_decay = 20160818
    integer, parameter :: rd_1gcx  = 20160819
    integer, parameter :: rd_yield = 20160820
    integer :: nuclid
    integer :: anum, mnum, inum, nnum !atomic number, mass number, isomer state, neutron number of nuclide
    integer :: nlb
    real(8) :: thalf, yyn
    integer :: i, j, k
    integer :: nsum = 1
    logical :: yes

    !==============================================================================
    !Burnup result files
!    if(icore==score) then
!        !open(prt_ntpy, file="dep_monitor",action="write",status="replace") 
!        !open(prt_dep,file="depletion1.out",action="write",status="replace")
!
!        ! burnup materials output file
!        if ( istep_burnup == 0 ) then
!        bufile = 'bumat1.out'
!        do
!        inquire(file=trim(bufile),exist=yes)
!        if ( yes ) then
!            nsum = nsum + 1
!            if ( nsum < 10 ) then
!                write(bufile,'(a,i1,a)'), 'bumat',nsum,'.out'
!            else
!                write(bufile,'(a,i2,a)'), 'bumat',nsum,'.out'
!            end if
!    
!        else
!            exit
!        end if
!        end do
!        ! open a new file
!        open(prt_bumat,file=trim(bufile))
!        close(prt_bumat)
!        end if
!    end if

    !==============================================================================
    !Initialization
    do inum = 0, 1
    do nnum = 0, 160
    do anum = 1, 100
      nuclide(inum,nnum,anum)%data_exist = .false.
      nuclide(inum,nnum,anum)%fy_idx = 0 
      nuclide(inum,nnum,anum)%idx = 0 
    end do
    end do
    end do


    !==============================================================================
    !Read decay library
    open(rd_decay, file=trim(dep_lib)//"decay_library",action="read")
    10 continue
    read(rd_decay,'(A80)',end=100) line0 !Title
    do
      read(rd_decay,'(A80)',end=100) line1 !First  line
      read(line1(1:4),'(i)') nlb; if(nlb==-1) go to 10 !library index (1=activation products, 2=actinides, 3=fission products)
      read(rd_decay,'(A80)',end=100) line2 !Second line

      read(line1(7:12),'(i)') nuclid     !six-digit nuclide identifier
      anum = nuclid/10000
      mnum = (nuclid - anum*10000)/10
      nnum = mnum - anum
      inum = nuclid - anum*10000 - mnum*10

      if(nuclide(inum,nnum,anum)%data_exist==.true.) cycle !Skip (duplicated case)

      nuclide(inum,nnum,anum)%data_exist = .true.
      nuclide(inum,nnum,anum)%sng = 0.d0
      nuclide(inum,nnum,anum)%sn2n = 0.d0
      nuclide(inum,nnum,anum)%sna = 0.d0
      nuclide(inum,nnum,anum)%snp = 0.d0
      nuclide(inum,nnum,anum)%sn3n = 0.d0
      nuclide(inum,nnum,anum)%snf = 0.d0
      nuclide(inum,nnum,anum)%sngx = 0.d0
      nuclide(inum,nnum,anum)%sn2nx = 0.d0
      nuclide(inum,nnum,anum)%fy_idx = 0

      read(line1(15:15),'(i)') nuclide(inum,nnum,anum)%iu  !Time unit designation
      read(line1(21:29),'(f)') thalf !half-life
      select case (nuclide(inum,nnum,anum)%iu)
        case(1) !second
          nuclide(inum,nnum,anum)%lambda = 6.931471806d-1/thalf
        case(2) !minutes
          nuclide(inum,nnum,anum)%lambda = 1.155245301d-2/thalf
        case(3) !hours
          nuclide(inum,nnum,anum)%lambda = 1.925408835d-4/thalf 
        case(4) !days
          nuclide(inum,nnum,anum)%lambda = 8.022536812d-6/thalf 
        case(5) !years
          nuclide(inum,nnum,anum)%lambda = 2.197955291d-8/thalf 
        case(6) !stable
          nuclide(inum,nnum,anum)%lambda = 0.d0
        case(7) !10^3 years
          nuclide(inum,nnum,anum)%lambda = 2.197955291d-11/thalf 
        case(8) !10^6 years
          nuclide(inum,nnum,anum)%lambda = 2.197955291d-14/thalf 
        case(9) !10^9 years
          nuclide(inum,nnum,anum)%lambda = 2.197955291d-17/thalf 
      end select
      read(line1(31:39),'(f)') nuclide(inum,nnum,anum)%fbx   !fraction of beta- decay to excited state
      read(line1(41:49),'(f)') nuclide(inum,nnum,anum)%fpec  !fraction of all beta+ decay and electron capture 
      read(line1(51:59),'(f)') nuclide(inum,nnum,anum)%fpecx !fraction of beta+ decay and electron capture to excited state
      read(line1(61:69),'(f)') nuclide(inum,nnum,anum)%fa  !fraction of all alpha decay
      read(line1(71:79),'(f)') nuclide(inum,nnum,anum)%fit   !fraction of gamma-ray release from excited state to ground state

      read(line2(21:29),'(f)') nuclide(inum,nnum,anum)%fsf   !fraction of decay events of spontaneous fission
      read(line2(31:39),'(f)') nuclide(inum,nnum,anum)%fn  !fraction of decay events that are (beta- & neutron) decays 
      read(line2(41:49),'(f)') nuclide(inum,nnum,anum)%qrec  !Average, total recoverable energy [MeV] for each decay event
      read(line2(51:59),'(f)') nuclide(inum,nnum,anum)%abund !Naturally occurring isotopic abundance in atom percent
      read(line2(61:69),'(f)') nuclide(inum,nnum,anum)%arcg  !Radioactivity concentration guide for continuous inhalation of nuclide
      read(line2(71:79),'(f)') nuclide(inum,nnum,anum)%wrcg  !Radioactivity concentration guide for continuous ingestion of nuclide
      nuclide(inum,nnum,anum)%fb = 1.d0 - nuclide(inum,nnum,anum)%fbx - nuclide(inum,nnum,anum)%fpec - nuclide(inum,nnum,anum)%fpecx - nuclide(inum,nnum,anum)%fa - nuclide(inum,nnum,anum)%fit - nuclide(inum,nnum,anum)%fsf - nuclide(inum,nnum,anum)%fn

      if(nuclide(inum,nnum,anum)%fb<0.d0) nuclide(inum,nnum,anum)%fb = 0.d0 !sum of decay fractions is greater than 1 

    end do
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
    !Read fission yield library
    open(rd_yield, file=trim(dep_lib)//"yield_library", action="read")
    read(rd_yield,*) nfssn
    allocate(fssn_zai(1:nfssn))
    read(rd_yield,*) fssn_zai(:)
      do i=1, nfssn
        anum = fssn_zai(i)/10000
        mnum = (fssn_zai(i) - anum*10000)/10
        nnum = mnum - anum
        inum = fssn_zai(i) - anum*10000 - mnum*10
        nuclide(inum,nnum,anum)%fy_idx = i  !fission yield data index in yield data array
      end do
    read(rd_yield,*) nfp
    allocate(fp_zai(1:nfp))
    allocate(yield_data(1:nfp,1:nfssn))
    allocate(Xe_cum_yield(1:nfssn))
    !Cumulative xenon yield : summation over independent yields of Sn-135, Sb-135, Te-135, I-135, Xe-135, Xe-135m
    Xe_cum_yield(1:nfssn) = (/5.28E+00, 6.09E+00, 6.26E+00, 6.13E+00, 7.34E+00, 6.76E+00, 6.65E+00, 6.65E+00/) 
    do i=1, nfp
      read(rd_yield,*) fp_zai(i), (yield_data(i,j),j=1,nfssn) !yield data (fission product index, fission nuclide index)
    end do
    close(rd_yield)


    !==============================================================================
    !Sort nuclide data
    nnuc = 0
    do anum=1,100
    do nnum=0,160
    do inum=0,1
      if(nuclide(inum,nnum,anum)%data_exist==.true.) then
      nnuc = nnuc + 1
      nuclide(inum,nnum,anum)%idx = nnuc
      zai_idx(nnuc) = anum*10000+(nnum+anum)*10+inum
      end if
    end do
    end do
    end do
    maxnnz = int(0.1 * nnuc**2) ! Assuming less than 10% sparsity  
    allocate(Acsr(1:maxnnz))
    allocate(jAcsr(1:maxnnz))
    allocate(iAcsr(1:nnuc+1))
    
end subroutine getdepletionlibrary



! ===================================================================
!         Initialize burnup tallies 
! ===================================================================
subroutine INIT_BURNUP
    integer :: i, imat
    
    if (istep_burnup == NSTEP_BURNUP) do_burn = .false. 
    if (.not. do_burn) return

    if (istep_burnup > 0) do_Xe_search = (Xe_search .and. bstep_size/86400.d0 > 5.00001d0)
    cram_init = .false.
    do imat = 1, n_materials
        materials(imat)%flux = 0.0d0 
        !materials(imat)%pwr = 0.0d0 
        materials(imat)%ogxs(:,:) = 0.0d0
        materials(imat)%kappa = 0.0D0
        !allocate(materials(imat)%full_numden(nnuc))
    enddo 
    bstep_size = burn_step(istep_burnup+1) - burn_step(istep_burnup) 

    ! initialization
    avg_power = 0D0
    
    !> Xe 
    if (.not. do_Xe_search) return 
    Xe_prod(:) = 0.0d0
    Xe_trans(:) = 0.0d0
    Xe_numden(:) = 0.0d0
    do imat = 1, n_materials 
        Xe_local: do i = 1, materials(imat)%n_iso 
            if (materials(imat)%ace_idx(i) == Xe135_iso) then
                Xe135_mt_iso(imat) = i 
                exit Xe_local
            endif 
        enddo Xe_local
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
        allocate(materials(imat)%ogxs(1:materials(imat)%n_iso,1:4)) 
        materials(imat)%ogxs(:,:) = 0.0d0
        materials(imat)%kappa = 0D0
        
        if (materials(imat)%depletable) then
            allocate(materials(imat)%full_numden(1:nnuc))
            if ( preco == 1 ) allocate(materials(imat)%full_numden1(1:nnuc))
            materials(imat)%full_numden = 0.d0     !Initialize number density
        endif
        
        do mt_iso = 1, materials(imat)%n_iso
            iso = materials(imat)%ace_idx(mt_iso)
            
            anum = ace(iso)%zaid/1000
            mnum = ace(iso)%zaid - anum*1000
            nnum = mnum - anum
            if(materials(imat)%depletable==.true. .and. do_burn) then 
                materials(imat)%full_numden(nuclide(0,nnum,anum)%idx) = materials(imat)%numden(mt_iso)
            endif 
        end do
        
    enddo 
    
    istep_burnup = 0 
    
    !> Xe 
    if (.not. Xe_search) return 
    
    !> Find Xe 135 index
    Xe135_iso = 0 
    Xe_search: do i = 1, num_iso 
        if (ace(i)%zaid == 54135) then     
            Xe135_iso = i 
            exit Xe_search
        endif 
    enddo Xe_search
    if (Xe135_iso == 0 .and. icore==score ) then 
        print *, "setmat() - WARNING :: NO Xe-135 LIBRARY IN THE INVENTORY"
    endif 
    
    allocate(Xe_prod(1:n_materials))
    allocate(Xe_trans(1:n_materials))
    allocate(Xe_numden(1:n_materials))
    allocate(Xe135_mt_iso(1:n_materials))
    
    Xe135_iso = find_ACE_iso_idx_zaid(541350) 
    do imat = 1, n_materials            
        do mt_iso = 1, materials(imat)%n_iso
            iso = materials(imat)%ace_idx(mt_iso)
            if(Xe135_iso == iso) Xe135_mt_iso(imat) = mt_iso
        end do
    enddo 
    
end subroutine setmat


! ===================================================================
!         Inline Xe-135  
! ===================================================================
subroutine inline_xenon
    integer :: imat
    integer :: loc1, loc2
    real(8) :: inv_nactive, temp, rcvbufXe(n_materials)
    type (Material_CE), pointer :: mat
    
    if (.not. do_Xe_search) return 
    
    call MPI_ALLREDUCE(Xe_prod,rcvbufXe,n_materials,MPI_REAL8,MPI_SUM,core,ierr)    
    call MPI_ALLREDUCE(Xe_trans,rcvbufXe,n_materials,MPI_REAL8,MPI_SUM,core,ierr)    
    
    inv_nactive = 1.0d0 / dble(n_act) 
    do imat = 1, n_materials
        mat => materials(imat)

        if(mat%fissionable == .false.) cycle

        !Obtain volume averaged Xe_tally
        Xe_prod(imat)   = Xe_prod(imat)/materials(imat)%vol   !Xenon production
        Xe_trans(imat) = Xe_trans(imat)/materials(imat)%vol !Xenon transmutation
    
        temp = Xe_prod(imat)/(Xe_trans(imat) + nuclide(0,81,54)%lambda)*1.d-24

        !Set Xenon equilibrium density
        materials(imat)%numden(Xe135_mt_iso(imat)) = temp
        if(curr_cyc > n_inact) then
            !Accumulate Xenon number density for active cycles
            Xe_numden(imat) = Xe_numden(imat) + temp
            !Update Xenon number density for depletion calculation
            if(curr_cyc == n_totcyc) then
                Xe_numden(imat) = Xe_numden(imat)*inv_nactive
                materials(imat)%numden(Xe135_mt_iso(imat)) =  Xe_numden(imat)
            end if
        end if
    end do
    
    !if(icore==score) print *, "inline xenon equilibrium updated"
    
end subroutine inline_xenon
    

! ===================================================================
!         Tally burnup parameters 
! ===================================================================
subroutine TALLY_BURNUP (imat, distance, wgt, erg, macro_xs)
    integer, intent(in) :: imat  ! material index (p%material) 
    real(8), intent(in) :: distance, wgt, erg, macro_xs
    integer :: iso, mt_iso
    integer :: i, ierg
    integer :: fy_midx, diff, tmp, Xe_ptr

    type (Material_CE), pointer :: mat
    real(8) :: micro_f, micro_d, micro_a, ipfac, micro_q
    integer :: anum, mnum, inum, nnum ! atomic number, mass number
                                      ! isomer state, neutron number of nuclide
    
    if(E_mode==0) return        !Multigroup  -> return
    if(do_burn==.false.) return !No burnup calculation -> return
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
            micro_d = ace(iso)%sigd(ierg) + ipfac*(ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
            micro_f = 0.d0
            ! Fissionable Material
            if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
                micro_f = ace(iso)%sigf(ierg) + ipfac*(ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
            endif
            micro_a = micro_d + micro_f
            
            anum = ace(iso)%zaid/1000
            mnum = (ace(iso)%zaid - anum*1000)
            nnum = mnum - anum
            inum = 0
            if(nuclide(inum,nnum,anum)%fy_idx>0) then
                fy_midx = nuclide(inum,nnum,anum)%fy_idx
                Xe_prod(imat) = Xe_prod(imat) + wgt*distance*mat%numden(mt_iso)*micro_f &
                    *Xe_cum_yield(fy_midx)*0.01d0
            else
                diff = 999999
                do i=1,nfssn
                    tmp = abs(fssn_zai(i) - zai_idx(nuclide(inum,nnum,anum)%idx))
                    if(diff > tmp) then
                    diff = tmp
                    fy_midx = i
                    end if
                end do
                Xe_prod(imat) = Xe_prod(imat) + wgt*distance*mat%numden(mt_iso)*micro_f &
                    *Xe_cum_yield(fy_midx)*0.01d0
            end if
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

    !%omp atomic
    mat%kappa = mat%kappa + wgt*distance*macro_xs
    
    do mt_iso=1, mat%n_iso
        iso = mat%ace_idx(mt_iso)
        call getierg(iso,ierg,erg)
        ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg)) &
            /(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
        micro_d   = ace(iso)%sigd(ierg) + ipfac &
             *(ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
        !$omp atomic
        mat%ogxs(mt_iso,1) = mat%ogxs(mt_iso,1) &
            + wgt*distance*micro_d*1.d-24 !(n,g)
        !$omp atomic
        mat%ogxs(mt_iso,2) = mat%ogxs(mt_iso,2) &
            + wgt*distance*getxs(16,iso,erg,ierg)*1.d-24 !(n,2n)
    
        if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then  !Fissionable
            micro_f   = ace(iso)%sigf(ierg) + ipfac &
                *(ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
            !print *,ace(iso)%library, erg, micro_f
            !$omp atomic
            mat%ogxs(mt_iso,3) = mat%ogxs(mt_iso,3) &
                + wgt*distance*getxs(17,iso,erg,ierg)*1.d-24 !(n,3n)
            !$omp atomic
            mat%ogxs(mt_iso,4) = mat%ogxs(mt_iso,4) &
                + wgt*distance*micro_f*1.d-24 !(n,f)
            
            !!$omp atomic
            !mat%pwr = mat%pwr + wgt*distance*micro_f*201*mat%numden(mt_iso)* barn
            
         else !Non-fissionable
            !$omp atomic
            mat%ogxs(mt_iso,3) = mat%ogxs(mt_iso,3) &
                + wgt*distance*getxs(107,iso,erg,ierg)*1.d-24 !(n,alpha)
            !$omp atomic
            mat%ogxs(mt_iso,4) = mat%ogxs(mt_iso,4) &
                + wgt*distance*getxs(103,iso,erg,ierg)*1.d-24 !(n,p)
         endif
     enddo
    
end subroutine tally_burnup


! ===================================================================
!         MPI reduce burnup tallies 
! ===================================================================
subroutine MPI_reduce_burnup 
    type(Material_CE), pointer:: mat
    integer :: iso, imat, i 
    real(8) :: rcvbuf0, rcvbuf, rcvbufarr(4)
    real(8), allocatable :: sndbufarrlong(:), rcvbufarrlong(:)
    
    if (.not. do_burn) return 

    do imat = 1, n_materials
        if (.not. materials(imat)%depletable) cycle 

        ! pointer
        mat => materials(imat)

        ! MPI
        call MPI_ALLREDUCE(mat%kappa,rcvbuf0,1,MPI_REAL8,MPI_SUM,core,ierr)
        call MPI_ALLREDUCE(mat%flux,rcvbuf,1,MPI_REAL8,MPI_SUM,core,ierr)

        if ( iscore ) mat%kappa = rcvbuf0/rcvbuf
        mat%flux = rcvbuf / (dble(n_act) * mat%vol)
        
        allocate(rcvbufarrlong(1:mat%n_iso*4)) 
        allocate(sndbufarrlong(1:mat%n_iso*4)) 
        
        do iso = 1, mat%n_iso
            sndbufarrlong((iso-1)*4+1:(iso-1)*4+4) = mat%ogxs(iso,:)
        enddo 
        call MPI_ALLREDUCE(sndbufarrlong,rcvbufarrlong,mat%n_iso*4, &
                            MPI_REAL8,MPI_SUM,core,ierr)            
                            
        rcvbufarrlong(:) = rcvbufarrlong(:) / rcvbuf
        !rcvbufarrlong(:) = rcvbufarrlong(:) / val
        
        do iso = 1, mat%n_iso
            mat%ogxs(iso,:) = rcvbufarrlong((iso-1)*4+1:(iso-1)*4+4) 
        enddo
        deallocate(rcvbufarrlong)
        deallocate(sndbufarrlong) 
    enddo


end subroutine MPI_reduce_burnup
    

! ===================================================================
!         Depletion :: Make depletion matrix and solve 
! ===================================================================
subroutine DEPLETION
use FMFD_HEADER, only: fmfdon, nrings, nfm
implicit none

integer :: imat, jmem, jnuc, knuc
integer :: mt_iso, iso, niso
integer :: anum, mnum, nnum, inum
integer :: anum1, mnum1, nnum1, inum1
integer :: i, j, k, ii, jj
logical :: exist_in_MC
real(8) :: real_flux
real(8) :: id(1:4)
integer :: fy_midx, diff, tmp, ierr
real(8) :: nxt_full_numden(1:nnuc)
type (Material_CE), pointer :: mat
integer :: iso_idx(nnuc)
real(8) :: weight
real(8) :: tt0, tt1, tt2, tt3


!call omp_set_num_threads(1)

if(do_burn==.false.) return

!if ( iscore ) print*, "power2", avg_power
avg_power = avg_power / dble(n_act)

call MPI_BCAST(avg_power, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)


! Initialize material independent burnup matrix
if( istep_burnup==0 .and. .not. allocated(bMat) ) then
allocate(bMat(1:nnuc,1:nnuc))   !2-D burnup matrix : row to column transition
allocate(bMat0(1:nnuc,1:nnuc)) !2-D burnup matrix : row to column transition (material independent)
bMat0 = 0.d0


!Build material independent burnup matrix
do jnuc = 1, nnuc
    anum = zai_idx(jnuc)/10000
    mnum = (zai_idx(jnuc) - anum*10000)/10
    nnum = mnum - anum
    inum = zai_idx(jnuc) - anum*10000 - mnum*10
    
    !if (zai_idx(jnuc) == 922350 ) jmem = jnuc 
    !Stable (iu=6) -> no radioactive decay
    if(nuclide(inum,nnum,anum)%iu==6) cycle

    !Radioactive decay : (I,N,Z) -> something else
    bMat0(jnuc,jnuc) = bMat0(jnuc,jnuc) - nuclide(inum,nnum,anum)%lambda

    !Negatron emission and goto ground state : (I,N,Z) -> (0,N-1,Z+1)
    if(nuclide(inum,nnum,anum)%fb > 0.d0) then
        knuc = nuclide(0,nnum-1,anum+1)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = & 
            bMat0(knuc,jnuc) + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fb
    end if

    !Negatron emission and goto excited state : (I,N,Z) -> (1,N-1,Z+1)
    if(nuclide(inum,nnum,anum)%fbx > 0.d0) then
        knuc = nuclide(1,nnum-1,anum+1)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = & 
            bMat0(knuc,jnuc) + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fbx
    end if

    !Positron emission and goto ground state : (I,N,Z) -> (0,N+1,Z-1)
    if(nuclide(inum,nnum,anum)%fpec > 0.d0) then
        knuc = nuclide(0,nnum+1,anum-1)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = & 
            bMat0(knuc,jnuc) + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fpec
    end if

    !Positron emission and goto excited state : (I,N,Z) -> (1,N+1,Z-1)
    if(nuclide(inum,nnum,anum)%fpecx > 0.d0) then
        knuc = nuclide(1,nnum+1,anum-1)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = &
            bMat0(knuc,jnuc) + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fpecx
    end if
    
    !Alpha particle emission and ground state : (I,N,Z) -> (0,N-2,Z-2) + (0,2,2)
    if(nuclide(inum,nnum,anum)%fa > 0.d0) then
    knuc = nuclide(0,nnum-2,anum-2)%idx
    if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) &
                + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fa
        !** consider alpha particle
        knuc = nuclide(0,2,2)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) & 
                + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fa
    end if

    !Isomeric transition from excited state to ground state : (1,N,Z) -> (0,N,Z)
    if(nuclide(inum,nnum,anum)%fit > 0.d0) then
        knuc = nuclide(0,nnum,anum)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) &
            + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fit
    end if
    
    !Negatron emssion and neutron emission : (I,N,Z) -> (0,N-2,Z+1)
    if(nuclide(inum,nnum,anum)%fn > 0.d0) then
        knuc = nuclide(0,nnum-2,anum+1)%idx
        if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) &
            + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fn
    end if

    !Spontaneous fission
    if(nuclide(inum,nnum,anum)%fsf > 0.d0) then
    if(nuclide(inum,nnum,anum)%fy_idx>0) then
        do j=1,nfp
            anum1 = fp_zai(j)/10000
            mnum1 = (fp_zai(j) - anum1*10000)/10
            nnum1 = mnum1 - anum1
            inum1 = fp_zai(j) - anum1*10000 - mnum1*10
            knuc = nuclide(inum1,nnum1,anum1)%idx
            if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) &
                + nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fsf &
                *yield_data(j,nuclide(inum,nnum,anum)%fy_idx)*0.01d0
        end do
    else
        if(anum>=89) then !For actinides
        !Find nearest isotopes which has fission yield data
        diff = 999999
        do i=1,nfssn
            tmp=abs(fssn_zai(i) - zai_idx(nuclide(inum,nnum,anum)%idx)) 
            if(diff > tmp) then 
                diff = tmp
                fy_midx = i
            end if
        end do
        do j=1,nfp
            anum1 = fp_zai(j)/10000
            mnum1 = (fp_zai(j) - anum1*10000)/10
            nnum1 = mnum1 - anum1
            inum1 = fp_zai(j) - anum1*10000 - mnum1*10
            knuc = nuclide(inum1,nnum1,anum1)%idx
            if(knuc/=0) bMat0(knuc,jnuc) = bMat0(knuc,jnuc) + &
                nuclide(inum,nnum,anum)%lambda*nuclide(inum,nnum,anum)%fsf &
                *yield_data(j,fy_midx)*0.01d0
        end do
        end if
    end if
    end if
end do
end if  ! initialization

!> inventory mass
if ( iscore ) then
tmass = 0
fmass = 0
end if

!if ( iscore .and. porc /= 2 ) then
!    open(prt_bumat,file=trim(bufile),access='append',status='old')
!    write(prt_bumat,1), '   =========================================='
!    !write(prt_dep,1),   '   =========================================='
!    write(prt_bumat,2), '      Burnup step', istep_burnup+1
!    !write(prt_dep,2),   '      Burnup step', istep_burnup+1
!    if ( preco == 1 ) then
!    if ( porc == 1 ) then
!    write(prt_bumat,3), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS (predictor)'
!    !write(prt_dep,3),   burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS (predictor)'
!    elseif ( porc == 2 ) then
!    write(prt_bumat,3), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS (corrector)'
!    !write(prt_dep,3),   burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS (corrector)'
!    end if
!    else
!    write(prt_bumat,3), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
!    !write(prt_dep,3),   burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
!    end if
!    write(prt_bumat,1), '   =========================================='
!    !write(prt_dep,1),   '   =========================================='
!    write(prt_bumat,*)
!    !write(prt_dep,*)
!    1 format(a45)
!    2 format(a17,i4)
!    3 format(f14.2,a16)
!end if

!Normalization constant to be real power
ULnorm = RealPower/(avg_power*eVtoJoule)

! -------------------------------------------------------------------------
! sharing flux
if ( DTMCBU ) &
call MPI_BCAST(materials(:)%flux,n_materials,MPI_REAL8,score,MPI_COMM_WORLD,ierr)


!Substitute burnup matrix element
do imat = mp1(icore), mp2(icore)
    !tt0 = MPI_WTIME()
    !material imat is not burned
    if(.not. materials(imat)%depletable) cycle

    mat => materials(imat)
!    if ( porc /= 2 ) then 
!        write(prt_bumat, *) '' 
!        write(prt_bumat, *) mat%mat_name 
!        write(prt_bumat, *) '' 
!        !write(prt_dep,*)
!        !write(prt_dep,*) mat%mat_name
!    endif         
    !Call the material independent burnup matrix 
    bMat = bMat0*bstep_size

    
    !Calculate real flux (volume-averaged)
    if ( DTMCBU ) then
        real_flux = mat%flux
        !print*, "flux", imat, real_flux
!        nnum = size(mat%dtmc) / 4
!        do inum = 1, nnum
!            id(1:4) = mat%dtmc(4*(inum-1)+1:4*(inum-1)+4)
!            real_flux = real_flux + buflux(id(1),id(2),id(3),id(4))
!        end do
!        print*, real_flux, buflux(id(1),id(2),id(3),id(4))
!        real_flux = real_flux/dble(mat%vol*dble(nnum))
!        real_flux = ULnorm*mat%flux
    else
        real_flux = ULnorm*mat%flux
    end if

    ! total mass calculation
    call TOTAL_MASS_CALC(mat)

    ! ======================================================================================
    if ( porc /= 2 ) then
!        if ( fmfdon ) then
!        write(prt_dep,*), " - POWER (MW)      : ", bukapa(id(1),id(2),id(3),id(4))
!        else
!        write(prt_dep,*), " - POWER (MW)      : ", real_flux*mat%kappa*mat%vol*eVtoJoule
!        endif
        !write(prt_dep,*), " - FLUX            : ", real_flux
        !write(prt_dep,*), " - NUMBER DENSITY : ", sum(mat%numden(:))
        !write(prt_dep,*), " - CROSS SECTION   : "
!         do mt_iso = 1, mat%n_iso
!         write(prt_dep,5), ace(mat%ace_idx(mt_iso))%library, (mat%ogxs(mt_iso,jj), jj = 1, 4)
!         end do
!        write(prt_dep,6), (sum(mat%ogxs(:,jj)*mat%numden(:)), jj = 1, 4)
!        write(prt_dep,*)
        5 format(a13,4es20.12)
        6 format(4es24.16)
    endif 

    ! ---------------------------------------------------------------------
    ! predictor-corrector
    ! CE/LI
    if ( preco == 1 ) then
    ! save (predictor)
    if ( porc == 1 ) then
        ! flux
        mat%flux1 = real_flux
        ! cross section
        if ( allocated(mat%ogxs1) ) deallocate(mat%ogxs1)
        if ( allocated(mat%ace_idx1) ) deallocate(mat%ace_idx1)
        allocate(mat%ogxs1(mat%n_iso,1:4))
        allocate(mat%ace_idx1(mat%n_iso))
        mat%ogxs1(:,:) = mat%ogxs(:,:)
        mat%ace_idx1(:) = mat%ace_idx(:)
        ! number density
        mat%full_numden1 = mat%full_numden
    ! load (corrector)
    else if ( porc == 2 ) then
        ! flux
        real_flux = (real_flux+mat%flux1)*5D-1
        ! cross section
        jj = 1
        do ii = 1, mat%n_iso
        if ( mat%ace_idx(ii) == mat%ace_idx1(jj) ) then
            mat%ogxs(ii,:) = (mat%ogxs(ii,:)+mat%ogxs1(jj,:))*5D-1
            jj = jj + 1
        else
            mat%ogxs(ii,:) = mat%ogxs(ii,:)*5D-1
        end if
        end do
        ! number density
        mat%full_numden = mat%full_numden1
    end if
    ! LE
    else if ( preco == 2 ) then
        ! --- data save ---
        ! flux
        mat%flux1 = real_flux
        ! cross section
        if ( allocated(mat%ogxs1) ) deallocate(mat%ogxs1)
        if ( allocated(mat%ace_idx1) ) deallocate(mat%ace_idx1)
        allocate(mat%ogxs1(mat%n_iso,1:4))
        allocate(mat%ace_idx1(mat%n_iso))
        mat%ogxs1 = mat%ogxs
        mat%ace_idx1 = mat%ace_idx

        ! --- data load ---
        if ( istep_burnup > 0 ) then
        weight = bstep_size/(stepsize0*2.D0)
        ! flux
        mat%flux = -weight*mat%flux0+(1D0+weight)*mat%flux1
        ! cross section
        if ( size(mat%ogxs0) == size(mat%ogxs) ) then
            mat%ogxs = -weight*mat%ogxs0+(1D0+weight)*mat%ogxs1
        else if ( size(mat%ogxs0) < size(mat%ogxs1) ) then
            jj = 1
            print*, "ace", allocated(mat%ace_idx0)
            do ii = 1, mat%n_iso
            if ( mat%ace_idx1(ii) == mat%ace_idx0(jj) ) then
                mat%ogxs(ii,:) = &
                    -weight*mat%ogxs0(ii,:)+(1+weight)*mat%ogxs1(jj,:)
                jj = jj + 1
            else
                mat%ogxs(ii,:) = mat%ogxs1(ii,:)*(1+weight)
            end if
            end do
        end if
        end if
    end if

    
    !Build burnup matrix with cross section obtained from MC calculation
    DO_ISO: do mt_iso = 1, mat%n_iso
    
        iso = mat%ace_idx(mt_iso)
        anum = ace(iso)%zaid/1000
        mnum = (ace(iso)%zaid - anum*1000)
        nnum = mnum - anum
        inum = 0
        jnuc = nuclide(inum,nnum,anum)%idx
        
        !(n,g)     leading to ground state : (I,N,Z) -> (0,N+1,Z)
        knuc = nuclide(0,nnum+1,anum)%idx
        !production of knuc
        if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                    + real_flux*mat%ogxs(mt_iso,1)*bstep_size
        !reduction of junc
        bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,1)*bstep_size
        
        if (nnum < 1) goto 20  ! Hydrogen
        
        !(n,2n)    leading to ground state : (I,N,Z) -> (0,N-1,Z)
        knuc = nuclide(0,nnum-1,anum)%idx
        !production of knuc
        if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                    + real_flux*mat%ogxs(mt_iso,2)*bstep_size
        !reduction of jnuc
        bMat(jnuc,jnuc) = bMat(jnuc,jnuc) &
                        - real_flux*mat%ogxs(mt_iso,2)*bstep_size
    
        20 continue
        
        if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
            !(n,3n)    leading to ground state : (I,N,Z) -> (0,N-2,Z) : fissionable
            knuc = nuclide(0,nnum-2,anum)%idx
            !production of knuc
            if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                        + real_flux*mat%ogxs(mt_iso,3)*bstep_size
            !reduction of jnuc
            bMat(jnuc,jnuc) = bMat(jnuc,jnuc) &
                            - real_flux*mat%ogxs(mt_iso,3)*bstep_size
    
            !(n,f) leading to ground state 
            ! : (I,N,Z) -> (i1,N1,Z1) + (i2,N2,Z2) : fissionable
            if(nuclide(inum,nnum,anum)%fy_idx>0) then
                do j=1,nfp
                    anum1 = fp_zai(j)/10000
                    mnum1 = (fp_zai(j) - anum1*10000)/10
                    nnum1 = mnum1 - anum1
                    inum1 = fp_zai(j) - anum1*10000 - mnum1*10
                    knuc = nuclide(inum1,nnum1,anum1)%idx
                    if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                        + real_flux*mat%ogxs(mt_iso,4) &
                        *yield_data(j,nuclide(inum,nnum,anum)%fy_idx)*0.01d0*bstep_size
                end do
            else
                if(anum>=89) then !For actinides
                    !Find nearest isotopes which has fission yield data
                    diff = 999999
                    do i=1,nfssn
                        tmp=abs(fssn_zai(i) - zai_idx(nuclide(inum,nnum,anum)%idx)) 
                        if(diff > tmp) then 
                        diff = tmp
                        fy_midx = i
                        end if
                    end do
                    do j=1,nfp
                        anum1 = fp_zai(j)/10000
                        mnum1 = (fp_zai(j) - anum1*10000)/10
                        nnum1 = mnum1 - anum1
                        inum1 = fp_zai(j) - anum1*10000 - mnum1*10
                        knuc = nuclide(inum1,nnum1,anum1)%idx
                        if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                            + real_flux*mat%ogxs(mt_iso,4)*yield_data(j,fy_midx) &
                            *0.01d0*bstep_size
                    end do
                end if
            end if
            !reduction of jnuc
            bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,4)*bstep_size
        else
            !(n,alpha) leading to ground state : (I,N,Z) -> (0,N-1,Z-2) + (0,2,2) : non-fission
            if (anum < 3) goto 10
            knuc = nuclide(0,nnum-1,anum-2)%idx
            ! production of knuc
            if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                        + real_flux*mat%ogxs(mt_iso,3)*bstep_size
            !reduction of jnuc
            bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,3)*bstep_size
            
            !** consider alpha production
            knuc = nuclide(0,2,2)%idx
            ! production of alpha particle
            if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                        + real_flux*mat%ogxs(mt_iso,3)*bstep_size
        
            10 continue
            if (anum < 2) cycle DO_ISO
            !(n,p)     leading to ground state : (I,N,Z) -> (0,N,Z-1) + (0,0,1) : non-fission
            knuc = nuclide(0,nnum,anum-1)%idx
            ! productino of knuc
            if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                        + real_flux*mat%ogxs(mt_iso,4)*bstep_size
            ! reduction of jnuc
            bMat(jnuc,jnuc) = bMat(jnuc,jnuc) - real_flux*mat%ogxs(mt_iso,4)*bstep_size
            !** consider proton production
            knuc = nuclide(0,0,1)%idx
            !production of proton
            if(knuc/=0) bMat(knuc,jnuc) = bMat(knuc,jnuc) &
                                        + real_flux*mat%ogxs(mt_iso,4)*bstep_size
        end if
    end do DO_ISO
    
    
    !Solve the burnup matrix equation 
    !tt2 = MPI_WTIME()
    select case(matrix_exponential_solver)
    case(0)
        call CRAM(bMat, materials(imat)%full_numden(:), nxt_full_numden(:)) 
        do jnuc=1, nnuc
            if(nxt_full_numden(jnuc)<0.d0) nxt_full_numden(jnuc) = 0.d0
        end do
        
    case default 
        print *, "ERROR :: No such matrix_exponential_solver option", matrix_exponential_solver
        stop 
    end select
    !tt3 = MPI_WTIME()
    !print*, "2", tt3-tt2
        
        
    !Update number density
    mat%full_numden = nxt_full_numden

    
    ! ======================================================================================
    ! Reset material isotope inventory
    knuc = 0 
    do jnuc=1, nnuc 
        tmp = find_ACE_iso_idx_zaid(zai_idx(jnuc))
        if(mat%full_numden(jnuc)>0.d0 .and. tmp > 0) then 
            knuc = knuc + 1
            iso_idx(knuc) = jnuc
        endif
    end do

    mat%n_iso = knuc
    deallocate(mat%ace_idx); allocate(mat%ace_idx(1:knuc)); mat%ace_idx(1:knuc) = 0
    deallocate(mat%numden); allocate(mat%numden(1:knuc)); mat%numden(1:knuc) = 0
    deallocate(mat%ogxs); allocate(mat%ogxs(1:knuc,1:4)); mat%ogxs(1:knuc,1:4) = 0
    
    i = 0
    do mt_iso=1, mat%n_iso
        ! find ace_idx
        tmp = find_ACE_iso_idx_zaid(zai_idx(iso_idx(mt_iso)))
        if ( tmp /= 0 ) then 
            i = i + 1
            mat%ace_idx(mt_iso) = tmp
            mat%numden(mt_iso)  = mat%full_numden(iso_idx(mt_iso))
            !write(prt_bumat,'(a15,e14.5)') ace(mat%ace_idx(mt_iso))%library, mat%numden(mt_iso)
        endif
    end do
    
!    write(prt_bumat, *) 'Num isotope', i
!    write(prt_bumat, *) '' 
!    close(prt_bumat)

    !tt1 = MPI_WTIME()
    !print*, "1", tt1-tt0
end do

if ( DTMCBU ) deallocate(buflux)
!if ( DTMCBU ) deallocate(buflux,bukapa)

!call omp_set_num_threads(28)

! Mass calculaiton
!if ( iscore ) then
!write(prt_dep,*)
!write(prt_dep,*), " - TOTAL & FUEL MASS  : ", tmass, fmass
!write(prt_dep,*)
!end if


! data sharing
do ii = 0, ncore-1
do imat = mp1(ii), mp2(ii)

    mat => materials(imat)
    call MPI_BCAST(mat%n_iso,1,MPI_INTEGER,ii,MPI_COMM_WORLD,ierr)
    niso = mat%n_iso

    if ( icore /= ii ) then
       deallocate(mat%ace_idx); allocate(mat%ace_idx(niso)); mat%ace_idx(:) = 0
       deallocate(mat%numden);  allocate(mat%numden(niso));  mat%numden(:)  = 0
       deallocate(mat%ogxs);    allocate(mat%ogxs(niso,4));  mat%ogxs(:,:)  = 0
   end if

   call MPI_BCAST(mat%ace_idx,niso,MPI_INTEGER,ii,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(mat%numden,niso,MPI_REAL8,ii,MPI_COMM_WORLD,ierr)

end do
end do


! predictor-corrector
if ( preco == 0 ) then
    istep_burnup = istep_burnup + 1
else if ( preco == 1 ) then
    if ( porc == 1 ) then
        if ( istep_burnup == nstep_burnup ) istep_burnup = istep_burnup + 1
        porc = 2
    else if ( porc == 2 ) then
        porc = 1
        istep_burnup = istep_burnup + 1
    end if
else if ( preco == 2 ) then
    istep_burnup = istep_burnup + 1
    ! time step
    stepsize0 = bstep_size
    ! flux
    mat%flux0 = mat%flux1
    ! cross section
    if ( allocated(mat%ogxs0) ) deallocate(mat%ogxs0)
    if ( allocated(mat%ace_idx0) ) deallocate(mat%ace_idx0)
    allocate(mat%ogxs0(size(mat%ogxs1,1),1:4))
    allocate(mat%ace_idx0(size(mat%ace_idx0)))
    mat%ogxs0 = mat%ogxs1
    mat%ace_idx0 = mat%ace_idx1
    print*, "ace_idx", allocated(mat%ace_idx0)
end if
 
end subroutine depletion 
    

subroutine TOTAL_MASS_CALC(mat)
    use ACE_HEADER, only: AceFormat
    implicit none
    type(Material_CE), intent(in), pointer:: mat
    real(8):: fmass_wise, tmass_wise, para
    integer:: ij

    fmass_wise = 0
    tmass_wise = 0
    do ij = 1, mat%n_iso
        if ( ace(mat%ace_idx(ij))%ZAID > 88000 ) &
        fmass_wise = fmass_wise + mat%numden(ij)*ace(mat%ace_idx(ij))%atn
        tmass_wise = tmass_wise + mat%numden(ij)*ace(mat%ace_idx(ij))%atn
    end do
    para = m_n *mat%vol / N_AVOGADRO
    fmass  = fmass + fmass_wise*para
    tmass  = tmass + tmass_wise*para

end subroutine



! ===================================================================
!         Chebyshev rational approximation method (CRAM) 
! ===================================================================
subroutine CRAM(burnupMat, nuc0, nuc1) 
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

    !Set control parameters
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
            call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
                idum, 1, iparm, msglvl, zdum, zdum, error)
        end if

        !Factorization
        phase = 22
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
            idum, 1, iparm, msglvl, zdum, zdum, error)

        !Back substitution and iterative refinement
        phase = 33
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
            idum, 1, iparm, msglvl, b, sol, error)
        nuc1 = nuc1 + 2.d0*dble(sol)
    end do

    !Termination and release memory used for factorization
    phase = -1 
    call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
        idum, 1, iparm, msglvl, zdum, zdum, error)

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
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
                idum, 1, iparm, msglvl, zdum, zdum, error)
        end if

        !Factorization
        phase = 22
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
                idum, 1, iparm, msglvl, zdum, zdum, error)

        !Back substitution and iterative refinement
        phase = 33
        call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
                idum, 1, iparm, msglvl, b, sol, error)
        nuc1 = nuc1 + 2.d0*dble(sol)
    end do

    !Termination and release memory used for factorization
    phase = -1 
    call pardiso(burn_pt, 1, 1, mtype, phase, nnuc, Acsr, iAcsr, jAcsr, &
            idum, 1, iparm, msglvl, zdum, zdum, error)


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
!  real ( kind = 8 ) r8_log_2
!  real ( kind = 8 ) r8mat_norm_li
!  integer ( kind = 4 ) s
!  real ( kind = 8 ) x(n,n)
!
!  a2(1:n,1:n) = a(1:n,1:n)
!
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
