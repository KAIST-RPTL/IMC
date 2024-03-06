module ace_xs

use constants, only : barn
use variables, only : E_mode, k_steady
use material_header 
use ace_header 
use ace_module 
use mpi

implicit none
    real(8), parameter:: inv_sqrt_pi = 0.564189583547756D0 ! 1/sqrt(pi)
    integer, allocatable :: mpiace(:,:)
    integer :: nace

contains

function getMacroXS_UEG(mat, erg, kT, urn) result (macro_xs)
    use constants, only : k_b
    implicit none
    type(Material_CE), intent(in) :: mat
    real(8), intent(in) :: erg, kT
    real(8), intent(in) :: urn(1:n_unr)
    real(8) :: macro_xs(5), macro_2(5)
    real(8) :: xs(4), micro_xs(6)

    real(8) :: dtemp ! OTF DB

    integer :: ierg, i, iso, iff, ierg0
    real(8) :: ipfac, ipfac0

    ! 1. Find ierg in UEG grid
    call getiueg(erg, ierg)

    ! 2. Find Interpolation factor
    if(ierg > nueg) ierg = nueg
    if(ierg < 1   ) ierg = 1
    ipfac = max(0D0,min(1D0,(erg-ueggrid(ierg))/(ueggrid(ierg+1)-ueggrid(ierg))))

    ! 3. Interpolate
    !print *, 'ALLOC?', allocated(mat % macro_ueg), trim(mat % mat_name)
    macro_xs(:) = (mat % macro_ueg(ierg,:) &
        + ipfac * (mat % macro_ueg(ierg+1,:) - mat % macro_ueg(ierg,:)))

!    macro_2 = getMacroXS(mat, erg, kT, urn)
!    print *, 'COMPARISON'
!    do i = 1, 5
!        print *, i, macro_xs(i), macro_2(i), erg
!    enddo

    ! 4. ADDITIONAL XS: URES
    if(n_unr == 0) return
    if(erg < Eumin .or. erg > Eumax) return

    do i = 1, n_unr
        iso = uresiso(i)
        !if(erg < ace(iso) % UNR % Emin .or. erg > ace(iso) % UNR % Emax) cycle
        ! MAC: Tot, Abs, Fis, NuF, QF 
        ! MIC: Tot, Ela, Sigd, Fis, Nuf, Qf
        if(.not. mat % ures(i)) cycle
        if(erg < ace(iso) % UNR % Emin .or. erg > ace(iso) % UNR % Emax) cycle

        micro_xs = 0D0
        micro_xs(1) = (ace(iso) % UEG % sigt(ierg) * (1D0-ipfac) &
            + ace(iso) % UEG % sigt(ierg+1) * ipfac)
        if(allocated(ace(iso)%UEG%sigel)) &
            micro_xs(2) = (ace(iso) % UEG % sigel(ierg) * (1D0-ipfac) &
            + ace(iso) % UEG % sigel(ierg+1) * ipfac)
        if(allocated(ace(iso)%UEG%sigd)) &
            micro_xs(3) = (ace(iso) % UEG % sigd(ierg) * (1D0-ipfac) &
            + ace(iso) % UEG % sigd(ierg+1) * ipfac)
        if(allocated(ace(iso)%UEG%sigf)) then
            micro_xs(4) = (ace(iso) % UEG % sigf(ierg) * (1D0-ipfac) &
                + ace(iso) % UEG % sigf(ierg+1) * ipfac)
            micro_xs(5) = micro_xs(4) * dble(getnu(iso,erg))
        endif

        if(erg >= ace(iso) % UNR % Emin .and. erg <= ace(iso) % UNR % Emax) then
            call URES_PTABLE(iso, erg, xs, iff, urn(ace(iso) % UNR % unridx))
            micro_xs(1) = xs(1)
            micro_xs(2) = xs(2)
            micro_xs(4) = xs(3)
            micro_xs(5) = xs(3) * getnu(iso,erg)
            micro_xs(3) = xs(4)
            !micro_xs(1) = xs(2) + xs(3) + xs(4)
        endif

        micro_xs(6) = micro_xs(4) * ace(iso) % qval

        macro_xs(1) = macro_xs(1) + micro_xs(1) * mat % numden( mat % uresidx(i) ) * barn
        macro_xs(2) = macro_xs(2) + micro_xs(3) * mat % numden( mat % uresidx(i) ) * barn
        if(micro_xs(4) == 0) cycle
        macro_xs(3) = macro_xs(3) + micro_xs(4) * mat % numden( mat % uresidx(i) ) * barn
        macro_xs(4) = macro_xs(4) + micro_xs(4) * mat % numden( mat % uresidx(i) ) * barn * getnu(iso, erg)
        macro_xs(5) = macro_xs(5) + micro_xs(4) * mat % numden( mat % uresidx(i) ) * barn * ace(iso) % qval

    enddo

    ! TODO: SAB, OTFDB, URES

end function

function getMacroXS (mat, erg,kT, urn) result (macro_xs)
    use CONSTANTS, only: K_B
    implicit none
    type(Material_CE), intent(in) :: mat
    real(8), intent(in) :: erg
    real(8), intent(in) :: kT
    real(8), intent(in) :: urn(1:n_unr)
    real(8) :: macro_xs(5)
    
    integer :: i 
    integer :: i_iso, iso_, ierg_
    integer :: pt1, pt2, pt3, pt4
    real(8) :: ipfac
    real(8) :: micro_t, micro_d, micro_f, micro_nuf, micro_a, micro_el, micro_xn
    real(8) :: macro_t, macro_f, macro_nuf, macro_a, macro_qf
    real(8) :: xn_xs(4)
    real(8) :: xs(5), xs_tmp
    integer :: isab, iff
    real(8) :: dtemp    ! temperautre difference | library - material |
    integer :: isab_l, isab_h
    macro_t   = 0.0d0
    macro_a   = 0.0d0
    macro_f   = 0.0d0
    macro_nuf = 0.0d0
    macro_qf  = 0.0d0
    xn_xs(:)  = 0.0d0

    !print *, mat%mat_name, mat%n_iso
    MAT_ISO_LOOP: do i_iso = 1, mat%n_iso     ! isotope number in the material
    
        iso_ = mat%ace_idx(i_iso)   ! isotope number in the inputfile

        ! =====================================================================
        ! S(a,b) treatment
        ! isab = ace(iso_)%sab_iso    ! isotope number for S(a,b)
        isab = mat % sablist(i_iso)
        if ( isab /= 0 .and. erg < 4D-6 ) then
            if(isab > 0) then
                call GET_SAB_MAC(mat%numden(i_iso),iso_,isab,erg,macro_t,macro_a)
                cycle MAT_ISO_LOOP
            elseif(isab < 0) then ! MODER
                if(.not.allocated(therm)) print *, 'NOTHERM XS'
                isab_l = therm(-isab) % iso_low
                isab_h = therm(-isab) % iso_high
                call GET_SAB_MAC(mat % numden(i_iso) * (1d0-therm(-isab) % f), &
                   iso_, isab_l, erg, macro_t, macro_a)
                call GET_SAB_MAC(mat % numden(i_iso) * (therm(-isab) % f), &
                   iso_, isab_h, erg, macro_t, macro_a)
                cycle MAT_ISO_LOOP
            endif
        end if

        ! =====================================================================
        ! On-the-fly Doppler broadening
        dtemp = abs(ace(iso_)%temp-kT)
        !@print *, 'TEMP:', iso_, trim(mat%mat_name), trim(ace(iso_)%xslib), ace(iso_)%temp/K_B, mat%temp/K_B
        !print *,  dtemp, K_B, kT, ace(iso_)%temp, ace(iso_)%zaid
        if ( mat%db .and. ( dtemp > K_B .and. erg < 1d0 ) ) then
            call GET_OTF_DB_MAC(mat%numden(i_iso), i_iso,iso_,erg,xs,dtemp)

            call getierg(iso_,ierg_,erg)
            
            ipfac = max(0.d0, min(1.d0,(erg-ace(iso_)%E(ierg_))/(ace(iso_)%E(ierg_+1)-ace(iso_)%E(ierg_))))
            
            !==============================================================
            ! Microscopic XS
            micro_t   = ace(iso_)%sigt(ierg_) + ipfac*(ace(iso_)%sigt(ierg_+1)-ace(iso_)%sigt(ierg_))
            micro_d   = ace(iso_)%sigd(ierg_) + ipfac*(ace(iso_)%sigd(ierg_+1)-ace(iso_)%sigd(ierg_))

!            print *, 'COMPARISON:', mat % temp / K_B, erg, trim(mat % mat_name), ace(iso_) % zaid
!            print *, 'TOT', xs(1), micro_t * mat % numden(i_iso) * barn  
!            print *, 'GAM', xs(2), micro_d * mat % numden(i_iso) * barn
            macro_t   = macro_t   + xs(1) 
            macro_a   = macro_a   + xs(2) 
            macro_f   = macro_f   + xs(3)
            macro_nuf = macro_nuf + xs(4)
            macro_qf  = macro_qf  + xs(5)
            cycle MAT_ISO_LOOP
        end if

        call getierg(iso_,ierg_,erg)
        
        ipfac = max(0.d0, min(1.d0,(erg-ace(iso_)%E(ierg_))/(ace(iso_)%E(ierg_+1)-ace(iso_)%E(ierg_))))
        
        !==============================================================
        ! Microscopic XS
        micro_t   = ace(iso_)%sigt(ierg_) + ipfac*(ace(iso_)%sigt(ierg_+1)-ace(iso_)%sigt(ierg_))
        micro_d   = ace(iso_)%sigd(ierg_) + ipfac*(ace(iso_)%sigd(ierg_+1)-ace(iso_)%sigd(ierg_))
        micro_f   = 0.d0
        micro_nuf = 0.d0
        micro_a   = micro_d
        
        
        ! Fissionable Material
        !if(ace(iso_)%jxs(21)/=0) then
        !if(allocated(ace(iso_)%sigf)) then 
        if(ace(iso_)%jxs(21)/=0 .or. allocated(ace(iso_)%sigf)) then
            micro_f   = ace(iso_)%sigf(ierg_) + ipfac*(ace(iso_)%sigf(ierg_+1)-ace(iso_)%sigf(ierg_))
            micro_nuf = getnu(iso_,erg)*micro_f
            micro_a   = micro_d + micro_f
        endif
        !micro_el = ace(iso_)%sigel(ierg_) + ipfac*(ace(iso_)%sigel(ierg_+1)-ace(iso_)%sigel(ierg_))
       
        ! Apply URES
        if  (ace(iso_) % UNR % URES .and. mat%numden(i_iso) > ures_cut) then
            ! Only for Energy between Emin and Emax
            if ( erg >= ace(iso_) % UNR % Emin .and. &
                    erg <= ace(iso_) % UNR % Emax ) then
                call URES_PTABLE(iso_, erg, xs, iff, urn(ace(iso_) % UNR % unridx))
                micro_t   = xs(1)
                if(allocated(ace(iso_)%sigf)) then
                    micro_f   = xs(3)
                    micro_nuf = micro_f * getnu(iso_,erg)
                else
                    micro_f = 0d0; micro_nuf = 0d0
                endif
                micro_d   = xs(4)
                micro_a   = micro_d + micro_f
            endif
        endif

        !>Summation for macroscopic cross sections
        macro_t   = macro_t   + mat%numden(i_iso) * micro_t   * barn
        macro_a   = macro_a   + mat%numden(i_iso) * micro_a   * barn
        macro_f   = macro_f   + mat%numden(i_iso) * micro_f   * barn
        macro_nuf = macro_nuf + mat%numden(i_iso) * micro_nuf * barn
        macro_qf  = macro_qf  + mat%numden(i_iso) * micro_f   * barn * ace(iso_)%qval

        !> Macro_xs of Sig_abs is only used for FMFD, (n,xn) XS is subtracted.
        do i = 1, ace(iso_)%NXS(5) !> through the reaction types...
            pt1 = abs(ace(iso_)%TY(i))
            if (pt1 > 1 .and. pt1 < 5) then 
                if ( dtemp > K_B * 1e-2 .and. mat % db ) then
                    call GET_OTF_DB_MT(kT, iso_, erg, i, micro_xn)
                else
                    micro_xn   = ace(iso_)%sig_MT(i)%cx(ierg_) & 
                                + ipfac*(ace(iso_)%sig_MT(i)%cx(ierg_+1) - ace(iso_)%sig_MT(i)%cx(ierg_))
                endif
                xn_xs(pt1) = xn_xs(pt1) + mat%numden(i_iso) * micro_xn * barn
            endif
        enddo

    enddo MAT_ISO_LOOP 

    do i = 2, 4 
        macro_a = macro_a - (dble(i)-1.0d0)*xn_xs(i)
    enddo     
    
    macro_xs(1) = macro_t  
    macro_xs(2) = macro_a  
    macro_xs(3) = macro_f  
    macro_xs(4) = macro_nuf
    macro_xs(5) = macro_qf
    
    
end function

function getMicroXS (iso, erg) result (micro_xs)
    integer, intent(in) :: iso
    real(8), intent(in) :: erg
    real(8) :: micro_xs(6)
    
    integer :: ierg_, i 
    integer :: pt1, pt2, pt3, pt4
    real(8) :: ipfac
    real(8) :: micro_t, micro_d, micro_f, micro_nuf, micro_a, micro_el

    real(8) :: xs(4)
    integer :: iff
    
    call getierg(iso,ierg_,erg)
    
    ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg_))/(ace(iso)%E(ierg_+1)-ace(iso)%E(ierg_))))
    !==============================================================
    ! Microscopic XS
    micro_t   = ace(iso)%sigt(ierg_) + ipfac*(ace(iso)%sigt(ierg_+1)-ace(iso)%sigt(ierg_))
    micro_d   = ace(iso)%sigd(ierg_) + ipfac*(ace(iso)%sigd(ierg_+1)-ace(iso)%sigd(ierg_))
    micro_f   = 0.d0
    micro_nuf = 0.d0
    micro_a   = micro_d
    
    ! Fissionable Material
    if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
        micro_f   = ace(iso)%sigf(ierg_) + ipfac*(ace(iso)%sigf(ierg_+1)-ace(iso)%sigf(ierg_))
        micro_nuf = getnu(iso,erg)*micro_f
        micro_a   = micro_d + micro_f
    endif
    micro_el = ace(iso)%sigel(ierg_) + ipfac*(ace(iso)%sigel(ierg_+1)-ace(iso)%sigel(ierg_))
    
    micro_xs(1) = micro_t
    micro_xs(2) = micro_el
    micro_xs(3) = micro_a
    micro_xs(4) = micro_f
    micro_xs(5) = micro_nuf
    micro_xs(6) = micro_d
    
    !print '(5F10.4)', micro_t, micro_el, micro_a, micro_f
    
    
end function

function getxs (mt_ENDF,iso, erg, ierg)
    integer, intent(in) :: mt_ENDF, iso
    real(8), intent(in) :: erg
    real(8) :: getxs
    real(8) :: ipfac
    integer :: i, iMT
    integer, optional :: ierg
    type (CrossSectionDataForm), pointer :: sigmt
    
    
    ! 1. Find MT index 
    iMT = 0; getxs = 0.0d0
    do i = 1, ace(iso)%NXS(4) 
        if (ace(iso)%MT(i) == mt_ENDF) then 
            iMT = i 
            exit
        endif 
    enddo 
    if (iMT == 0) return  ! no such reaction 
    
    ! 2. Find erg grid index 
    if (.not. present(ierg)) then 
        call getierg(iso,ierg,erg)
    endif
    
    
    ! 3. Calculate xs_MT
    ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
    sigmt => ace(iso)%sig_MT(iMT)
    if (ierg < (sigmt%IE+sigmt%NE-1) .and. ierg >= sigmt%IE) then 
        getxs = sigmt%cx(ierg) + ipfac*(sigmt%cx(ierg+1)-sigmt%cx(ierg))
    endif
    
    return 

end function



function getnu (iso_,erg0) result (nu)
    integer, intent(in) :: iso_
    real(8), intent(in) :: erg0
    real(8) :: nu
    integer :: nublock, n, i, NC, ierg, pt1, pt2, pt3 
    real(8) :: ipfac
    type (AceFormat), pointer :: ac
    
    nu = 0 
    ac => ace(iso_)
    if(ac % nu_block_exist == .false.) return
    
    select case( ac % nu_tot_flag )
    case(1) !> polynomial function form
        nu = 0.0d0
        do i = 1, ac % nu_tot % NC
            nu = nu + ac%nu_tot%C(i)*erg0**(i-1)
        enddo
        
    case(2) !> tabular data form
        !ac%nu_tot%E(:) !> nu energy grid 
        !ac%nu_tot%F(:) !> corresponding nu
        
        !if (ac%nu_tot%NR /= 0) print *, "WARNING :: nu is not lin-lin for", ace(iso_)%library
        ! 1. binary search to find ierg of erg in E(:) 
        pt1 = 1
        pt2 = ac % nu_tot % NE
        Do
            if(pt2 - pt1 == 1) exit
            pt3 = (pt2 + pt1)/2
            if(erg0 >= ac%nu_tot%E(pt3)) then
              pt1 = pt3
            else
              pt2 = pt3
            endif
        Enddo
        ierg = pt1 !store low bound energy index 
        
        !if (ac % nu_tot % NR /= 0) print *, "ENDF interpolation required for Nu "
        ! 2. calculate interpolation factor
        ipfac = max(0.d0, min(1.d0,(erg0-ac%nu_tot%E(ierg))/(ac%nu_tot%E(ierg+1)-ac%nu_tot%E(ierg))))
        
        ! 3. linear-linear interpolation
        nu = ac%nu_tot%F(ierg) + ipfac*(ac%nu_tot%F(ierg+1)-ac%nu_tot%F(ierg))

    end select
    
	if (do_transient) nu = nu / k_steady 
	
end function

function getnudel (iso_,erg0) result (nu)
    integer, intent(in) :: iso_
    real(8), intent(in) :: erg0
    real(8) :: nu
    integer :: ierg, pt1, pt2, pt3 
    real(8) :: ipfac
    type (AceFormat), pointer :: ac
    
    nu = 0 
    ac => ace(iso_)
    if(ac % nu_del_block_exist == .false.) return
    
	!if (ac%nu_del%NR /= 0) then 
	!	print *, "WARNING :: nu_del is not lin-lin for", ace(iso_)%library
	!	print *, ace(iso_)%nu_del%NR
	!	print *, ac%nu_del%NR
	!	print *, iso_
	!endif 
	! 1. binary search to find ierg of erg in E(:) 
	pt1 = 1
	pt2 = ac % nu_del % NE
	Do
		if(pt2 - pt1 == 1) exit
		pt3 = (pt2 + pt1)/2
		if(erg0 >= ac%nu_del%E(pt3)) then
		pt1 = pt3
		else
		pt2 = pt3
		endif
	Enddo
	ierg = pt1 !store low bound energy index 
	
	if (ac % nu_del % NR /= 0) print *, "ENDF interpolation required for Nu "
	! 2. calculate interpolation factor
	ipfac = max(0.d0, min(1.d0,(erg0-ac%nu_del%E(ierg))/(ac%nu_del%E(ierg+1)-ac%nu_del%E(ierg))))
	
	! 3. linear-linear interpolation
	nu = ac%nu_del%F(ierg) + ipfac*(ac%nu_del%F(ierg+1)-ac%nu_del%F(ierg))
	
end function




subroutine getierg(iso_,ierg_,erg0)
    implicit none
    real(8), intent(in) :: erg0
    integer, intent(in) :: iso_
    integer, intent(out) :: ierg_

    integer :: pt1, pt2, pt3, uidx


    !Energy grid search algorithm
    !erg0 = incidient neutron energy in lab system
    !iso_ = collision isotope index
    !ierg_ = low bound energy grid index
    
    
    uidx = 1 + int(log10(erg0/Emin)/udelta)
	
	if (uidx > 8192) then 
		print *, uidx
		print *, erg0, Emin, udelta
        uidx = 8192
	endif 
    if(uidx<1) uidx = 1
    pt1 = ugrid(uidx-1,iso_)
    pt2 = min(ugrid(uidx,iso_)+1,ace(iso_)%nxs(3))

    !pt1 = 1 
    !pt2 = size(ace(iso_)%E)
    
    if(pt1==pt2) then
      pt1 = pt1 - 1
    else
      Do
        if(pt2 - pt1 == 1) exit
        pt3 = (pt2 + pt1)/2
        if(erg0 > ace(iso_)%E(pt3)) then
          pt1 = pt3
        elseif(erg0 < ace(iso_)%E(pt3)) then
          pt2 = pt3
        else
        ierg_ = pt3; return
        endif
      Enddo
    endif
    ierg_ = pt1 !store low bound energy index 
    
    !if(erg0 < 1.d-11) print *, erg0, ierg_ !, ace(iso_)%E(ierg_), ace(iso_)%E(ierg_+1)
end subroutine 

subroutine getiueg(erg, ierg)
    implicit none
    real(8), intent(in) :: erg
    integer, intent(out):: ierg

    integer :: pt1, pt2, pt3
    integer :: uidx

    uidx = 1 + int(log10(erg/ueggrid(1))/unidel)

    if ( uidx > nuni ) then
        print *, 'EXCEED', erg, uidx
        uidx = nuni; !erg = ueggrid(nueg)
        ierg = nueg;
        return
    elseif ( uidx < 1 ) then
        !print *, 'BELOW', erg, uidx
        !uidx = 1; !erg = ueggrid(1)
        ierg  = 0;
        return
    endif

    pt1 = max(unigrid(uidx-1), 0)
    pt2 = min(unigrid(uidx)+1,nueg)

    if(ueggrid(pt1) > erg) print *,'LO', ueggrid(pt1), int(log10(erg/ueggrid(1))/unidel), erg
    if(ueggrid(pt2) < erg) print *,'HI', ueggrid(pt2), erg

    if(pt1==pt2) then
        pt1 = pt1 - 1
    else
        do
        if(pt2-pt1==1) exit
        pt3 = (pt2 + pt1) / 2
        if(erg > ueggrid(pt3)) then
            pt1 = pt3
        elseif(erg < ueggrid(pt3)) then
            pt2 = pt3
        else
            ierg = pt3
            return
        endif
        enddo
    endif
    ierg = pt1
    if(ueggrid(ierg+1) <= erg) print *, 'WTF?', ierg, ueggrid(ierg), ueggrid(ierg+1)
end subroutine

subroutine setugrid
    implicit none
    real(8) :: Etmp
    integer :: totngrid
    integer :: i, j, k, iso_, idx
    integer :: pt1
    real(8), allocatable :: tmpgrid(:)

    if(E_mode==0) return
    !Set ugrid to accelerate energy-grid search
    allocate(ugrid(0:nugrid,1:num_iso))
    !print *, 'ugrid size', nugrid, num_iso
    do iso_ = 1, num_iso
      ugrid(0,iso_) = 1
      ugrid(nugrid,iso_) = ace(iso_)%nxs(3)
    end do

    Emin = 1d-11
    udelta = log10((Emax+Emin*5d-1)/Emin)/dble(nugrid)

    do iso_ = 1, num_iso
      idx = 1
      do i=1, nugrid-1
        Etmp = Emin*10.d0**(dble(i)*udelta)
        if(Etmp > ace(iso_)%E( ace(iso_)%NXS(3) )) then 
          idx = ace(iso_)%nxs(3)
          go to 10  
        end if
        do
          if(Etmp < ace(iso_)%E(idx)) go to 10
          idx = idx + 1
        end do
10        ugrid(i,iso_) = idx - 1
      end do
    enddo

    if(icore==score) print *, "   Setting ugrid..."
end subroutine 

subroutine setuegrid
    implicit none
    real(8) :: Etmp, tolerance
    integer :: totngrid
    integer :: i, j, k, iso_, idx
    integer :: pt, pt1, pt2, pt3, pt4
    real(8), allocatable :: tmpgrid(:), tmpgrid_2(:), heaps(:), tmpgrid_3(:), tmpgrid_4(:), sabpts(:), thresh(:)

    if(E_mode==0) return
    totngrid = 0
    !Set ugrid to accelerate energy-grid search
    Emin = 1.d-11
    !print *, 'ugrid size', nugrid, num_iso
    do iso_ = 1, num_iso
      totngrid = totngrid + ace(iso_)%NXS(3)
      if(ace(iso_) % UNR % URES) & ! URR case
          totngrid = totngrid + ace(iso_) % UNR % N
    end do

    ! SAB case !TODO
    if(sab_iso /= 0) then
        do iso_ = 1, sab_iso
            totngrid = totngrid + sab(iso_) % NXS(3)
        enddo
    endif
    allocate(tmpgrid(1:totngrid)); pt1 = 1
    allocate(heaps(1:totngrid)); pt2 = 0
    allocate(sabpts(1:totngrid)); pt3 = 0
    allocate(thresh(1:totngrid)); pt4 = 0

    udelta = log10((Emax+1E-9)/Emin)/dble(nugrid)

    do iso_ = 1, num_iso
      tmpgrid(pt1:pt1-1+ace(iso_)%NXS(3)) = ace(iso_) % E(:)
      pt1 = pt1 + ace(iso_)%NXS(3)

      do i = 1, ace(iso_) % NXS(3) ! Conserve Heaps
        if( i == 1 .or. i == ace(iso_)%NXS(3) ) then
            pt2 = pt2 + 1
            heaps(pt2) = ace(iso_) % E(i)
        elseif ( ace(iso_) % sigd(i-1) < ace(iso_) % sigd(i) .and. &
                ace(iso_) % sigd(i+1) < ace(iso_) % sigd(i) ) then
            pt2 = pt2 + 1
            heaps(pt2) = ace(iso_) % E(i)
        elseif ( ace(iso_) % sigd(i-1) > ace(iso_) % sigd(i) .and. &
                ace(iso_) % sigd(i+1) > ace(iso_) % sigd(i) ) then
            pt2 = pt2 + 1
            heaps(pt2) = ace(iso_) % E(i)
        elseif ( allocated(ace(iso_) % sigf) ) then
            if ( ace(iso_) % sigf(i-1) < ace(iso_) % sigf(i) .and. &
                    ace(iso_) % sigf(i+1) < ace(iso_) % sigf(i) ) then
                pt2 = pt2 + 1
                heaps(pt2) = ace(iso_) % E(i)
            elseif ( ace(iso_) % sigf(i-1) > ace(iso_) % sigf(i) .and. &
                    ace(iso_) % sigf(i+1) > ace(iso_) % sigf(i) ) then
                pt2 = pt2 + 1
                heaps(pt2) = ace(iso_) % E(i)
            endif
        endif
      enddo

      do i = 1, ace(iso_) % NXS(4) ! Number of RX: Threshold
          pt4 = pt4 + 1
          if( ace(iso_) % sig_MT(i) % IE > 1 ) &
              thresh(pt4) = ace(iso_) % E( ace(iso_) % sig_MT(i) % IE )
      enddo
      
      !URES
      if(ace(iso_) % UNR % URES) then
          tmpgrid(pt1:pt1-1+ace(iso_)%UNR%N) = ace(iso_) % UNR % E(:)
          pt1 = pt1 + ace(iso_) % UNR % N
      endif

      !SAB !TODO
    enddo

    do iso_ = 1, sab_iso
        sabpts(pt3:pt3-1+sab(iso_) % itie % ne) = sab(iso_) % itie % erg(:)
        pt3 = pt3 + sab(iso_) % itie % ne

        if ( sab(iso_) % jxs (4) /= 0 ) then
            sabpts(pt3 : pt3-1 + sab(iso_) % itce % ne) = sab(iso_) % itce % erg(:)
            pt3 = pt3 + sab(iso_) % itce % ne
        endif
    enddo

    if(icore==score) print *, "   Setting UNIONIZED GRID..."
    if(icore==score) print *, 'HEAP #:', pt2
    if(icore==score) print *, 'SAB  #:', pt3
    if(icore==score) print *, 'THRS #:', pt4

    ! SORT and COLLIDE UEGGRID 
    ! 1. SORT
    call QUICKSORT(tmpgrid,1,totngrid)

    ! 2. COLLECT UNIQUEs
    allocate(tmpgrid_2(0:totngrid))
    idx = 1
    tmpgrid_2(0)   = 0d0
    tmpgrid_2(idx) = tmpgrid(1)
    tolerance = 0d0
    do i = 2, totngrid
        !if(tmpgrid(i)/=tmpgrid(i-1) .and. tmpgrid(i)<Emax &
        !    )then
        if(tmpgrid(i) < Emax) then
        if(abs(tmpgrid(i)-tmpgrid(i-1))>tmpgrid(i-1) * tolerance) then
        
            idx = idx + 1
            tmpgrid_2(idx) = tmpgrid(i)
        endif
        endif
    enddo

    pt = idx + pt2 + pt3 + pt4
    allocate(tmpgrid_3(1:idx+pt2))
    tmpgrid_3(1:idx) = tmpgrid_2(1:idx)
    tmpgrid_3(idx+1:idx+pt2) = heaps(1:pt2)
    tmpgrid_3(idx+pt2+1:idx+pt2+pt3) = sabpts(1:pt3)
    tmpgrid_3(idx+pt2+pt3+1:pt) = thresh(:)

    open(502, file='ueg.out', action='write', status='unknown')
    do i = 1, idx
        write(502, *) i, tmpgrid_2(i)
    enddo

    deallocate(tmpgrid, tmpgrid_2)
    call quicksort(tmpgrid_3, 1, idx+pt2)
    if(icore==score) print *, 'TMPGRID', tmpgrid_3(1:10)

    allocate(tmpgrid_4(0:idx+pt2))
    idx = 0
    tmpgrid_4   = 0d0
    do i = 1, pt
        if( i == 1 ) then
            if( tmpgrid_3(i) > 0 ) then
                idx = idx + 1
                tmpgrid_4(idx) = tmpgrid_3(i)
            endif
        elseif(tmpgrid_3(i)/=tmpgrid_3(i-1) .and. tmpgrid_3(i)<Emax .and. tmpgrid_3(i) > 0d0 ) then
        
            idx = idx + 1
            tmpgrid_4(idx) = tmpgrid_3(i)
        endif
    enddo

    nueg    = idx
    allocate(ueggrid(0:nueg))
    ueggrid(0:nueg) = tmpgrid_4(0:nueg)

    deallocate(tmpgrid_3, tmpgrid_4)

    if(icore==score)print *, 'NUEG', nueg
    do i = 1, nueg
        write(502, *) i, ueggrid(i), log(ueggrid(i))
    enddo

    ! 3. SETUP HASH TABLE
    !   NUNI = len(UNIGRID)
    Emin = ueggrid(1); Emax = ueggrid(nueg)
    unidel = log10((Emax+Emin*1d-1)/Emin)/dble(nuni)

    allocate(unigrid(0:nuni))

    idx = 1
    unigrid(0) = 0
    do i = 1, nuni-1
        Etmp = Emin * 1d1 ** (dble(i) * unidel)
        if(Etmp > ueggrid(nueg)) then
            idx = nuni
            goto 22
        endif
        do
            if(Etmp < ueggrid(idx)) goto 22
            idx = idx + 1
        enddo
22      unigrid(i) = idx - 1
    enddo
    unigrid(nuni) = Emax

    if(icore==score) print *, 'NUNI:',nuni,'UNIDEL:',unidel,Emin,Emax

    if(icore==score) write(*,'(A,I8,A)') '   UNIONIZED GRID SET: ', nueg, 'grids'
end subroutine

!subroutine setuegrid
!    implicit none
!    real(8) :: Etmp
!    integer :: totngrid
!    integer :: i, j, k, iso_, idx
!    integer :: pt1
!    real(8), allocatable :: tmpgrid(:)
!
!    if(E_mode==0) return
!    totngrid = 0
!    !Set ugrid to accelerate energy-grid search
!    Emin = 1.d-11
!    !print *, 'ugrid size', nugrid, num_iso
!    do iso_ = 1, num_iso
!      totngrid = totngrid + ace(iso_)%NXS(3)
!      if(ace(iso_) % UNR % URES) & ! URR case
!          totngrid = totngrid + ace(iso_) % UNR % N
!    end do
!
!    ! SAB case !TODO
!!    if(sab_iso /= 0) then
!!        do iso_ = 1, sab_iso
!!            totngrid = totngrid + sab(iso_) % NXS(3)
!!        enddo
!!    endif
!    allocate(tmpgrid(1:totngrid)); pt1 = 1
!
!    udelta = log10((Emax+1E-9)/Emin)/dble(nugrid)
!
!    do iso_ = 1, num_iso
!      tmpgrid(pt1:pt1-1+ace(iso_)%NXS(3)) = ace(iso_) % E(:)
!      pt1 = pt1 + ace(iso_)%NXS(3)
!      
!      !URES
!      if(ace(iso_) % UNR % URES) then
!          tmpgrid(pt1:pt1-1+ace(iso_)%UNR%N) = ace(iso_) % UNR % E(:)
!          pt1 = pt1 + ace(iso_) % UNR % N
!      endif
!
!      !SAB !TODO
!    enddo
!
!    if(icore==score) print *, "   Setting UNIONIZED GRID..."
!
!    ! SORT and COLLIDE UEGGRID 
!    ! 1. SORT
!    call QUICKSORT(tmpgrid,1,totngrid)
!
!    ! 2. COLLECT UNIQUEs
!    allocate(ueggrid(0:totngrid))
!    idx = 1
!    ueggrid(0)   = 0d0
!    ueggrid(idx) = tmpgrid(1)
!    do i = 2, totngrid
!        if(tmpgrid(i)/=tmpgrid(i-1) .and. tmpgrid(i)<UEGMAX &
!            )then
!            !.and. tmpgrid(i)-tmpgrid(i-1)>=5E-5*tmpgrid(i-1)) then
!            idx = idx + 1
!            ueggrid(idx) = tmpgrid(i)
!        endif
!    enddo
!    nueg    = idx
!    ueggrid = ueggrid(1:nueg)
!
!    ! 3. SETUP HASH TABLE
!    !   NUNI = len(UNIGRID)
!    Emin = ueggrid(1); Emax = ueggrid(nueg)
!    unidel = log10((Emax+Emin*1d-1)/Emin)/dble(nuni)
!
!    allocate(unigrid(0:nuni))
!
!    idx = 1
!    unigrid(0) = 0
!    do i = 1, nuni-1
!        Etmp = Emin * 1d1 ** (dble(i) * unidel)
!        if(Etmp > nueg) then
!            idx = nuni
!            goto 22
!        endif
!        do
!            if(Etmp < ueggrid(idx)) goto 22
!            idx = idx + 1
!        enddo
!22      unigrid(i) = idx - 1
!        if(icore==score) print *, 'hash', i, unigrid(i), Etmp
!    enddo
!    unigrid(nuni) = UEGMAX
!
!    if(icore==score) print *, 'NUNI:',nuni,'UNIDEL:',unidel,Emin,Emax
!
!    if(icore==score) write(*,'(A,I8,A)') '   UNIONIZED GRID SET: ', nueg, 'grids'
!end subroutine

! =============================================================================
! GET_SAB_MAC
! =============================================================================
subroutine GET_SAB_MAC(nd,iiso,isab,erg,xs_t,xs_a)
    real(8), intent(in):: nd            ! number density
    integer, intent(in):: iiso, isab    ! index for isotope & S(a,b)
    real(8), intent(in):: erg           ! energy
    real(8), intent(inout):: xs_t, xs_a ! cross section
    type(SAB_INEL_XS), pointer:: abi
    type(SAB_EL_XS), pointer:: abe
    real(8):: micro_t, micro_i, micro_e, micro_a
    integer:: ierg
    real(8):: ipfac

    abi => sab(isab)%itie
    abe => sab(isab)%itce

    ! absorption 
    call getierg(iiso,ierg,erg)
    ipfac = max(0D0,min(1D0,(erg-ace(iiso)%e(ierg)) &
        /(ace(iiso)%E(ierg+1)-ace(iiso)%E(ierg))))
    micro_a = ace(iiso)%sigd(ierg) + & 
        ipfac*(ace(iiso)%sigd(ierg+1)-ace(iiso)%sigd(ierg))

    ! inelastic
    call GET_IERG_SABI(isab,ierg,erg)
    ipfac = max(0D0,min(1D0,(erg-abi%erg(ierg)) &
        /(abi%erg(ierg+1)-abi%erg(ierg))))
    micro_i = abi%xs(ierg) + ipfac*(abi%xs(ierg+1)-abi%xs(ierg))

    ! elastic
    micro_e = 0D0
    if ( sab(isab)%jxs(4) /= 0 ) then
    call GET_IERG_SABE(isab,ierg,erg)
    ipfac = max(0D0,min(1D0,(erg-abe%erg(ierg)) &
        /(abe%erg(ierg+1)-abe%erg(ierg))))
    micro_e = abe%xs(ierg) + ipfac*(abe%xs(ierg+1)-abe%xs(ierg))
    if ( sab(isab)%nxs(5) == 4 ) micro_e = micro_e / abe%erg(ierg)
    if ( abe % erg(ierg) > erg ) micro_e = 0d0
    end if

    if ( associated(abi) ) nullify(abi)
    if ( associated(abe) ) nullify(abe)

    micro_t = micro_i + micro_e + micro_a

    xs_t = xs_t + nd * micro_t * barn
    xs_a = xs_a + nd * micro_a * barn

end subroutine


! =============================================================================
! GET_SAB_MIC
! =============================================================================
subroutine GET_SAB_MIC(mat,imat,erg,xs)
    type(Material_CE), intent(in):: mat
    integer, intent(in):: imat
    real(8), intent(in):: erg
    real(8), intent(inout):: xs(:)
    type(SAB_INEL_XS), pointer:: abi
    type(SAB_EL_XS), pointer:: abe
    integer:: iiso, isab, ierg, isab_l, isab_h
    real(8) :: f
    real(8):: ipfac
    real(8) :: xs2l, xs2h, xs6l, xs6h

    if ( erg > 4E-6 ) return
    iiso = mat%ace_idx(imat)
    isab = mat%sablist(imat)
    if ( isab == 0 ) return
    ! total / elastic / absorption
    call getierg(iiso,ierg,erg)
    ipfac = max(0D0,min(1D0,(erg-ace(iiso)%e(ierg)) &
        /(ace(iiso)%E(ierg+1)-ace(iiso)%E(ierg))))
    xs(1) = ace(iiso)%sigt(ierg) + & 
        ipfac*(ace(iiso)%sigt(ierg+1)-ace(iiso)%sigt(ierg))
    xs(2) = ace(iiso)%sigel(ierg) + & 
        ipfac*(ace(iiso)%sigel(ierg+1)-ace(iiso)%sigel(ierg))
    xs(3) = ace(iiso)%sigd(ierg) + & 
        ipfac*(ace(iiso)%sigd(ierg+1)-ace(iiso)%sigd(ierg))
    xs(4:5) = 0D0

    xs(1) = xs(1) - xs(2)

    if( isab > 0 ) then ! S(a,b)
        abi => sab(isab)%itie
        abe => sab(isab)%itce
    
        ! thermal inelastic
        call GET_IERG_SABI(isab,ierg,erg)
        ipfac = max(0D0,min(1D0,(erg-abi%erg(ierg)) &
            /(abi%erg(ierg+1)-abi%erg(ierg))))
        xs(2) = abi%xs(ierg) + ipfac*(abi%xs(ierg+1)-abi%xs(ierg))
    
        ! thermal elastic
        xs(6) = 0D0
        if ( sab(isab)%jxs(4) /= 0 ) then
        call GET_IERG_SABE(isab,ierg,erg)
        ipfac = max(0D0,min(1D0,(erg-abe%erg(ierg)) &
            /(abe%erg(ierg+1)-abe%erg(ierg))))
        xs(6) = abe%xs(ierg) + ipfac*(abe%xs(ierg+1)-abe%xs(ierg))
        if ( sab(isab)%nxs(5) == 4 ) xs(6) = xs(6) / abe%erg(ierg)
        if ( abe % erg ( ierg ) > erg ) xs(6) = 0d0
        end if
    
        xs(2) = xs(2) + xs(6)  ! thermal scattering = inelastic + elastic
        xs(1) = xs(1) + xs(2)  ! total += thermal scattering
    
        if ( associated(abi) ) nullify(abi)
        if ( associated(abe) ) nullify(abe)
    else ! MODER
        if(.not.allocated(therm)) print *, 'NOTHERM XS2', isab
        isab_l = therm(-isab) % iso_low
        isab_h = therm(-isab) % iso_high
        f      = therm(-isab) % f

        abi => sab(isab_l)%itie
        abe => sab(isab_l)%itce

        ! LOW TEMP:
        ! thermal inelastic
        call GET_IERG_SABI(isab_l,ierg,erg)
        ipfac = max(0D0,min(1D0,(erg-abi%erg(ierg)) &
            /(abi%erg(ierg+1)-abi%erg(ierg))))
        xs2l = abi%xs(ierg) + ipfac*(abi%xs(ierg+1)-abi%xs(ierg))
    
        ! thermal elastic
        xs6l = 0D0
        if ( sab(isab_l)%jxs(4) /= 0 ) then
            call GET_IERG_SABE(isab_l,ierg,erg)
            ipfac = max(0D0,min(1D0,(erg-abe%erg(ierg)) &
                /(abe%erg(ierg+1)-abe%erg(ierg))))
            xs6l = abe%xs(ierg) + ipfac*(abe%xs(ierg+1)-abe%xs(ierg))
            if ( sab(isab_l)%nxs(5) == 4 ) xs6l = xs6l / abe%erg(ierg)
            if ( abe % erg ( ierg ) > erg ) xs6l = 0d0
        end if
        if ( associated(abi) ) nullify(abi)
        if ( associated(abe) ) nullify(abe)
    
        ! HIGH TEMP:
        ! thermal inelastic
        abi => sab(isab_h)%itie
        abe => sab(isab_h)%itce

        call GET_IERG_SABI(isab_h,ierg,erg)
        ipfac = max(0D0,min(1D0,(erg-abi%erg(ierg)) &
            /(abi%erg(ierg+1)-abi%erg(ierg))))
        xs2h = abi%xs(ierg) + ipfac*(abi%xs(ierg+1)-abi%xs(ierg))
    
        ! thermal elastic
        xs6h = 0D0
        if ( sab(isab_h)%jxs(4) /= 0 ) then
            call GET_IERG_SABE(isab_h,ierg,erg)
            ipfac = max(0D0,min(1D0,(erg-abe%erg(ierg)) &
                /(abe%erg(ierg+1)-abe%erg(ierg))))
            xs6h = abe%xs(ierg) + ipfac*(abe%xs(ierg+1)-abe%xs(ierg))
            if ( sab(isab_h)%nxs(5) == 4 ) xs6h = xs6h / abe%erg(ierg)
            if ( abe % erg ( ierg ) > erg ) xs6h = 0d0
        end if
        if ( associated(abi) ) nullify(abi)
        if ( associated(abe) ) nullify(abe)

        xs(2) = (1d0-f) * xs2l + f * xs2h
        xs(6) = (1d0-f) * xs6l + f * xs6h
        xs(2) = xs(2) + xs(6)  ! thermal scattering = inelastic + elastic
        xs(1) = xs(1) + xs(2)  ! total += thermal scattering
    
    endif

end subroutine




! =============================================================================
! GET_IERG_SABI
! =============================================================================
subroutine GET_IERG_SABI(iso_,ierg_,erg)
    integer, intent(in)::  iso_
    integer, intent(out):: ierg_
    real(8), intent(in)::  erg
    type(SAB_INEL_XS), pointer :: ab
    integer:: low, high, mid

    ab => sab(iso_)%itie
    
!    ! linear search
!    if ( erg > ab%erg(ab%ne) ) then
!        ierg_ = ab%ne-1
!    else
!    do ii = 2, ab%ne
!        if ( erg < ab%erg(ii) ) then
!            ierg_ = ii-1
!            exit
!        end if
!    end do
!    end if

    ! binary search
    low = 1
    high = ab%ne
    if ( erg > ab%erg(ab%ne) ) then
        ierg_ = ab%ne-1
    else
    do while ( low+1 /= high ) 
        mid = (low+high)/2
        if ( erg < ab%erg(mid) ) then
            high = mid
        else
            low = mid
        end if
    end do
    ierg_ = low
    end if

    if ( associated(ab) ) nullify(ab)

end subroutine

! =============================================================================
! GET_IERG_SABE
! =============================================================================
subroutine GET_IERG_SABE(iso_,ierg_,erg)
    integer, intent(in)::  iso_
    integer, intent(out):: ierg_
    real(8), intent(in)::  erg
    type(SAB_EL_XS), pointer:: ab
    integer:: low, mid, high

    ab => sab(iso_)%itce

    ! binary search
    low = 1
    high = ab%ne
    if ( erg > ab%erg(ab%ne) ) then
        ierg_ = ab%ne-1
    else
    do while ( low+1 /= high ) 
        mid = (low+high)/2
        if ( erg < ab%erg(mid) ) then
            high = mid
        else
            low = mid
        end if
    end do
    ierg_ = low
    end if

    if ( associated(ab) ) nullify(ab)

end subroutine


! =============================================================================
! GET_OTF_DB
! =============================================================================
subroutine GET_OTF_DB_MAC(nd, i_iso,iso,E0,xs1,dtemp)
    use ACE_HEADER, only: ace, ghq, wghq, ghq2, xghq2, wghq2
    use FMFD_HEADER, only: fmfdon
    use constants, only: k_b
    implicit none
    real(8), intent(in)   :: nd
    integer, intent(in)   :: i_iso, iso ! index for material and ACE
    real(8), intent(in)   :: E0
    real(8), intent(inout):: xs1(5)
    real(8), intent(in)   :: dtemp
    real(8):: xs0(5)
    real(8):: xn(2:4)
    real(8):: erg_l, erg_u, E1
    integer:: ierg0, ierg1
    real(8):: bb, yy, inv_b, inv_y, inv_y2  ! parameters 1
    real(8):: xx, x2, wx2  ! parameters 2
    real(8):: p1, p2       ! parameters 3
    integer:: ii

    ! parameters
    bb    = ace(iso)%atn/dtemp
    yy    = sqrt(bb*E0)
    inv_b = 1D0/bb
    inv_y = 1D0/yy
    inv_y2 = inv_y*inv_y

    ! initialization
    xs1(:) = 0D0

    if ( yy > ghq(16) ) then
        ! print *, 'CASE 1:', E0*1E6, bb
        erg_l = (yy+ghq(1))*(yy+ghq(1))*inv_b
        erg_u = (yy+ghq(16))*(yy+ghq(16))*inv_b
        call getierg(iso,ierg0,erg_l)
        call getierg(iso,ierg1,erg_u)

        do ii = 1, 16
            xx = ghq(ii) + yy
            x2 = xx*xx
            E1 = x2*inv_b
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB1(iso,ierg0,E1,xs0,xn(2:4))
            wx2 = wghq(ii) * x2
            xs1(:) = xs1(:) + wx2 * xs0(:) * (1d0-exp(-4d0*yy*xx))
        end do
        p1 = inv_sqrt_pi*inv_y2
        xs1(:) = xs1(:) * p1

    else
        erg_l = xghq2(1)*inv_b
        erg_u = xghq2(16)*inv_b
        call getierg(iso,ierg0,erg_l)
        call getierg(iso,ierg1,erg_u)

        ! print *, 'CASE 2:', inv_b, erg_l, erg_u
        do ii = 1, 16
            E1 = xghq2(ii)*inv_b
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB1(iso,ierg0,E1,xs0,xn(2:4))
            p1 = exp(2D0*ghq2(ii)*yy)
            p2 = wghq2(ii)*(p1-1D0/p1)
            xs1(:) = xs1(:) + p2 * xs0(:)
        end do
        p1 = inv_sqrt_pi*inv_y2*exp(-yy*yy)
        xs1(:) = xs1(:) * p1

    end if

    if ( fmfdon ) then
    do ii = 2, 4
        xs1(2) = xs1(2) - (dble(ii)-1D0)*xn(ii)
    end do
    end if

    xs1(:) = xs1(:) * nd * barn
end subroutine

! =============================================================================
! GET_OTF_DB_MT: MT-based OTF DB
! =============================================================================
subroutine GET_OTF_DB_MT(temp1,iso,E0,mt,xs1)
    use ACE_HEADER, only: ace, ghq, wghq, ghq2, xghq2, wghq2
    use CONSTANTS, only: k_b
    implicit none
    real(8), intent(in)   :: temp1  ! material tempearture where particle is
    integer, intent(in)   :: iso    ! index for MAT and ACE library
    real(8), intent(in)   :: E0
    integer, intent(in)   :: mt
    real(8), intent(inout):: xs1
    real(8):: xs0
    real(8):: erg_l, erg_u, E1
    integer:: ierg0, ierg1
    real(8):: bb, yy, inv_b, inv_y, inv_y2  ! parameters 1
    real(8):: xx, x2, wx2  ! parameters 2
    real(8):: p1, p2       ! parameters 3
    real(8):: nd
    integer:: ii
	
	
    ! parameters
    bb    = ace(iso)%atn/abs(temp1-ace(iso)%temp)
    yy    = sqrt(bb*E0)
    inv_b = 1D0/bb
    inv_y = 1D0/yy
    inv_y2 = inv_y*inv_y

    ! initialization
    xs1 = 0D0

    if ( yy > ghq(16) ) then
        erg_l = (yy+ghq(1))*(yy+ghq(1))*inv_b
        erg_u = (yy+ghq(16))*(yy+ghq(16))*inv_b
        call getierg(iso,ierg0,erg_l)
        call getierg(iso,ierg1,erg_u)

        do ii = 1, 16
            xx = ghq(ii) + yy
            x2 = xx*xx
            E1 = x2*inv_b
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB3(iso,ierg0,E1,mt,xs0)
            wx2 = wghq(ii) * x2
            xs1 = xs1 + wx2 * xs0 !* (1d0-exp(-4d0*yy*xx))
        end do

        p1 = inv_sqrt_pi*inv_y2
        xs1 = xs1 * p1

    else
        erg_l = xghq2(1)*inv_b
        erg_u = xghq2(16)*inv_b
        call getierg(iso,ierg0,erg_l)
        call getierg(iso,ierg1,erg_u)

        do ii = 1, 16
            E1 = xghq2(ii)*inv_b
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB3(iso,ierg0,E1,mt,xs0)
            p1 = exp(2D0*ghq2(ii)*yy)
            p2 = wghq2(ii)*(p1-1D0/p1)
            xs1 = xs1 + p2 * xs0
        end do

        p1 = inv_sqrt_pi*inv_y2*exp(-yy*yy)
        xs1 = xs1 * p1

    end if

end subroutine

! =============================================================================
! GET_OTF_DB
! =============================================================================
subroutine GET_OTF_DB_MIC(temp1,iso,E0,xs1)
    use ACE_HEADER, only: ace, ghq, wghq, ghq2, xghq2, wghq2
    use CONSTANTS, only: k_b
    implicit none
    real(8), intent(in)   :: temp1  ! material tempearture where particle is
    integer, intent(in)   :: iso    ! index for MAT and ACE library
    real(8), intent(in)   :: E0
    real(8), intent(inout):: xs1(6)
    real(8):: xs0(6)
    real(8):: erg_l, erg_u, E1
    integer:: ierg0, ierg1
    real(8):: bb, yy, inv_b, inv_y, inv_y2  ! parameters 1
    real(8):: xx, x2, wx2  ! parameters 2
    real(8):: p1, p2       ! parameters 3
    real(8):: nd
    integer:: ii
	
	
    ! parameters
    bb    = ace(iso)%atn/abs(temp1-ace(iso)%temp)
    yy    = sqrt(bb*E0)
    inv_b = 1D0/bb
    inv_y = 1D0/yy
    inv_y2 = inv_y*inv_y

    ! initialization
    xs1(:) = 0D0

    if ( yy > ghq(16) ) then
        erg_l = (yy+ghq(1))*(yy+ghq(1))*inv_b
        erg_u = (yy+ghq(16))*(yy+ghq(16))*inv_b
        call getierg(iso,ierg0,erg_l)
        call getierg(iso,ierg1,erg_u)

        do ii = 1, 16
            xx = ghq(ii) + yy
            x2 = xx*xx
            E1 = x2*inv_b
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB2(iso,ierg0,E1,xs0)
            wx2 = wghq(ii) * x2
            xs1(1:6) = xs1(1:6) + wx2 * xs0(1:6)! * (1d0-exp(-4d0*yy*xx))
        end do

        p1 = inv_sqrt_pi*inv_y2
        xs1(:) = xs1(:) * p1

    else
        erg_l = xghq2(1)*inv_b
        erg_u = xghq2(16)*inv_b
        call getierg(iso,ierg0,erg_l)
        call getierg(iso,ierg1,erg_u)

        do ii = 1, 16
            E1 = xghq2(ii)*inv_b
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB2(iso,ierg0,E1,xs0)
            p1 = exp(2D0*ghq2(ii)*yy)
            p2 = wghq2(ii)*(p1-1D0/p1)
            xs1(1:6) = xs1(1:6) + p2 * xs0(1:6)
        end do

        p1 = inv_sqrt_pi*inv_y2*exp(-yy*yy)
        xs1(:) = xs1(:) * p1

    end if

end subroutine

! =============================================================================
! GET_OTF_DB
! =============================================================================

function get_OTF_DB_UEG(mat, ierg0, dtemp) result (macro_xs)
    use ace_header, only: ace, ghq, wghq, ghq2, xghq2, wghq2
    use constants, only: k_b
    implicit none
    type(Material_CE), intent(in) :: mat
    integer, intent(in) :: ierg0
    real(8), intent(in) :: dtemp
    real(8) :: macro_xs(3) !Total, NuF, QF

end function

! =============================================================================
!
! =============================================================================
function EFF_IERG(E0,iso,p1,p2) result(pt4)
    use ACE_HEADER, only: ace
    implicit none
    real(8), intent(in):: E0
    integer, intent(in):: iso
    integer, intent(in):: p1, p2
    integer:: pt1, pt2, pt3, pt4

    pt1 = p1
    pt2 = p2

    !pt2 = pt2 + 1

    if ( pt1 == pt2 ) then
        pt4 = pt1 - 1
    else
        do
            if ( (pt2-pt1) == 1 ) exit
            pt3 = (pt2+pt1)/2
            if ( E0 >= ace(iso)%E(pt3) ) then
                pt1 = pt3
            else
                pt2 = pt3
            end if
        end do
        pt4 = pt1
    end if

end function

! =============================================================================
!
! =============================================================================
subroutine GET_MIC_DB1(iso,ierg,E1,xs,xn)
    use FMFD_HEADER, only: fmfdon
    implicit none
    integer, intent(in):: iso, ierg
    real(8), intent(in):: E1
    real(8):: xs(5)
    real(8):: xn(2:4)
    real(8):: slope
    integer:: ii, jj
    real(8):: xs_xn

    slope = max(0D0,min(1D0,(E1-ace(iso)%E(ierg)) &
        /(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))

    xs(1) = ace(iso)%sigt(ierg) &
          + slope * (ace(iso)%sigt(ierg+1)-ace(iso)%sigt(ierg))
    xs(2) = ace(iso)%sigd(ierg) &
          + slope * (ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
    xs(3:5) = 0D0

    ! fissionable material
    if ( ace(iso)%jxs(21) /= 0 .or. allocated(ace(iso)%sigf) ) then
    xs(3) = ace(iso)%sigf(ierg) &
          + slope * (ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
    xs(4) = xs(3) * getnu(iso,E1)
    xs(5) = xs(3) * ace(iso)%qval
    xs(2) = xs(2) + xs(3)
    end if

    ! (n,xn) cross-section
    ! seperately considered and then collapsed? when FMFD is on
    xn(:) = 0D0
    if ( fmfdon ) then
    do ii = 1, ace(iso)%nxs(5)
        jj = abs(ace(iso)%TY(ii))
        if ( jj > 1 .and. jj < 5 ) then
            xn(jj) = xn(jj) + ace(iso)%sig_MT(ii)%cx(ierg) + slope * &
                (ace(iso)%sig_MT(ii)%cx(ierg+1)-ace(iso)%sig_MT(ii)%cx(ierg))
        end if
    end do
    end if

end subroutine

! =============================================================================
!
! =============================================================================
subroutine GET_MIC_DB2(iso,ierg,E1,xs)
    use FMFD_HEADER, only: fmfdon
    implicit none
    integer, intent(in):: iso, ierg
    real(8), intent(in):: E1
    real(8):: xs(6)
    real(8):: slope
    integer:: ii, jj
    real(8):: xs_xn(4)
    ! 1 : total
    ! 2 : elastic scattering
    ! 3 : absorption
    ! 4 : fission
    ! 5 : nu-fission
    ! 6 : disapperance

    slope = max(0D0,min(1D0,(E1-ace(iso)%E(ierg)) &
        /(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))

    xs(1) = ace(iso)%sigt(ierg) &
          + slope * (ace(iso)%sigt(ierg+1)-ace(iso)%sigt(ierg))
    xs(2) = ace(iso)%sigel(ierg) &
          + slope * (ace(iso)%sigel(ierg+1)-ace(iso)%sigel(ierg))
    xs(3) = ace(iso)%sigd(ierg) &
          + slope * (ace(iso)%sigd(ierg+1)-ace(iso)%sigd(ierg))
    xs(4) = 0D0
    xs(5) = 0D0
    xs(6) = xs(3)

    ! fissionable material
    if ( ace(iso)%jxs(21) /= 0 .or. allocated(ace(iso)%sigf) ) then
    xs(4) = ace(iso)%sigf(ierg) &
          + slope * (ace(iso)%sigf(ierg+1)-ace(iso)%sigf(ierg))
    xs(5) = xs(4)*getnu(iso,E1)
    xs(3) = xs(3) + xs(4)
    end if

end subroutine
! =============================================================================
!
! =============================================================================
subroutine GET_MIC_DB3(iso,ierg,E1,mt,xs)
    use FMFD_HEADER, only: fmfdon
    implicit none
    integer, intent(in):: iso, ierg, mt
    real(8), intent(in):: E1
    real(8):: xs
    real(8):: slope

    slope = max(0D0,min(1D0,(E1-ace(iso)%E(ierg)) &
        /(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))

    xs = ace(iso) % sig_MT(mt) % cx(ierg) + &
        slope * (ace(iso) % sig_MT(mt) % cx(ierg+1) - &
        ace(iso) % sig_MT(mt) % cx(ierg))

end subroutine

subroutine URES_PTABLE(iso, erg, xs, IFF, r)
    use RANDOMS, only: rang
    integer, intent(in) :: iso
    real(8), intent(in) :: erg
    real(8), intent(in) :: r
    integer, intent(out):: IFF

    real(8) :: xs(4)
    real(8) :: xs0(4), xs1(4)
    real(8) :: xsarr0(4), xsarr1(4)
    integer :: col0, col1
    type(UNRtype), pointer :: ac
    real(8), pointer :: C0(:), C1(:)
    integer :: N, M, ierg, mm, ii, j1, j2, iMT
    integer :: iergu, ierg0, ierg1
    real(8) :: ipfac, ipfac0, ipfac1
    real(8) :: f
    real(8) :: xst, xse, xsie, xsa, xsf, xsd

    ac => ace(iso) % UNR
    N = ac % N; M = ac % M
    IFF = ac % IFF
    
    ! 1-1. Select PTABLE with erg
    do iergu = 1, N-1
        if(ac % E(iergu) <= erg .and. erg < ac % E(iergu+1)) exit
    enddo
    call getiueg(ac%E(iergu), j1)
    call getiueg(ac%E(iergu+1), j2)

    !==== NOTES ====
    ! 1: CDF
    ! 2: Total XS
    ! 3: Elastic XS
    ! 4: Fission XS
    ! 5: (n,gamma) XS

    ! 2. Call Random to choose column
    C0 => ac % P(iergu  ,1,:)

    col0 = 1; col1 = 1 ! To avoid crash...
    do mm = 1, M
        ! ACE format implies that CDF monotonically increases
        if(r < C0(mm)) then
            col0 = mm
            exit
        endif
    enddo

    
    !col0 = nint(r)
    !col0 = r; col1 = r

    !print *, 'T', erg, r, col0
    
    xs0 = ac % P(iergu  ,2:5,col0)
    xs1 = ac % P(iergu+1,2:5,col0)

    ! 5. Interpolate Probability Table
    select case(ac % INT)
    case(2) ! Linear
        f  = (erg-ac % E(iergu))/(ac % E(iergu+1) - ac % E(iergu))
        xs = xs0 *  (1D0-f) + xs1  * f
    case(5) ! Log-Log
        f  = (log(erg/ac % E(iergu))) / (log(ac%E(iergu+1)/ac%E(iergu)))
        !xs = exp(log(xs0) + f * log(xs1/xs0))
        xs(2) = exp(log(xs0(2)) + f * log(xs1(2)/xs0(2)))
        xs(3) = exp(log(xs0(3)) + f * log(xs1(3)/xs0(3)))
        xs(4) = exp(log(xs0(4)) + f * log(xs1(4)/xs0(4)))
    case default
        write(*,*) 'INVALID INTERPOLATION...', ac % INT
    end select

    call getierg(iso, ierg, erg)
    ipfac = max(0D0,min(1D0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))

    if(ac % IFF > 0) then
        xst = 0d0; xse = 0d0; xsa = 0d0; xsf = 0d0; xsd = 0d0
        xst  = xst  + (1d0-ipfac) * ace(iso) %  sigt(ierg) + ipfac * ace(iso) %  sigt(ierg+1)
        if(allocated(ace(iso)%sigel)) &
            xse  = xse  + (1d0-ipfac) * ace(iso) % sigel(ierg) + ipfac * ace(iso) % sigel(ierg+1)
        if(allocated(ace(iso)%sigf)) &
            xsf  = xsf  + (1d0-ipfac) * ace(iso) %  sigf(ierg) + ipfac * ace(iso) %  sigf(ierg+1)
        if(allocated(ace(iso)%sigel)) &    
            xsd  = xsd  + (1d0-ipfac) * ace(iso) %  sigd(ierg) + ipfac * ace(iso) %  sigd(ierg+1)
        xs(2) = xs(2) * xse
        xs(3) = xs(3) * xsf
        xs(4) = xs(4) * xsd
    endif
    xs(1) = sum(xs(2:4))
    ! 
    ! 4. Additional XS from contribution factors; ILF and IOA
    if(ac % ILF > 0) then
        xs(1) = xs(1) + &
            ace(iso) % sig_MT(ac % ILFidx) % cx(ierg) * (1D0-ipfac) + & 
            ace(iso) % sig_MT(ac % ILFidx) % cx(ierg+1) * ipfac
    elseif(ac % ILF == 0) then ! USE BALANCE: ILF = TOT-ABS-ELAS
        xs(1) = xs(1) + &
            xst - xsf - xsd - xse
    endif
    
    if(ac % IOA > 0) then
        xs(1) = xs(1) + &
            ace(iso) % sig_MT(ac % IOAidx) % cx(ierg) * (1D0-ipfac) + &
            ace(iso) % sig_MT(ac % IOAidx) % cx(ierg+1) * ipfac
        xs(4) = xs(4) + &
            ace(iso) % sig_MT(ac % IOAidx) % cx(ierg) * (1d0-ipfac) + &
            ace(iso) % sig_MT(ac % IOAidx) % cx(ierg+1) * ipfac
    elseif(ac % IOA == 0) then
    endif
    

    if(associated(C0)) nullify(C0)
    if(associated(ac)) nullify(ac)

end subroutine

subroutine GET_URR_MICRO(iso, erg, xs, urn)
    integer, intent(in) :: iso
    real(8), intent(in) :: erg
    real(8), intent(in) :: urn(:)
    real(8), intent(inout) :: xs(:)
    real(8) :: xs0(4)
    integer :: iff, ierg, idx
    real(8) :: ipfac
    ! 1: Total, 2: Elastic, 3: Abs, 4: Fiss, 5: NuFis, 6: (n,gamma)
    if(.not. ace(iso) % UNR % URES ) return
    if(erg < ace(iso) % UNR % Emin .or. erg > ace(iso) % UNR % Emax) return

    call URES_PTABLE(iso, erg, xs0, iff, urn(ace(iso) % UNR % unridx))

    xs(1) = xs0(1)
    xs(2) = xs0(2)
    xs(4) = xs0(3)
    xs(5) = xs0(3) * getnu(iso,erg)
    xs(6) = xs0(4)
    xs(3) = xs(4) + xs(6)

    !xs(1) = xs(2) + xs(3) + xs(4)

end subroutine

recursive subroutine quicksort(arr, first, last)
  implicit none
  real(8) :: arr(:), x, tmp
  integer :: first, last
  integer :: pt1, pt2

  x = arr( (first+last) / 2 )
  pt1 = first;  pt2 = last
  do
     do while (arr(pt1) < x)
        pt1 = pt1 + 1
     end do
     do while (x < arr(pt2))
        pt2 = pt2 - 1
     end do
     if (pt1 >= pt2) exit
     tmp = arr(pt1)
     arr(pt1) = arr(pt2)
     arr(pt2) = tmp
     pt1 = pt1 + 1
     pt2 = pt2 - 1
  end do
  if (first < pt1 - 1) call quicksort(arr, first, pt1 - 1)
  if (pt2 + 1 < last)  call quicksort(arr, pt2 + 1, last)
end subroutine quicksort

subroutine setueg
use constants 
integer :: iso, i, r, ierg, idx
real(8) :: ipfac
type(AceFormat), pointer :: ac
    !$OMP PARALLEL DO PRIVATE(i, r, iso, ac, idx, ipfac)
    do iso = 1, num_iso
        ! 1. Initialize Egrid: corresponding E points for each iso.
        ac => ace(iso)
        allocate(ac%UEG%Egrid(1:nueg))
        !if(icore==score) print *, iso, nueg
        !   * Set Egrid: Closest lower E point of iso for each ueggrid
        ac % UEG % Egrid = ac % NXS(3)+1
        idx = 0
        do i = 1, nueg
            !if(iso==1 .and. ueggrid(i) >= ac % E(idx+1)) print *, ac % zaid, idx, ueggrid(i), ac % E(idx), ac % E(idx+1)
            if(ueggrid(i) >= ac % E(idx+1)) idx = idx + 1
            ac % UEG % Egrid(i) = idx
            if(idx >= ac % NXS(3)) exit
        enddo

        allocate(ac%UEG%sigt(1:nueg))
        ac % UEG % sigt = (ac % sigt ( ac % NXS(3) ))
        do i = 1, nueg
            if( ac % UEG % Egrid(i) == 0 ) then ! For Non-defined XS...
                ac % UEG % sigt(i) = (ac % sigt(1))
                cycle
            elseif ( ac % UEG % Egrid(i) == ac % NXS(3) ) then ! For Upper
                exit
            endif
            ipfac = max(0D0,min(1D0,(ueggrid(i)-ac % E ( ac % UEG % Egrid(i)) )/( ac % E(ac % UEG % Egrid(i)+1) - ac % E(ac % UEG % Egrid(i)) )))
            ac % UEG % sigt(i) = (ac % sigt( ac % UEG % Egrid(i)) + (ac % sigt( ac % UEG % Egrid(i)+1 ) - ac % sigt(ac % UEG % Egrid(i) )) * ipfac)
        enddo
        !if(icore==score) print *, 'Total XS for ', ace(iso) % xslib
            
        allocate(ac%UEG%sigd(1:nueg))
        ac % UEG % sigd = (ac % sigd ( ac % NXS(3) ))
        do i = 1, nueg-1
            if( ac % UEG % Egrid(i) == 0 ) then ! For Non-defined XS...
                ac % UEG % sigd(i) = (ac % sigd(1))
                cycle
            elseif ( ac % UEG % Egrid(i) == ac % NXS(3) ) then ! For Upper
                exit
            endif
            ipfac = max(0D0,min(1D0,(ueggrid(i)-ac % E ( ac % UEG % Egrid(i)) )/( ac % E(ac % UEG % Egrid(i)+1) - ac % E(ac % UEG % Egrid(i)) )))
            ac % UEG % sigd(i) = ( ac % sigd( ac % UEG % Egrid(i)) + (ac % sigd( ac % UEG % Egrid(i)+1 ) - ac % sigd(ac % UEG % Egrid(i) )) * ipfac)
        enddo
        !if(icore==score) print *, '(n,gamma) XS for ', ace(iso) % xslib
        
        if(allocated(ace(iso)%sigf)) then
        allocate(ac%UEG%sigf(1:nueg))
        allocate(ac%UEG%signuf(1:nueg))
!        allocate(ac%UEG%sigqf(1:nueg))
        ! 1) Fission XS
        ac % UEG % sigf = (ac % sigf ( ac % NXS(3) ))
        do i = 1, nueg-1
            if( ac % UEG % Egrid(i) == 0 ) then ! For Non-defined XS...
                ac % UEG % sigf(i) = (ac % sigf(1))
                cycle
            elseif ( ac % UEG % Egrid(i) == ac % NXS(3) ) then ! For Upper
                exit
            endif
            ipfac = max(0D0,min(1D0,(ueggrid(i)-ac % E ( ac % UEG % Egrid(i)) )/( ac % E(ac % UEG % Egrid(i)+1) - ac % E(ac % UEG % Egrid(i)) )))
            ac % UEG % sigf(i) = (ac % sigf( ac % UEG % Egrid(i)) + (ac % sigf( ac % UEG % Egrid(i)+1 ) - ac % sigf(ac % UEG % Egrid(i) )) * ipfac)
        enddo

!        ! 2) NuFiss XS and Q Fiss XS
        do i = 1, nueg
            ac % UEG % signuf(i) = ac % UEG % sigf (i) * getnu(iso, ueggrid(i))
        enddo
        endif
        
        if(allocated(ac % sigel)) then
        allocate(ac%UEG%sigel(1:nueg))
        ac % UEG % sigel = (ac % sigel( ac % NXS(3) ))
        do i = 1, nueg-1
            if( ac % UEG % Egrid(i) == 0 ) then ! For Non-defined XS...
                ac % UEG % sigel(i) = (ac % sigel(1))
                cycle
            elseif ( ac % UEG % Egrid(i) == ac % NXS(3) ) then ! For Upper
                exit
            endif
            ipfac = max(0D0,min(1D0,(ueggrid(i)-ac % E ( ac % UEG % Egrid(i)) )/( ac % E(ac % UEG % Egrid(i)+1) - ac % E(ac % UEG % Egrid(i)) )))
            ac % UEG % sigel(i) = (ac % sigel( ac % UEG % Egrid(i)) + (ac % sigel( ac % UEG % Egrid(i)+1 ) - ac % sigel(ac % UEG % Egrid(i) )) * ipfac)
        enddo
        endif
        if(icore==score) write(*,'(A,A,I4,A,I4)') '   UEG XS set for ', trim(ace(iso) % xslib), iso, '/', num_iso
    enddo
    !$OMP END PARALLEL DO
end subroutine

subroutine setures
    ! Additional Cross-section
    ! For URES isotopes, exclude from setMacro
    ! Instead, add ac % UEG % sigX_ures(:,:)
    ! sig_ures: nueg x M 
    ! X: sigt, sigel, siga, sigf, signuf
    integer :: iso, i, j, ierg
    real(8) :: erg, ipfac
    
    type(AceFormat), pointer :: ac

    ! Currently, using Linear interpolation
    ! NEEDS TO CONSIDER IFF
    if(.not. do_ures) return
    allocate(uresiso(1:n_unr));
    !$OMP PARALLEL DO PRIVATE(iso, ac, erg, j, ierg, ipfac)
    do iso = 1, num_iso
        ! 1. Initialize sig_ures
        ac => ace(iso)
        if(.not. ac % UNR % URES ) cycle
        uresiso(ac % UNR % unridx) = iso
    enddo

end subroutine

subroutine setMacroXS(BU)
    ! OBJECTIVE: Obtain UEG-wise XS
    ! 1 >> Total XS
    ! 2 >> Absorption XS
    ! 3 >> Fission XS
    ! 4 >> NuFiss. XS
    ! 5 >> QFiss. XS
    type(Material_CE), pointer :: mat
    logical, intent(in) :: BU
    integer :: i, j, iso, ii, imat
    integer :: nprod
    real(8) :: micro(5), xs(6), ipfac, xs1(5)

    integer :: ierg, cnt
    integer :: ierg_sab, ierg_otf, epoint, i_low, i_high
    real(8) :: ff

    cnt = 0
    !$OMP PARALLEL DO PRIVATE(iso, mat, micro, nprod, xs, xs1, ierg_sab, ierg_otf, i, j, epoint)
    do imat = 1, n_materials
    !if(BU .and. .not. materials(imat)%depletable) cycle
    mat => materials(imat)
    if(.not. allocated(mat % macro_ueg)) allocate(mat % macro_ueg(1:nueg, 5))
    
    if(do_ures .and. n_unr>0 .and. .not. allocated(mat % ures)) &
        allocate(mat % ures(1:n_unr))
    mat % ures = .false.
    if(do_ures .and. n_unr>0 .and. .not. allocated(mat % uresidx)) &
        allocate(mat % uresidx(1:n_unr))
    mat % uresidx = 0
    
    mat % macro_ueg = 0D0
    do i = 1, mat % n_iso
        iso = mat % ace_idx(i)
        if(ace(iso) % UNR % URES .and. mat % numden(i) > ures_cut) then
            mat % ures(ace(iso) % UNR % unridx) = .true.
            mat % uresidx(ace(iso)%UNR% unridx) = i
        endif

        if(do_ures) then
            do j = 1, nueg
                if(ace(iso) % UNR % URES) then
                !if(mat % ures(ace(iso) % UNR % unridx)) cycle
                if(mat % ures(ace(iso) % UNR % unridx) .and. &
                    ueggrid(j) >= ace(iso) % UNR % Emin .and. &
                    ueggrid(j) <= ace(iso) % UNR % Emax) cycle
                endif
                micro    = 0D0
                micro(1) = (ace(iso) % UEG % sigt(j))
                if(allocated(ace(iso)%UEG%sigd)) &
                    micro(2) = (ace(iso) % UEG % sigd(j))
    
                if(allocated(ace(iso)%UEG%sigf)) then
                    micro(2) = micro(2) + ace(iso)%UEG%sigf(j)
                    micro(3) = (ace(iso) % UEG % sigf(j))
                    micro(4) = (getnu(iso,ueggrid(j)) * ace(iso) % UEG % sigf(j))
                    micro(5) = (ace(iso) % qval * ace(iso) % UEG % sigf(j))
                    !micro(2) = (micro(2) + micro(3))
                endif
    
                mat % macro_ueg(j,:) = mat % macro_ueg(j,:) &
                   + (micro(:) * mat % numden(i) * barn)
            enddo

        elseif ( mat % sablist(i) > 0 ) then ! S(a,b), no interp.
            call getiueg( 4d-6, ierg_sab )
            do epoint = 1, ierg_sab
            call GET_SAB_MAC(mat % numden(i), iso, mat%sablist(i), ueggrid(epoint), &
                mat % macro_ueg(epoint,1), mat % macro_ueg(epoint,2))
            enddo
            mat % macro_ueg(ierg_sab + 1:, 1) =  mat % macro_ueg(ierg_sab + 1:, 1) + &
                ace(iso) % UEG % sigt(ierg_sab+1:) * mat % numden(i) * barn
            mat % macro_ueg(ierg_sab + 1:, 2) =  mat % macro_ueg(ierg_sab + 1:, 2) + &
                ace(iso) % UEG % sigd(ierg_sab+1:) * mat % numden(i) * barn

        elseif ( mat % sablist(i) < 0 ) then ! S(a,b), interp.
            call getiueg( 4d-6, ierg_sab )
            i_low = therm(-mat%sablist(i)) % iso_low
            i_high= therm(-mat%sablist(i)) % iso_high
            ff    = therm(-mat%sablist(i)) % f
            do epoint = 1, ierg_sab
                call GET_SAB_MAC(mat % numden(i) * (1d0-ff), iso, i_low, ueggrid(epoint), &
                    mat % macro_ueg(epoint,1), mat % macro_ueg(epoint,2))
                call GET_SAB_MAC(mat % numden(i) * ff, iso, i_high, ueggrid(epoint), &
                    mat % macro_ueg(epoint,1), mat % macro_ueg(epoint,2))
            enddo
            mat % macro_ueg(ierg_sab + 1:, 1) =  mat % macro_ueg(ierg_sab + 1:, 1) + &
                ace(iso) % UEG % sigt(ierg_sab+1:) * mat % numden(i) * barn
            mat % macro_ueg(ierg_sab + 1:, 2) =  mat % macro_ueg(ierg_sab + 1:, 2) + &
                ace(iso) % UEG % sigd(ierg_sab+1:) * mat % numden(i) * barn

        elseif ( mat % db .and. ((mat%temp-ace(iso)%temp) > K_B * 1e-2 ) ) then
            call getiueg( 1d0, ierg_otf )
            ierg_otf = ierg_otf - 1
            do epoint = 1, ierg_otf
                call GET_OTF_DB_MAC( mat % numden(i), i, iso, ueggrid(epoint), xs1, (mat%temp-ace(iso)%temp))
                mat % macro_ueg(epoint, 1) = mat % macro_ueg(epoint, 1) + &
                    xs1(1)
                mat % macro_ueg(epoint, 2) = mat % macro_ueg(epoint, 2) + &
                    xs1(2)
                mat % macro_ueg(epoint, 3) = mat % macro_ueg(epoint, 3) + &
                    xs1(3)
                mat % macro_ueg(epoint, 4) = mat % macro_ueg(epoint, 4) + &
                    xs1(4)
                mat % macro_ueg(epoint, 5) = mat % macro_ueg(epoint, 5) + &
                    xs1(5)
            enddo
            mat % macro_ueg(ierg_otf+1:,1) = mat % macro_ueg(ierg_otf+1:,1) + &
                ace(iso) % UEG % sigt(ierg_otf+1:) * mat % numden(i) * barn

            mat % macro_ueg(ierg_otf+1:,2) = mat % macro_ueg(ierg_otf+1:,2) + &
                ace(iso) % UEG % sigd(ierg_otf+1:) * mat % numden(i) * barn

            if(allocated(ace(iso) % UEG % sigf)) then
                mat % macro_ueg(ierg_otf+1:,3) = mat % macro_ueg(ierg_otf+1:,3) + &
                    ace(iso) % UEG % sigf(ierg_otf+1:) * mat % numden(i) * barn
                mat % macro_ueg(ierg_otf+1:,4) = mat % macro_ueg(ierg_otf+1:,4) + &
                    ace(iso) % UEG % signuf(ierg_otf+1:) * mat % numden(i) * barn
                mat % macro_ueg(ierg_otf+1:,5) = mat % macro_ueg(ierg_otf+1:,5) + &
                    ace(iso) % UEG % sigf(ierg_otf+1:) * mat % numden(i) * barn * ace(iso) % qval
            endif               
        else
            mat % macro_ueg(:,1) = mat % macro_ueg(:,1) + &
                ace(iso) % UEG % sigt(:) * mat % numden(i) * barn

            mat % macro_ueg(:,2) = mat % macro_ueg(:,2) + &
                ace(iso) % UEG % sigd(:) * mat % numden(i) * barn

            if(allocated(ace(iso) % UEG % sigf)) then
                mat % macro_ueg(:,3) = mat % macro_ueg(:,3) + &
                    ace(iso) % UEG % sigf(:) * mat % numden(i) * barn
                mat % macro_ueg(:,4) = mat % macro_ueg(:,4) + &
                    ace(iso) % UEG % signuf(:) * mat % numden(i) * barn
                mat % macro_ueg(:,5) = mat % macro_ueg(:,5) + &
                    ace(iso) % UEG % sigf(:) * mat % numden(i) * barn * ace(iso) % qval
            endif               
        endif
    enddo
    if(n_materials < 10) then
        if(icore==score) write(*,*) 'XS SET FOR MAT #:', imat
    else ! n_materials >= 100
        !$OMP ATOMIC
        cnt = cnt + 1
        if(icore/=score) cycle
        if(cnt == n_materials/4) then
            print *, 'MATERIAL XS CALC ( 25%/100%)'
        elseif(cnt == n_materials/4*2) then
            print *, 'MATERIAL XS CALC ( 50%/100%)'
        elseif(cnt == n_materials/4*3) then
            print *, 'MATERIAL XS CALC ( 75%/100%)'
        elseif(cnt == n_materials) then
            print *, 'MATERIAL XS CALC (100%/100%)'
        endif
    endif
    enddo

end subroutine

subroutine setDBPP(BU)
implicit none
integer, intent(in) :: BU
integer :: i, j, iso, iso_, rx
real(8) :: xs1(6), urn(n_unr), xs2
type(AceFormat), pointer :: ac
logical :: found
real(8) :: base_tmp
! 23/12/04 : Preprocessor
if (.not. allocated(ace_base)) then
    if(icore==score) print *, 'Base ACE format is not allocated'
    return
endif

do i = 1, n_materials
    if ( materials(i) % db ) cycle
        do iso = 1, materials(i) % n_iso
            if( abs(materials(i) % temp - ace(materials(i)%ace_idx(iso)) % temp) > 1E-3*K_B ) then 
                found = .false.
                ISO_LOOP: do iso_ = 1, num_iso
                    if( abs(materials(i) % temp - ace(iso_) % temp) < 1E-3 * K_B .and. &
                        ace(materials(i)%ace_idx(iso)) % zaid == ace(iso_) % zaid) then
                        materials(i) % ace_idx(iso) = iso_
                        found = .true.
                    exit ISO_LOOP
                    endif
                enddo ISO_LOOP
                    
                if( .not. found ) then
                    num_iso = num_iso + 1
                    ace(num_iso) = ace(materials(i)%ace_idx(iso))
                    ac => ace(num_iso)
                    ac % temp = materials(i) % temp
                    if(icore==score .and. ac % zaid == 92238) print *, 'U238 TST', ac % E(1), ac % E(ac % NXS(3)), ac % NXS(3), allocated(ugrid)
                    do j = 1, ac % NXS(3)
                        if ( ac % E(j) > 1d0 ) exit
                        call GET_OTF_DB_MIC(materials(i)%temp, materials(i)%ace_idx(iso), ac % E(j), xs1)
                        ac % sigt(j) = xs1(1)
                        ac % sigel(j) = xs1(2)
                        ac % sigd(j) = xs1(3)-xs1(4)
                        if( allocated (ac % sigf) ) ac % sigf(j) = xs1(4)
                        do rx = 1, ac % NXS(5)
                            call GET_OTF_DB_MT(materials(i)%temp, materials(i)%ace_idx(iso), ac % E(j), rx, xs2)
                            ac % sig_MT(rx) % cx(j) = xs2
                        enddo
                    end do
                    nullify(ac)
                    materials(i) % ace_idx(iso) = num_iso
                    if(icore==score) print *, trim(materials(i)%mat_name), ': Adjusted XS for ', trim(ace(num_iso) % xslib), ' to', ace(num_iso) % temp/K_B
                else
                    if(icore==score) print *, trim(materials(i)%mat_name), ': Linked XS to ', trim(ace(materials(i)%ace_idx(iso)) % xslib), ' with T:', ace(materials(i)%ace_idx(iso)) % temp / K_B
                endif
            elseif( abs(materials(i) % temp - ace(materials(i)% ace_idx(iso)) % temp) > 1E-3 * K_B) then
                if(icore==score) print *, 'WARNING: Invalid Temperature for ', trim(materials(i)%mat_name), materials(i)%temp/K_B, ace(materials(i)%ace_idx(iso))%temp/K_B
            else
                if(icore==score) print *, trim(materials(i)%mat_name), ': no adjust required for ', trim(ace(materials(i)%ace_idx(iso))%xslib)
            endif
        enddo
enddo
end subroutine

end module
