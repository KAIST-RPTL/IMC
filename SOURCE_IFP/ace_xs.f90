module ace_xs

use constants, only : barn
use variables, only : E_mode, k_steady
use material_header 
use ace_header 
use ace_module 

implicit none
    real(8), parameter:: inv_sqrt_pi = 0.564189583547756D0 ! 1/sqrt(pi)

contains

function getMacroXS (mat, erg,kT,urn) result (macro_xs)
    use CONSTANTS, only: K_B
    implicit none
    type(Material_CE), intent(in) :: mat
    real(8), intent(in) :: erg
    real(8), intent(in) :: kT
    real(8) :: macro_xs(5), urn(1:n_unr)
    
    integer :: i, iff
    integer :: i_iso, iso_, ierg_
    integer :: pt1, pt2, pt3, pt4
    real(8) :: ipfac
    real(8) :: micro_t, micro_d, micro_f, micro_nuf, micro_a, micro_el, micro_xn
    real(8) :: macro_t, macro_f, macro_nuf, macro_a, macro_qf
    real(8) :: xn_xs(4)
    real(8) :: xs(5)
    integer :: isab
    real(8) :: dtemp    ! temperautre difference | library - material |
    macro_t   = 0.0d0
    macro_a   = 0.0d0
    macro_f   = 0.0d0
    macro_nuf = 0.0d0
    macro_qf  = 0.0d0
    xn_xs(:)  = 0.0d0

    !print *, mat%mat_name, mat%n_iso
    do i_iso = 1, mat%n_iso     ! isotope number in the material
    
        iso_ = mat%ace_idx(i_iso)   ! isotope number in the inputfile

        ! =====================================================================
        ! S(a,b) treatment
        isab = ace(iso_)%sab_iso    ! isotope number for S(a,b)
        if ( mat%sab .and. isab /= 0 .and. erg < 4D-6 ) then
            call GET_SAB_MAC(mat%numden(i_iso),iso_,isab,erg,macro_t,macro_a)
            cycle
        end if

        ! =====================================================================
        ! On-the-fly Doppler broadening
        dtemp = abs(ace(iso_)%temp-kT)
        if ( mat%db .and. ( dtemp > K_B .and. erg < 1d0 ) ) then
            call GET_OTF_DB_MAC(mat%numden(i_iso),i_iso,iso_,erg,xs,dtemp)
            macro_t   = macro_t   + xs(1)
            macro_a   = macro_a   + xs(2)
            macro_f   = macro_f   + xs(3)
            macro_nuf = macro_nuf + xs(4)
            macro_qf  = macro_qf  + xs(5)
            cycle
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
        if  (ace(iso_) % UNR % URES) then
            ! Only for Energy between Emin and Emax
            if ( erg >= ace(iso_) % UNR % Emin .and. &
                    erg <= ace(iso_) % UNR % Emax ) then
                call URES_PTABLE(iso_, erg, xs, iff, urn(ace(iso_) % UNR % unridx))
                
!                if(iff == 1) then ! IFF = 1: Ratio
!                    micro_t   = micro_t   * xs(1)
!                    micro_f   = micro_f   * xs(3)
!                    micro_nuf = micro_nuf * xs(3)
!                    micro_d   = micro_d   * xs(4)
!                    micro_a   = micro_d + micro_f
!                else ! IFF = 0: Absolute
                    micro_t   = xs(1)
                    micro_f   = xs(3)
                    micro_nuf = micro_f * getnu(iso_,erg)
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
                micro_xn   = ace(iso_)%sig_MT(i)%cx(ierg_) & 
                            + ipfac*(ace(iso_)%sig_MT(i)%cx(ierg_+1) - ace(iso_)%sig_MT(i)%cx(ierg_))
                xn_xs(pt1) = xn_xs(pt1) + mat%numden(i_iso) * micro_xn * barn
            endif
        enddo

    enddo 

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
        if(erg0 >= ace(iso_)%E(pt3)) then
          pt1 = pt3
        else
          pt2 = pt3
        endif
      Enddo
    endif
    ierg_ = pt1 !store low bound energy index 
    
    !if(erg0 < 1.d-11) print *, erg0, ierg_ !, ace(iso_)%E(ierg_), ace(iso_)%E(ierg_+1)
end subroutine 


subroutine setugrid
    implicit none
    real(8) :: Etmp
    integer :: i, j, k, iso_, idx

    if(E_mode==0) return
    !Set ugrid to accelerate energy-grid search
    Emin = 1.d-11
    allocate(ugrid(0:nugrid,1:num_iso))
    !print *, 'ugrid size', nugrid, num_iso
    do iso_ = 1, num_iso
      ugrid(0,iso_) = 1
      ugrid(nugrid,iso_) = ace(iso_)%nxs(3)
    end do
    udelta = log10(Emax/Emin)/dble(nugrid)

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


!> Only if OTF-DB is Enabled

!    if(num_iso0K) then
!      allocate(ugrid0K(0:nugrid,1:num_iso0K))
!      do iso_ = 1, num_iso0K
!        ugrid0K(0,iso_) = 1
!        ugrid0K(nugrid,iso_) = ace0K(iso_)%nxs(3)
!      end do
!      
!      do iso_ = 1, num_iso0K
!        idx = 1
!        do i=1, nugrid-1
!          Etmp = Emin*10.d0**(dble(i)*udelta)
!          if(Etmp > ace0K(iso_)%xss(ace0K(iso_)%nxs(3))) then
!            idx = ace0K(iso_)%nxs(3)
!            go to 20  
!          end if
!          do
!            if(Etmp < ace0K(iso_)%xss(idx)) go to 20
!            idx = idx + 1
!          end do
!20        ugrid0K(i,iso_) = idx - 1
!        end do
!      end do
!    end if


    if(icore==score) print *, "   Setting ugrid..."
    
    
    !do i=1, nugrid-1
    !    write(99, *) i, Emin*10.d0**(dble(i)*udelta)
    !enddo 
    !stop
    
end subroutine 


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
    end if

    micro_t = micro_i + micro_e + micro_a

    xs_t = xs_t + nd * micro_t * barn
    xs_a = xs_a + nd * micro_a * barn

    if ( associated(abi) ) nullify(abi)
    if ( associated(abe) ) nullify(abe)

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
    integer:: iiso, isab, ierg
    real(8):: ipfac

    if ( .not. mat%sab .or. erg > 4E-6 ) return
    iiso = mat%ace_idx(imat)
    isab = ace(iiso)%sab_iso
    if ( isab == 0 ) return

    abi => sab(isab)%itie
    abe => sab(isab)%itce

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
    end if

    xs(2) = xs(2) + xs(6)  ! thermal scattering = inelastic + elastic
    xs(1) = xs(1) + xs(2)  ! total += thermal scattering

    if ( associated(abi) ) nullify(abi)
    if ( associated(abe) ) nullify(abe)

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
subroutine GET_OTF_DB_MAC(nd,i_iso,iso,E0,xs1,dtemp)
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
            xs1(:) = xs1(:) + wx2 * xs0(:)
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
!    real(8):: ee0(16), ee1(16), ee2(16)
	
	
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
!            ee0(ii) = xx
            x2 = xx*xx
            E1 = x2*inv_b
!            ee1(ii) = E1
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB2(iso,ierg0,E1,xs0)
            wx2 = wghq(ii) * x2
            xs1(1:6) = xs1(1:6) + wx2 * xs0(1:6)
!            ee2(ii) = xs0(1)
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
!            ee0(ii) = sqrt(xghq2(ii))
!            ee1(ii) = E1
            ierg0 = EFF_IERG(E1,iso,ierg0,ierg1+1)
            call GET_MIC_DB2(iso,ierg0,E1,xs0)
            p1 = exp(2D0*ghq2(ii)*yy)
            p2 = wghq2(ii)*(p1-1D0/p1)
            xs1(1:6) = xs1(1:6) + p2 * xs0(1:6)
!            ee2(ii) = xs0(1)
        end do

        p1 = inv_sqrt_pi*inv_y2*exp(-yy*yy)
        xs1(:) = xs1(:) * p1

    end if
!    write(8,1), E0, (ee0(ii), ii = 1, 16), yy, (ee1(ii), ii = 1, 16), &
!                (ee2(ii), ii = 1, 16)
!    1 format(100es15.7)

end subroutine

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

    pt2 = pt2 + 1

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

subroutine URES_PTABLE(iso, erg, xs, IFF, r)
    use RANDOMS, only: rang
    integer, intent(in) :: iso
    real(8), intent(in) :: erg
    real(8), intent(in) :: r
    integer, intent(out):: IFF

    real(8) :: xs(4)
    real(8) :: xs0(4), xs1(4)
    integer :: col0, col1
    type(UNRtype), pointer :: ac
    real(8), pointer :: C0(:), C1(:)
    integer :: N, M, ierg, mm
    real(8) :: f

    integer :: ierg0, ierg1
    real(8) :: ipfac0, ipfac1

    ac => ace(iso) % UNR

    N = ac % N; M = ac % M
    IFF     = ac % IFF
    
    ! 1-1. Select PTABLE with erg
    do ierg = 1, N-1
        if(ac % E(ierg) <= erg .and. erg < ac % E(ierg+1)) exit
    enddo

    !==== NOTES ====
    ! 1: CDF
    ! 2: Total XS
    ! 3: Elastic XS
    ! 4: Fission XS
    ! 5: (n,gamma) XS

    ! 2. Call Random to choose column

    C0 => ac % P(ierg  ,1,:)
    C1 => ac % P(ierg+1,1,:)

    col0 = M; col1 = M ! To avoid crash...
    do mm = 1, M-1
        ! ACE format implies that CDF monotonically increases
        if( r >= C0(mm) .and. r < C0(mm+1)) col0 = mm
        if( r >= C1(mm) .and. r < C1(mm+1)) col1 = mm
    enddo
    
    xs0 = ac % P(ierg  ,2:5,col0)
    xs1 = ac % P(ierg+1,2:5,col1)

    if(IFF == 1) then
        call getierg(iso, ierg0, ac % E(ierg)) 
        call getierg(iso, ierg1, ac % E(ierg+1))
        ipfac0 = max(0D0,min(1D0,(erg-ace(iso)%E(ierg0))/(ace(iso)%E(ierg0+1)-ace(iso)%E(ierg0))))
        ipfac1 = max(0D0,min(1D0,(erg-ace(iso)%E(ierg1))/(ace(iso)%E(ierg1+1)-ace(iso)%E(ierg1))))
        xs0(1) = xs0(1) * (ace(iso) % sigt(ierg0) * (1D0-ipfac0) + ace(iso) % sigt(ierg0+1) * ipfac0)
        xs1(1) = xs1(1) * (ace(iso) % sigt(ierg1) * (1D0-ipfac1) + ace(iso) % sigt(ierg1+1) * ipfac1)

        if(allocated(ace(iso)%sigel)) then
            xs0(2) = xs0(2) * (ace(iso) % sigel(ierg0) * (1D0-ipfac0) + ace(iso) % sigel(ierg0+1) * ipfac0)
            xs1(2) = xs1(2) * (ace(iso) % sigel(ierg1) * (1D0-ipfac1) + ace(iso) % sigel(ierg1+1) * ipfac1)
        else
            xs0(2) = 0; xs1(2) = 0
        endif
        if(allocated(ace(iso)%sigf)) then
            xs0(3) = xs0(3) * (ace(iso) % sigf(ierg0) * (1D0-ipfac0) + ace(iso) % sigf(ierg0+1) * ipfac0)
            xs1(3) = xs1(3) * (ace(iso) % sigf(ierg1) * (1D0-ipfac1) + ace(iso) % sigf(ierg1+1) * ipfac1)
        else
            xs0(3) = 0; xs1(3) = 0
        endif
        if(allocated(ace(iso)%sigd)) then
            xs0(4) = xs0(4) * (ace(iso) % sigd(ierg0) * (1D0-ipfac0) + ace(iso) % sigd(ierg0+1) * ipfac0)
            xs1(4) = xs1(4) * (ace(iso) % sigd(ierg1) * (1D0-ipfac1) + ace(iso) % sigd(ierg1+1) * ipfac1)
        else
            xs0(4) = 0; xs1(4) = 0
        endif
    endif

    ! 3. Interpolate Probability Table
    select case(ac % INT)
    case(2) ! Linear
        f  = (erg-ac % E(ierg))/(ac % E(ierg+1) - ac % E(ierg))
        xs = xs0 *  (1D0-f) + xs1  * f
    case(5) ! Log-Log
        f  = (log(erg)-log(ac % E(ierg))) / (log(ac%E(ierg+1))-log(ac%E(ierg)))
        xs = xs0 ** (1D0-f) * xs1 ** f
    case default
        write(*,*) 'INVALID INTERPOLATION...', ac % INT
    end select
    
    if(associated(C0)) nullify(C0)
    if(associated(C1)) nullify(C1)
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

!    if(iff == 1) then
!        xs(1) = xs(1) * xs0(1)
!        xs(2) = xs(2) * xs0(2)
!        xs(4) = xs(4) * xs0(3)
!        xs(5) = xs(5) * xs0(3)
!        xs(6) = xs(6) * xs0(4)
!        xs(3) = xs(6) + xs(4)
!    else
        xs(1) = xs0(1)
        xs(2) = xs0(2)
        xs(4) = xs0(3)
        xs(5) = xs(4) * getnu(iso,erg)
        xs(6) = xs0(4)
        xs(3) = xs(6) + xs(4)
!    idx = ace(iso) % UNR % unridx
!    call getiueg(erg, ierg)
!    ipfac = max(0D0,min(1D0,(erg-ueggrid(ierg))/(ueggrid(ierg+1)-ueggrid(ierg))))
!    xs = 0D0
!    xs(1) = ace(iso) % UEG % sigt_ures(ierg, urn(idx)) * (1D0-ipfac) + &
!            ace(iso) % UEG % sigt_ures(ierg+1, urn(idx)) * ipfac
!    if(allocated(ace(iso)%UEG%sigel)) &
!    xs(2) = ace(iso) % UEG % sigel_ures(ierg,urn(idx)) * (1D0-ipfac) + &
!            ace(iso) % UEG % sigel_ures(ierg+1,urn(idx)) * ipfac
!    xs(6) = ace(iso) % UEG % sigd_ures(ierg, urn(idx)) * (1D0-ipfac) + &
!           ace(iso) % UEG % sigd_ures(ierg, urn(idx)) * ipfac
!    xs(3) = xs(6)
!    if(allocated(ace(iso)%UEG%sigf)) then
!    xs(4) = ace(iso) % UEG % sigf_ures(ierg, urn(idx)) * (1D0-ipfac) + &
!           ace(iso) % UEG % sigf_ures(ierg, urn(idx)) * ipfac
!    xs(5) = ace(iso) % UEG % signuf_ures(ierg, urn(idx)) * (1D0-ipfac) + &
!           ace(iso) % UEG % signuf_ures(ierg, urn(idx)) * ipfac
!    xs(3) = xs(6) + xs(4)
!    endif

end subroutine
!function getnudel (iso_,erg0) result (nu)
!    integer, intent(in) :: iso_
!    real(8), intent(in) :: erg0
!    real(8) :: nu
!    integer :: ierg, pt1, pt2, pt3
!    real(8) :: ipfac
!    type (AceFormat), pointer :: ac
!
!    nu = 0
!    ac => ace(iso_)
!    if(ac % nu_block_exist == .false.) return
!
!   !if (ac%nu_del%NR /= 0) then
!   !   print *, "WARNING :: nu_del is not lin-lin for", ace(iso_)%library
!   !   print *, ace(iso_)%nu_del%NR
!   !   print *, ac%nu_del%NR
!   !   print *, iso_
!   !endif
!   ! 1. binary search to find ierg of erg in E(:)
!   pt1 = 1
!   pt2 = ac % nu_del % NE
!   Do
!      if(pt2 - pt1 == 1) exit
!      pt3 = (pt2 + pt1)/2
!      if(erg0 >= ac%nu_del%E(pt3)) then
!      pt1 = pt3
!      else
!      pt2 = pt3
!      endif
!   Enddo
!   ierg = pt1 !store low bound energy index
!
!   if (ac % nu_del % NR /= 0) print *, "ENDF interpolation required for Nu "
!   ! 2. calculate interpolation factor
!   ipfac = max(0.d0, min(1.d0,(erg0-ac%nu_del%E(ierg))/(ac%nu_del%E(ierg+1)-ac%nu_del%E(ierg))))
!
!   ! 3. linear-linear interpolation
!   nu = ac%nu_del%F(ierg) + ipfac*(ac%nu_del%F(ierg+1)-ac%nu_del%F(ierg))
!
!end function

end module
