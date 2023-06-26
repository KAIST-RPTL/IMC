module HGC

	use FMFD_header 
    ! ==== surface numbering ====
    !  1    2    3    4    5    6
    ! x0 / x1 / y0 / y1 / z0 / z1
    ! ===========================
	implicit none

    ! =========================================================================
    ! HGC type 
    type :: HGC_parameters
        real(8) :: phi
        real(8) :: sig_t
        real(8) :: sig_a 
		real(8) :: chi
        real(8) :: nusig_f 
        real(8) :: kasig_f 
        real(8) :: Jn(6) 
        real(8) :: J0(6)
        real(8) :: J1(6)
    end type
    real(8), allocatable :: HGC_scat(:,:,:,:,:)
	
	logical :: do_HGC = .false.
	integer :: g_HGC
    real(8),allocatable :: E_grid(:)
	
    type(HGC_parameters), allocatable :: HGC_dat(:,:,:,:)

	

contains 


! =============================================================================
! SELECT_GROUP determines the HGC group of the particle 
! =============================================================================

function select_group(E) result(g) 
	real(8), intent(in) :: E 
	real(8) :: g 
	integer :: i 
	g = g_HGC
	do i = 1, g_HGC-1
		if (E < E_grid(i)) then 
			g = i
			return 
		endif 
	enddo 
end function 



! =============================================================================
! FMFD_TRK calculates FMFD parameters such as flux, group contstans by 
! track-length estiamtor
! =============================================================================
subroutine HGC_TRK(wgt,E,distance,macro_xs,id)
    implicit none
    real(8), intent(in) :: wgt
    real(8), intent(in) :: E
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(1:3)
    real(8) :: flux
	integer :: g
    
	g = select_group(E)
    flux = wgt * distance
    
	!$omp atomic
    HGC_dat(id(1),id(2),id(3),g) % phi = & 
    HGC_dat(id(1),id(2),id(3),g) % phi + flux
	!$omp atomic
    HGC_dat(id(1),id(2),id(3),g) % sig_t = &
    HGC_dat(id(1),id(2),id(3),g) % sig_t + flux*macro_xs(1)
	!$omp atomic
    HGC_dat(id(1),id(2),id(3),g) % sig_a = &
    HGC_dat(id(1),id(2),id(3),g) % sig_a + flux*macro_xs(2)
	!!$omp atomic
    !HGC_dat(id(1),id(2),id(3),g) % chi = &
    !HGC_dat(id(1),id(2),id(3),g) % chi 
	
	!$omp atomic
    HGC_dat(id(1),id(2),id(3),g) % nusig_f = &
    HGC_dat(id(1),id(2),id(3),g) % nusig_f + flux*macro_xs(4)
	!$omp atomic
    HGC_dat(id(1),id(2),id(3),g) % kasig_f = &
    HGC_dat(id(1),id(2),id(3),g) % kasig_f + flux*macro_xs(5)
	
end subroutine


! =============================================================================
! HGC_COL calculates HGC parameters such as flux, group contstans by 
! collision estiamtor
! =============================================================================
subroutine HGC_COL(wgt,Ein,Eout,distance,macro_xs,id)
    implicit none
    real(8), intent(in) :: wgt
    real(8), intent(in) :: Ein,Eout
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(1:3)
    real(8) :: flux
	integer :: gin,gout
    
	gin  = select_group(Ein)
	gout = select_group(Eout)
    flux = wgt / macro_xs(1)
	
	!$omp atomic
	HGC_scat(id(1),id(2),id(3),gin,gout) = &
	HGC_scat(id(1),id(2),id(3),gin,gout) + flux*(macro_xs(1)-macro_xs(2))

end subroutine


! =============================================================================
! HGC_SURF calculates HGC surface parameters like net and particle current
! =============================================================================
subroutine HGC_SURF (inside,income, is, id, uvw, wgt,E, bc)
    logical, intent(in) :: inside
    integer, intent(in) :: income
    integer, intent(in) :: is, id(3)
    real(8), intent(in) :: uvw(3)
    real(8), intent(in) :: wgt
    real(8), intent(in) :: E
    integer, intent(in) :: bc
    integer :: dir
	integer :: g
    
	g = select_group(E)
	
    ! inner nodes
    if ( inside ) then 
        ! surface partial current
        select case(is)
        case(1,3,5)
			!$omp atomic
            HGC_dat(id(1),id(2),id(3),g)%J0(is) = &
            HGC_dat(id(1),id(2),id(3),g)%J0(is) + wgt
        case(2,4,6)
			!$omp atomic
            HGC_dat(id(1),id(2),id(3),g)%J1(is) = &
            HGC_dat(id(1),id(2),id(3),g)%J1(is) + wgt
        end select

        ! boundary condition
        if ( bc == 2 ) then
        select case(is)
        case(1,3,5)
			!$omp atomic
            HGC_dat(id(1),id(2),id(3),g)%J1(is) = &
            HGC_dat(id(1),id(2),id(3),g)%J1(is) + wgt
        case(2,4,6)
			!$omp atomic
            HGC_dat(id(1),id(2),id(3),g)%J0(is) = &
            HGC_dat(id(1),id(2),id(3),g)%J0(is) + wgt
        end select
        end if
        return
    end if

    ! boundary nodes
    select case(income)
    case(1)
		!$omp atomic
        HGC_dat(1,id(2),id(3),g)%J1(1) = &
        HGC_dat(1,id(2),id(3),g)%J1(1) + wgt
    case(2)
		!$omp atomic
        HGC_dat(nfm(1),id(2),id(3),g)%J0(2) = &
        HGC_dat(nfm(1),id(2),id(3),g)%J0(2) + wgt
    case(3)
		!$omp atomic
        HGC_dat(id(1),1,id(3),g)%J1(3) = &
        HGC_dat(id(1),1,id(3),g)%J1(3) + wgt
    case(4)
		!$omp atomic
        HGC_dat(id(1),nfm(2),id(3),g)%J0(4) = &
        HGC_dat(id(1),nfm(2),id(3),g)%J0(4) + wgt
    case(5)
		!$omp atomic
        HGC_dat(id(1),id(2),1,g)%J1(5) = &
        HGC_dat(id(1),id(2),1,g)%J1(5) + wgt
    case(6)
		!$omp atomic
        HGC_dat(id(1),id(2),nfm(3),g)%J0(6) = &
        HGC_dat(id(1),id(2),nfm(3),g)%J0(6) + wgt
    end select
            
end subroutine


end module 
