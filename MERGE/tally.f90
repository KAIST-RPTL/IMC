module tally
    use mpi
    use geometry_header 
    use constants
    use particle_header, only : Particle
    use FMFD_HEADER, only: dfm, fm0, fm1, dcm, p_dep_mc, p_dep_dt
    use VARIABLES
    USE XS_header, ONLY: n_group
	
    implicit none
    
    
    type :: LocalCoord
        integer :: cell      = NONE
        integer :: universe  = NONE
        integer :: lattice   = NONE
        integer :: lattice_x = NONE
        integer :: lattice_y = NONE
        integer :: lattice_z = NONE
      contains 
        procedure :: reset => reset_coord
    end type LocalCoord

    type :: CoordStruct
        integer          :: n_coord
        type(LocalCoord) :: coord(MAX_COORD) 
        real(8)          :: vol
        integer          :: flag ! 0 if pin 
                                 ! 1 if cell
    end type
    
    !> Variables
    type(CoordStruct), allocatable, target :: TallyCoord(:)
    real(8), allocatable :: TallyFlux(:), TallyPower(:), tally_buf(:)
    real(8), allocatable :: tally1(:), tally2(:)
    
    ! tally parameters
    logical:: tallyon = .false.
    logical:: meshon  = .false.
    integer:: n_type    ! # of tally types
    integer:: n_tgroup  ! # of tally groups
    integer:: n_tcycle  ! # of tally cycles
    real(8), allocatable:: k_eff(:,:)                 ! (n_batch,n_tot)
    real(8), allocatable:: MC_tally(:,:,:,:,:,:,:)    ! (batches,cycles,types,groups,x,y,z)
    real(8), allocatable:: MC_thread(:,:,:,:,:)       ! (types,groups,x,y,z)
    integer, allocatable:: ttally(:)                  ! (types)
    real(8), allocatable:: tgroup(:)                  ! (groups)
    integer, allocatable:: tcycle(:)                  ! (cycles)
    real(8), allocatable:: MC_stally(:,:,:,:,:,:,:,:) ! (batches,cycles,2,groups,x,y,z,6)
    real(8), allocatable:: MC_sthread(:,:,:,:,:,:)    ! (2,groups,x,y,z,6)
    real(8), allocatable:: MC_scat(:,:,:,:,:,:,:)     ! (batches,cycles,groups,groups,x,y,z)
    real(8), allocatable:: MC_scatth(:,:,:,:,:)       ! (groups,groups,x,y,z)
    logical :: meshon_tet_vrc = .false. 
	real(8) :: mesh_power 
	
    contains
    
    elemental subroutine reset_coord(this)
        class(LocalCoord), intent(inout) :: this
        
        this % cell     = NONE
        this % universe = NONE
        this % lattice  = NONE
        this % lattice_x = NONE
        this % lattice_y = NONE
        this % lattice_z = NONE
        
    end subroutine reset_coord


    function FindTallyBin(p) result (idx)
        type(particle), intent(in) :: p 
        integer :: idx(4)
        integer :: i_bin, i_coord, i, j
        integer, dimension(6) :: A, B, C 
        
        idx = -1
        Bin: do i_bin = 1, size(TallyCoord)
            !if (p%n_coord /= TallyCoord(i_bin)%n_coord) cycle Bin
            do i_coord = 1, TallyCoord(i_bin)%n_coord-1! p%n_coord-1
                A(1) = p%coord(i_coord)%cell
                A(2) = p%coord(i_coord)%universe
                A(3) = p%coord(i_coord)%lattice
                A(4) = p%coord(i_coord)%lattice_x
                A(5) = p%coord(i_coord)%lattice_y
                A(6) = p%coord(i_coord)%lattice_z
                
                B(1) = TallyCoord(i_bin)%coord(i_coord)%cell
                B(2) = TallyCoord(i_bin)%coord(i_coord)%universe
                B(3) = TallyCoord(i_bin)%coord(i_coord)%lattice
                B(4) = TallyCoord(i_bin)%coord(i_coord)%lattice_x
                B(5) = TallyCoord(i_bin)%coord(i_coord)%lattice_y
                B(6) = TallyCoord(i_bin)%coord(i_coord)%lattice_z
                C(:) = abs(A(:)-B(:))
                if (sum(C) /= 0) cycle Bin
                
            enddo 
            i_coord = TallyCoord(i_bin)%n_coord !p%n_coord
            A(1) = p%coord(i_coord)%cell     * TallyCoord(i_bin)%flag
            A(2) = p%coord(i_coord)%universe * TallyCoord(i_bin)%flag
			if (do_gmsh) then 
				A(1) = 0
				A(2) = p%coord(i_coord)%universe
			endif 
            A(3) = p%coord(i_coord)%lattice
            A(4) = p%coord(i_coord)%lattice_x
            A(5) = p%coord(i_coord)%lattice_y
            A(6) = p%coord(i_coord)%lattice_z
            
            B(1) = TallyCoord(i_bin)%coord(i_coord)%cell
            B(2) = TallyCoord(i_bin)%coord(i_coord)%universe
            B(3) = TallyCoord(i_bin)%coord(i_coord)%lattice
            B(4) = TallyCoord(i_bin)%coord(i_coord)%lattice_x
            B(5) = TallyCoord(i_bin)%coord(i_coord)%lattice_y
            B(6) = TallyCoord(i_bin)%coord(i_coord)%lattice_z
			
            C(:) = abs(A(:)-B(:))
            if (sum(C) /= 0) cycle Bin
            idx(1) = i_bin
            idx(2:4) = B(4:6)
            exit Bin 

        enddo Bin

    end function 


! =============================================================================
! SET_MC_TALLY
! =============================================================================
subroutine SET_MC_TALLY
    use VARIABLES, only: icore, score, n_batch, n_act, t_totcyc, n_totcyc
    use FMFD_HEADER, only: nfm
    implicit none

    allocate(k_eff(n_batch,n_totcyc))
    if ( tallyon ) then
		! --- MESH_GRID information not specified / Alter tallyon -> .FALSE.
		IF(nfm(1) == 0 .OR. nfm(2) == 0 .OR. nfm(3) == 0) THEN
			IF(icore == score) WRITE(*,'(A)') '(WARNING) TALLY_MESH is .TRUE. w/o provision of MESH_GRID information'
			tallyon = .FALSE.
		! --- MESH GRID INFORMATION given
		ELSE
			allocate(MC_tally (n_batch,n_act,n_type,n_tgroup,nfm(1),nfm(2),nfm(3)))
			allocate(MC_thread(              n_type,n_tgroup,nfm(1),nfm(2),nfm(3)))
			MC_tally  = 0
			MC_thread = 0
		END IF
    end if
end subroutine SET_MC_TALLY

! =============================================================================
! SET_MC_TALLY_ADJ
! =============================================================================
SUBROUTINE SET_MC_TALLY_ADJ
	use VARIABLES, only: icore, score, do_IFP, do_IFP_LONG, nfm_adj, E_mode
	IMPLICIT NONE
	IF((tally_adj_flux .EQ. .TRUE.) .AND. (tally_adj_spect .EQ. .TRUE.)) THEN
		IF(icore == score) WRITE(*,'(A)') '(WARNING) tally_adj_flux and tally_adj_spect cannot be supported at the same time (at the moment...)' 
		STOP
	END IF
	!> RELATED TO TALLYING ADJ_FLUX
	IF(tally_adj_flux) THEN
		IF(nfm_adj(1) == 0 .OR. nfm_adj(2) == 0 .OR. nfm_adj(3) == 0) THEN
			WRITE(*,'(A)') '(WARNING) tally_adj_flux is .TRUE. w/o provision of MESH_ADJ_GRID information'
			STOP
		END IF
		IF(E_mode == 1) THEN
			IF(icore == score) WRITE(*,'(A)') '(WARNING) AT THE MOMENT, ADJOINT FLUX DISTRIBUTION CALCULATION FOR C.E. Monte-Carlo NOT SUPPORTED'
			STOP
		ELSE
			num_adj_group = n_group
		END IF
		if(icore==score) print '(A30)', '    ADJOINT DISTRIBUTION TALLY PERFORMED...' 
		do_IFP_LONG = .TRUE.
	!> RELATED TO TALLYING ADJ SPECTRUM
	ELSE IF(tally_adj_spect) THEN
		if(icore==score) print '(A30)', '    ADJOINT SPECTRUM     TALLY PERFORMED...' 
		do_IFP_LONG = .TRUE.
	ELSE
		do_IFP_LONG = .FALSE.
	END IF
	IF(do_IFP_LONG) THEN
		IF(NOT(do_IFP)) THEN
			IF(icore == score) WRITE(*,'(A)') '(WARNING) TALLYING ADJOINT SPECTRUM DEMANDS IFP OPTION TO BE TRUE'
			STOP
		END IF
	END IF	
END SUBROUTINE SET_MC_TALLY_ADJ

! =============================================================================
! MESH_DISTANCE
! =============================================================================
subroutine MESH_DISTANCE (p,i_xyz,d_mesh,inside_mesh,income_mesh,i_surf)
    type(particle), intent(in) :: p
    integer, intent(inout) :: i_xyz(3)
    real(8), intent(inout) :: d_mesh
    logical, intent(inout) :: inside_mesh 
    integer, intent(inout) :: income_mesh
    integer, intent(inout) :: i_surf
    real(8) :: xyz(3), uvw(3), xyz1(3)
    real(8) :: d_temp(6)
    integer :: ij
    
    ! Find lattice index in FM grid
    xyz(:) = p%coord(1)%xyz(:)
    uvw(:) = p%coord(1)%uvw(:)
    i_xyz  = FM_ID(p%coord(1)%xyz(:))
    
    ! the particle is inside the FM grid
    inside_mesh = INSIDE(xyz)
    income_mesh = 0

    if ( inside_mesh ) then
        d_temp(1) = ((dfm(1)*(i_xyz(1)-1)+fm0(1))-xyz(1))/uvw(1)   ! x0
        d_temp(2) = ((dfm(1)*(i_xyz(1)  )+fm0(1))-xyz(1))/uvw(1)   ! x1
        d_temp(3) = ((dfm(2)*(i_xyz(2)-1)+fm0(2))-xyz(2))/uvw(2)   ! y0
        d_temp(4) = ((dfm(2)*(i_xyz(2)  )+fm0(2))-xyz(2))/uvw(2)   ! y1
        d_temp(5) = ((dfm(3)*(i_xyz(3)-1)+fm0(3))-xyz(3))/uvw(3)   ! z0
        d_temp(6) = ((dfm(3)*(i_xyz(3)  )+fm0(3))-xyz(3))/uvw(3)   ! z1
        
        do ij = 1, 6
        if ( d_temp(ij) > 0 .and. d_mesh > d_temp(ij) ) then
            d_mesh = d_temp(ij)
            i_surf = ij
        end if
        end do

    ! the particle is outside the FM grid
    else
        d_temp(1) = (fm0(1)-xyz(1))/uvw(1)   ! x0
        d_temp(2) = (fm1(1)-xyz(1))/uvw(1)   ! x1
        d_temp(3) = (fm0(2)-xyz(2))/uvw(2)   ! x0
        d_temp(4) = (fm1(2)-xyz(2))/uvw(2)   ! y1
        d_temp(5) = (fm0(3)-xyz(3))/uvw(3)   ! z0
        d_temp(6) = (fm1(3)-xyz(3))/uvw(3)   ! z1

        do ij = 1, 6
        if ( d_temp(ij) > 0 .and. d_mesh > d_temp(ij) ) then
            xyz1(:) = xyz(:) + d_temp(ij)*uvw(:)
            select case(ij)
            case(1,2)
                if ( xyz1(2) < fm0(2) .or. fm1(2) < xyz1(2) ) cycle
                if ( xyz1(3) < fm0(3) .or. fm1(3) < xyz1(3) ) cycle
            case(3,4)
                if ( xyz1(3) < fm0(3) .or. fm1(3) < xyz1(3) ) cycle
                if ( xyz1(1) < fm0(1) .or. fm1(1) < xyz1(1) ) cycle
            case(5,6)
                if ( xyz1(1) < fm0(1) .or. fm1(1) < xyz1(1) ) cycle
                if ( xyz1(2) < fm0(2) .or. fm1(2) < xyz1(2) ) cycle
            end select
            i_xyz = FM_ID(xyz1)
            d_mesh = d_temp(ij)
            income_mesh = ij
        end if
        end do
    end if

end subroutine

! =============================================================================
! MESH_DISTANCE ADJOINT
! =============================================================================
subroutine MESH_DISTANCE_ADJ (p,i_xyz,d_mesh,inside_mesh,income_mesh,i_surf)
    type(particle), intent(in) :: p
    integer, intent(inout) :: i_xyz(3)
    real(8), intent(inout) :: d_mesh
    logical, intent(inout) :: inside_mesh 
    integer, intent(inout) :: income_mesh
    integer, intent(inout) :: i_surf
    real(8) :: xyz(3), uvw(3), xyz1(3)
    real(8) :: d_temp(6)
    integer :: ij
    
    ! Find lattice index in FM grid
    xyz(:) = p%coord(1)%xyz(:)
    uvw(:) = p%coord(1)%uvw(:)
    i_xyz  = AM_ID(p%coord(1)%xyz(:))
    
    ! the particle is inside the FM grid
    inside_mesh = INSIDE_ADJ(xyz)
    income_mesh = 0

    if ( inside_mesh ) then
        d_temp(1) = ((dfm_adj(1)*(i_xyz(1)-1)+fm0_adj(1))-xyz(1))/uvw(1)   ! x0
        d_temp(2) = ((dfm_adj(1)*(i_xyz(1)  )+fm0_adj(1))-xyz(1))/uvw(1)   ! x1
        d_temp(3) = ((dfm_adj(2)*(i_xyz(2)-1)+fm0_adj(2))-xyz(2))/uvw(2)   ! y0
        d_temp(4) = ((dfm_adj(2)*(i_xyz(2)  )+fm0_adj(2))-xyz(2))/uvw(2)   ! y1
        d_temp(5) = ((dfm_adj(3)*(i_xyz(3)-1)+fm0_adj(3))-xyz(3))/uvw(3)   ! z0
        d_temp(6) = ((dfm_adj(3)*(i_xyz(3)  )+fm0_adj(3))-xyz(3))/uvw(3)   ! z1
        
        do ij = 1, 6
        if ( d_temp(ij) > 0 .and. d_mesh > d_temp(ij) ) then
            d_mesh = d_temp(ij)
            i_surf = ij
        end if
        end do

    ! the particle is outside the FM grid
    else
        d_temp(1) = (fm0_adj(1)-xyz(1))/uvw(1)   ! x0
        d_temp(2) = (fm1_adj(1)-xyz(1))/uvw(1)   ! x1
        d_temp(3) = (fm0_adj(2)-xyz(2))/uvw(2)   ! x0
        d_temp(4) = (fm1_adj(2)-xyz(2))/uvw(2)   ! y1
        d_temp(5) = (fm0_adj(3)-xyz(3))/uvw(3)   ! z0
        d_temp(6) = (fm1_adj(3)-xyz(3))/uvw(3)   ! z1

        do ij = 1, 6
        if ( d_temp(ij) > 0 .and. d_mesh > d_temp(ij) ) then
            xyz1(:) = xyz(:) + d_temp(ij)*uvw(:)
            select case(ij)
            case(1,2)
                if ( xyz1(2) < fm0_adj(2) .or. fm1_adj(2) < xyz1(2) ) cycle
                if ( xyz1(3) < fm0_adj(3) .or. fm1_adj(3) < xyz1(3) ) cycle
            case(3,4)
                if ( xyz1(3) < fm0_adj(3) .or. fm1_adj(3) < xyz1(3) ) cycle
                if ( xyz1(1) < fm0_adj(1) .or. fm1_adj(1) < xyz1(1) ) cycle
            case(5,6)
                if ( xyz1(1) < fm0_adj(1) .or. fm1_adj(1) < xyz1(1) ) cycle
                if ( xyz1(2) < fm0_adj(2) .or. fm1_adj(2) < xyz1(2) ) cycle
            end select
            i_xyz = AM_ID(xyz1)
            d_mesh = d_temp(ij)
            income_mesh = ij
        end if
        end do
    end if
end subroutine MESH_DISTANCE_ADJ

! =============================================================================
! FM_ID finds the x, y, z indice in the fine mesh grid
! =============================================================================
function FM_ID(fmxyz) result(fmid)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: fmid(3)               ! indice
       
    fmid(:) = floor((fmxyz(:)-fm0(:))/dfm(:))+1

end function

! =============================================================================
! AM_ID finds the x, y, z indice in the ADJOINT MESH grid
! =============================================================================
function AM_ID(fmxyz) result(fmid)
	IMPLICIT NONE
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: fmid(3)               ! indice
    fmid(:) = floor((fmxyz(:)-fm0_adj(:))/dfm_adj(:))+1
end function AM_ID

! =============================================================================
! CM_ID finds the x, y, z indice in the fine mesh grid
! =============================================================================
function CM_ID(fmxyz) result(fmid)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: fmid(3)               ! indice
       
    fmid(:) = floor((fmxyz(:)-fm0(:))/dcm(:))+1

end function

! =============================================================================
! FMFD_ID finds the x, y, z indice in the FMFD mesh grid
! =============================================================================
function INSIDE(fmxyz) result(inside_mesh)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: inside_mesh
    integer:: ij
       
    inside_mesh = .true.
    do ij = 1, 3
        if ( fmxyz(ij) < fm0(ij) .or. fmxyz(ij) >= fm1(ij) ) then
            inside_mesh = .false.
            exit
        end if
    end do

end function

! =============================================================================
! INSED_ADJ finds the x, y, z indice in the FMFD mesh grid
! =============================================================================
function INSIDE_ADJ(fmxyz) result(inside_mesh_ADJ)
	IMPLICIT NONE
    real(8), intent(in):: fmxyz(:)  ! coordinate
    LOGICAL:: inside_mesh_ADJ
    integer:: ij
    inside_mesh_ADJ = .true.
    do ij = 1, 3
        if ( fmxyz(ij) < fm0_adj(ij) .or. fmxyz(ij) >= fm1_adj(ij) ) then
            inside_mesh_ADJ = .false.
            exit
        end if
    end do
end function INSIDE_ADJ

! =============================================================================
! OUT_OF_ZZ_ADJ (Determines whether the node index is out of ADJ ZZ-boundary)
! =============================================================================
function OUT_OF_ZZ_ADJ(io,jo)
    logical:: OUT_OF_ZZ_ADJ
    integer, intent(in):: io, jo
    integer:: mo, no
    
    if ( .not. zigzag_adjflux ) then
        OUT_OF_ZZ_ADJ = .false.
        return
    end if
	IF((io <= 0) .OR. (io > nfm_adj(1)) .OR. (jo <= 0) .OR. (jo > nfm_adj(2))) THEN
		OUT_OF_ZZ_ADJ = .TRUE.
		RETURN
	END IF
    OUT_OF_ZZ_ADJ = .false.
    do mo = 1, zz_div_adj
    if ( zzf0_adj(mo) < io .and. io <= zzf0_adj(mo+1) ) then
        no = mo
        exit
    end if
    end do

    if ( zzf1_adj(no) < jo .and. jo <= zzf2_adj(no) ) then
        OUT_OF_ZZ_ADJ = .false.
    else
        OUT_OF_ZZ_ADJ = .true.
    end if

end function OUT_OF_ZZ_ADJ

! =============================================================================
! INITIALIZE MC_THREAD 
! =============================================================================
subroutine TALLY_THREAD_INITIAL(cyc)
    implicit none
    integer, intent(in):: cyc
    integer:: acyc

    if ( cyc > n_inact ) then
    acyc = cyc - n_inact
    MC_thread = 0D0
    ! (MUTED) MC_tally  = 0D0 --- 안그러면 사이클마다 MC_tally 값이 0이 되버림 
    end if

end subroutine

! =============================================================================
! MC_TRK calculates MC parameters such as flux, group contstans by 
! track-length estiamtor
! =============================================================================
subroutine MC_TRK(E0,wgt,distance,macro_xs,id)
    implicit none
    type(Particle):: p
    real(8), intent(in) :: E0
    real(8), intent(in) :: wgt
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(1:3)
    real(8) :: flux, xs
    integer :: ii

    if ( curr_cyc <= n_inact ) return
    flux = wgt * distance
    do ii = 1, n_type
		select case(ttally(ii))
			case(1,11); xs = macro_xs(1)                ! total
			case(2,12); xs = macro_xs(2)                ! absorption
			case(3,13); xs = macro_xs(3)                ! fission
			case(4,14); xs = macro_xs(4)                ! nu fission
			case(5,15); xs = macro_xs(5)                ! k fission
			case(6,16); xs = 0                          ! scattering
			case default; xs = 1D0                      ! flux (type = 0)
		end select
		MC_thread(ii,E2G(E0),id(1),id(2),id(3)) = & 
		MC_thread(ii,E2G(E0),id(1),id(2),id(3)) + flux * xs
    end do
    
end subroutine

! =============================================================================
! MC_TRK calculates MC parameters such as flux, group contstans by 
! track-length estiamtor
! =============================================================================
subroutine MC_TRK_S(cyc,E0,E1,wgt,distance,macro_xs,id)
    implicit none
    type(Particle):: p
    integer, intent(in) :: cyc
    real(8), intent(in) :: E0
    real(8), intent(in) :: E1
    real(8), intent(in) :: wgt
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(1:3)
    real(8) :: flux, xs
    integer :: ii
    
    if ( .not. tallyon ) return
    if ( cyc <= n_inact ) return
    flux = wgt * distance
    xs = macro_xs(1) - macro_xs(2)
    MC_scatth(E2G(E0),E2G(E1),id(1),id(2),id(3)) = & 
    MC_scatth(E2G(E0),E2G(E1),id(1),id(2),id(3)) + flux * xs
    
end subroutine

! =============================================================================
! MC_COL calculates MC parameters such as flux, group contstans by 
! collision estimator
! =============================================================================
subroutine MC_COL(E0,wgt,macro_xs,id)
    real(8), intent(in) :: E0
    real(8), intent(in) :: wgt
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(3)
    real(8) :: flux, xs
    integer:: ii
    
    flux = wgt / macro_xs(1)

    do ii = 1, n_type
    select case(ttally(ii))
    case(1,11); xs = macro_xs(1)                ! total
    case(2,12); xs = macro_xs(2)                ! absorption
    case(3,13); xs = macro_xs(3)                ! fission
    case(4,14); xs = macro_xs(4)                ! nu fission
    case(5,15); xs = macro_xs(5)                ! k fission
    case(6,16); xs = macro_xs(1) - macro_xs(2)  ! scattering
    case default; xs = 1D0
    end select
    MC_thread(ii,E2G(E0),id(1),id(2),id(3)) = & 
    MC_thread(ii,E2G(E0),id(1),id(2),id(3)) + flux * xs
    end do
    
end subroutine

! =============================================================================
! TALLY_SURF calculates FMFD surface parameters like net and particle current
! =============================================================================
subroutine TALLY_SURF (inside,income, is, id, E0, uvw, wgt, bc)
    use FMFD_HEADER, only: nfm
    implicit none
    logical, intent(in) :: inside
    integer, intent(in) :: income
    integer, intent(in) :: is, id(3)
    real(8), intent(in) :: E0
    real(8), intent(in) :: uvw(3)
    real(8), intent(in) :: wgt
    integer, intent(in) :: bc
    integer:: dir

    if ( .not. tallyon ) return
    
    ! inner nodes
    if ( inside ) then 
        ! surface partial current
        select case(is)
        case(1,3,5)
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) - wgt
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) + wgt / abs(uvw(is/2+1))
        case(2,4,6)
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) + wgt
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) + wgt / abs(uvw(is/2))
        end select

        ! boundary condition
        if ( bc == 2 ) then
        select case(is)
        case(1,3,5)
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) + wgt
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) + wgt / abs(uvw(is/2+1))
        case(2,4,6)
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(1,E2G(E0),id(1),id(2),id(3),is) - wgt 
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) = &
            MC_sthread(2,E2G(E0),id(1),id(2),id(3),is) + wgt / abs(uvw(is/2))
        end select
        end if
        return
    end if

    ! boundary nodes
    select case(income)
    case(1)
        MC_sthread(1,E2G(E0),1,id(2),id(3),1) = &
        MC_sthread(1,E2G(E0),1,id(2),id(3),1) + wgt
        MC_sthread(2,E2G(E0),1,id(2),id(3),1) = &
        MC_sthread(2,E2G(E0),1,id(2),id(3),1) + wgt / abs(uvw(1))
    case(2)
        MC_sthread(1,E2G(E0),nfm(1),id(2),id(3),2) = &
        MC_sthread(1,E2G(E0),nfm(1),id(2),id(3),2) - wgt
        MC_sthread(2,E2G(E0),nfm(1),id(2),id(3),2) = &
        MC_sthread(2,E2G(E0),nfm(1),id(2),id(3),2) + wgt / abs(uvw(1))
    case(3)
        MC_sthread(1,E2G(E0),id(1),1,id(3),3) = &
        MC_sthread(1,E2G(E0),id(1),1,id(3),3) + wgt
        MC_sthread(2,E2G(E0),id(1),1,id(3),3) = &
        MC_sthread(2,E2G(E0),id(1),1,id(3),3) + wgt / abs(uvw(2))
    case(4)
        MC_sthread(1,E2G(E0),id(1),nfm(2),id(3),4) = &
        MC_sthread(1,E2G(E0),id(1),nfm(2),id(3),4) - wgt
        MC_sthread(2,E2G(E0),id(1),nfm(2),id(3),4) = &
        MC_sthread(2,E2G(E0),id(1),nfm(2),id(3),4) + wgt / abs(uvw(2))
    case(5)
        MC_sthread(1,E2G(E0),id(1),id(2),1,5) = &
        MC_sthread(1,E2G(E0),id(1),id(2),1,5) + wgt
        MC_sthread(2,E2G(E0),id(1),id(2),1,5) = &
        MC_sthread(2,E2G(E0),id(1),id(2),1,5) + wgt / abs(uvw(3))
    case(6)
        MC_sthread(1,E2G(E0),id(1),id(2),nfm(3),6) = &
        MC_sthread(1,E2G(E0),id(1),id(2),nfm(3),6) - wgt
        MC_sthread(2,E2G(E0),id(1),id(2),nfm(3),6) = &
        MC_sthread(2,E2G(E0),id(1),id(2),nfm(3),6) + wgt / abs(uvw(3))
    end select
            
end subroutine

! =============================================================================
! 
! =============================================================================
function E2G(ee)
    integer:: E2G
    real(8), intent(in):: ee
    integer:: ii

    E2G = 1
    do ii = 1, n_tgroup-1, 1
    if ( ee < tgroup(ii) ) then
        E2G = n_tgroup-ii+1
        return
    end if
    end do

end function E2G

! =============================================================================
! FUNCTION FOR MAPPING ENERGY TO GROUP (FOR SPECTRUM TALLY)
! =============================================================================
FUNCTION E2G_SPECTRUM(ee)
    integer:: E2G_SPECTRUM
    real(8), intent(in):: ee
    integer:: ii
	! --- IF THE GIVEN ENERGY EXCEEDS THE LARGEST ENERGY RANGE VALUE
	IF(ee > egroup_spect(1)) THEN
		E2G_SPECTRUM = 0
		RETURN
	END IF
	! --- IF THE GIVEN ENERGY IS BELOW THE SMALLEST ENERGY RANGE VALUE
	IF(ee < egroup_spect(n_egroup_spect)) THEN
		E2G_SPECTRUM = 0
		RETURN		
	END IF
	! --- LOWEST ENERGY RANGE BIN INDEX (제공된 에너지 군의 최소값 보다 작을 시)
    E2G_SPECTRUM = n_egroup_spect - 1
	! --- egroup_spect descends in its value
    do ii = 1, n_egroup_spect-1, 1
    if ( ee > egroup_spect(ii) ) then
        E2G_SPECTRUM = ii-1
        return
    end if
    end do
END FUNCTION E2G_SPECTRUM

! =============================================================================
! NORM_TALLY normalizes cycle-wise FMFD parameters
! =============================================================================
subroutine NORM_TALLY(bat,cyc)
    use VARIABLES, only: n_inact
    implicit none
    integer, intent(in):: bat, cyc
    integer:: acyc

    if ( cyc <= n_inact ) return
    if ( bat == 0 ) return
    acyc = cyc-n_inact

    !> gather thread tally parameters
    MC_tally(bat,acyc,:,:,:,:,:) = &
    MC_tally(bat,acyc,:,:,:,:,:) + MC_thread(:,:,:,:,:)

end subroutine

! =============================================================================
! PROCESS_FMFD deals with MPI process and average quantities
! =============================================================================
subroutine PROCESS_TALLY(bat,cyc)
    use VARIABLES, only: n_inact, icore, score, ngen, ierr
    use FMFD_HEADER, only: nfm, v_fm, a_fm
    use MPI, only: MPI_SUM
    implicit none
    !> MPI derived type reduce parameters 
    integer, intent(in):: bat, cyc
    real(8), allocatable:: MC_temp0(:,:,:,:,:), MC_temp1(:,:,:,:,:)
    real(8), allocatable:: MC_stemp0(:,:,:,:,:,:), MC_stemp1(:,:,:,:,:,:)
    real(8), allocatable:: MC_temps0(:,:,:,:,:), MC_temps1(:,:,:,:,:)
    integer:: acyc, dsize
    integer:: ii

    if ( cyc <= n_inact ) return
    if ( bat == 0 ) return
    acyc = cyc-n_inact

    ! data gathering
    allocate(MC_temp0(n_type,n_tgroup,nfm(1),nfm(2),nfm(3)))
    allocate(MC_temp1(n_type,n_tgroup,nfm(1),nfm(2),nfm(3)))
    MC_temp0  = MC_tally(bat,acyc,:,:,:,:,:)
    MC_temp1  = 0D0

    dsize = nfm(1)*nfm(2)*nfm(3)*n_type*n_tgroup
    call MPI_REDUCE(MC_temp0(:,:,:,:,:),MC_temp1(:,:,:,:,:), &
                    dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)

    if ( icore == score ) then
		MC_tally(bat,acyc,:,:,:,:,:) = MC_temp1(:,:,:,:,:)/(dble(ngen)*v_fm)
    end if
    deallocate(MC_temp0,MC_temp1)

    if ( icore == score ) then
		do ii = 1, n_type
		if ( ttally(ii) > 10 ) &
		MC_tally(bat,acyc,ii,:,:,:,:) = &
		MC_tally(bat,acyc,ii,:,:,:,:) / MC_tally(bat,acyc,n_type,:,:,:,:)
		end do
    end if

    if ( do_burn .and. icore == score .and. ttally(1) == 5 ) &
        p_dep_mc(acyc,:,:,:) = MC_tally(bat,acyc,1,1,:,:,:)
end subroutine 

end module 
