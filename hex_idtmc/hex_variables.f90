module hex_variables
    use CONSTANTS,          only: eVtoJoule
    use omp_lib 
    use mpi 
    use geometry_header,    only : lattices, find_lat_idx, lattice_coord
    use variables
    use particle_header,    only : particle
    use FMFD_HEADER
    use TH_HEADER, only: th_on, pp
	
    use FMFD_header
	
	implicit none
	
    real(8) :: dduct ! distance between the centreline of an edge pin and the corresponding assembly boundary
	real(8), allocatable :: s_hf(:,:,:,:) ! area of each fine-mesh surface
	real(8), allocatable :: v_hf(:,:,:) ! volume of each fine-mesh cell
	integer :: x_max, y_max, z_max, n_cells_hf ! maximum global fine-mesh x, y, z coordinates and resulting number of cells
	
	! hex zigzag
	integer, allocatable :: hc_zmin(:)
	integer, allocatable :: hc_zmax(:)
    
    !! CMFD VARIABLES !!
    real(8) :: hc_v, hc_sr, hc_sa ! coarse mesh node sizes
    
    real(8), allocatable :: hc_t(:,:,:), hcD(:,:,:), hc_a(:,:,:), hc_nf(:,:,:) ! cross-sections
    real(8), allocatable :: hcphi0(:,:,:), hcphi1(:,:,:), hcJ0(:,:,:,:), hcJ1(:,:,:,:), hcJn(:,:,:,:) ! fluxes and correction factors
    real(8), allocatable :: hcDt(:,:,:,:), hcDh(:,:,:,:), Mhc(:,:,:,:) ! FMFD elements
    
    real(8), allocatable :: s_hc(:,:,:,:), v_hc(:,:,:) ! fine-mesh cell sizes
    
	real(8), allocatable :: A_hc(:,:) ! linear equations solver variables
	real(8), allocatable :: b_hc(:), x_hc(:), b_hc_old(:)
	integer :: diags_hc(9)
	
    contains
	
	! function hc_in_zz(i, j) to check if coarse mesh cell (i, j) falls within the reactor area
    logical function hc_in_zz(i, j)
        integer, intent(in) :: i, j
	    if (i < 1 .or. i > size(hc_zmin)) then
		    hc_in_zz = .false.
			return
		end if
        if (j < hc_zmin(i) .or. j > hc_zmax(i)) then
		    hc_in_zz = .false.
		else
		    hc_in_zz = .true.
		end if
    end function
	
	logical function hc_in_zz0(i, j, k)
	    integer, intent(in) :: i, j, k
		
	    if (.not. hc_in_zz(i, j)) then
		    hc_in_zz0 = .false.
			return
		else if (.not. hc_in_zz(i-1, j)) then
		    hc_in_zz0 = .false.
			return
		else if (.not. hc_in_zz(i+1, j)) then
		    hc_in_zz0 = .false.
			return
		else if (.not. hc_in_zz(i-1, j+1)) then
		    hc_in_zz0 = .false.
			return
		else if (.not. hc_in_zz(i+1, j-1)) then
		    hc_in_zz0 = .false.
			return
		else if (.not. hc_in_zz(i, j-1)) then
		    hc_in_zz0 = .false.
			return
		else if (.not. hc_in_zz(i, j+1)) then
		    hc_in_zz0 = .false.
			return
		else if (k == 1 .or. k == z_max) then
		    hc_in_zz0 = .false.
			return
		end if
	    hc_in_zz0 = .true.
	    return
	end function
	
	subroutine hf_initialise ! equivalent to FMFD_initialize()
        implicit none
        integer:: ij
	    
        ! initialization
        fm(:,:,:) % phi     = 0 
        fm(:,:,:) % sig_t   = 0 
        fm(:,:,:) % sig_a   = 0 
        fm(:,:,:) % nusig_f = 0 
        fm(:,:,:) % kappa   = 0 
        do ij = 1, 8
            fm(:,:,:) % J0(ij)   = 0 
            fm(:,:,:) % J1(ij)   = 0 
        end do 
	end subroutine hf_initialise
	
	subroutine hf_initialise_thread ! equivalent to FMFD_initialize_thread()
	    integer :: i
		
	    if (.not. allocated(fm_thread)) allocate(fm_thread(x_max, y_max, z_max))
    
        fm_thread(:,:,:) % phi     = 0 
        fm_thread(:,:,:) % sig_t   = 0 
        fm_thread(:,:,:) % sig_a   = 0 
        fm_thread(:,:,:) % nusig_f = 0 
        fm_thread(:,:,:) % kappa   = 0 
        do i = 1, 8
            fm_thread(:,:,:) % J0(i)  = 0 
            fm_thread(:,:,:) % J1(i)  = 0 
        end do 
	end subroutine hf_initialise_thread
	
	subroutine hf_norm(cyc) ! equivalent to NORM_FMFD()
        implicit none
        integer, intent(in):: cyc
        integer:: i, j, k, ij

        !> gather thread FMFD parameters
        do i = 1, nfm(1)
            do j = 1, nfm(2)
                do k = 1, nfm(3)
                    fm(i,j,k)%phi     = fm(i,j,k)%phi     + fm_thread(i,j,k)%phi
                    fm(i,j,k)%sig_t   = fm(i,j,k)%sig_t   + fm_thread(i,j,k)%sig_t 
                    fm(i,j,k)%sig_a   = fm(i,j,k)%sig_a   + fm_thread(i,j,k)%sig_a 
                    fm(i,j,k)%nusig_f = fm(i,j,k)%nusig_f + fm_thread(i,j,k)%nusig_f 
                    fm(i,j,k)%kappa   = fm(i,j,k)%kappa   + fm_thread(i,j,k)%kappa 
                    do mm = 1, 8
                        fm(i,j,k)%J0(mm)  = fm(i,j,k)%J0(mm)  + fm_thread(i,j,k)%J0(mm)
                        fm(i,j,k)%J1(mm)  = fm(i,j,k)%J1(mm)  + fm_thread(i,j,k)%J1(mm)
                    end do
                end do
            end do
        end do
        !print '(A,8F15.5)', 'CHK', fm(5,6,1)%J0(1:4), fm(5,6,1)%J1(5:8)
        !print '(A,8F15.5)', 'THD', fm_thread(5,6,1)%J0(1:4), fm_thread(5,6,1)%J1(5:8)
	end subroutine hf_norm
	
    subroutine hf_SET_ZERO_FLUX(mat)
        implicit none
        real(8), intent(out):: mat(:,:,:)   ! for given matrix
		
		integer :: i, j, k
		
		do i = 1, x_max
		    do j = 1, y_max
			    do k = 1, z_max
				    if (v_hf(i,j,k) /= 0.0) cycle
					mat(i,j,k) = 0.0
				end do
			end do
		end do
    end subroutine
end module hex_variables

