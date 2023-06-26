module hex_geometry
	
	use FMFD_header, only: ncm, fmfdon, fm_thread
	use hex_variables
    implicit none
	
    contains
	
	! note: assumes row-aligned hex grid; for a column-aligned grid, flip the first two indices
	function hex_coords(cart_coords, ori, sca)
	    real(8), intent(in) :: cart_coords(3) ! cartesian coordinates
		real(8), intent(in) :: ori(3) ! centre of the (0, 0) hex
		real(8), intent(in) :: sca(3) ! scale; distance from hex centre to next
		integer :: hex_coords(3) ! hex coordinate indices
		
		real(8) :: norm_coords(3) ! normalised cartesian coordinates
		integer :: i, j ! temporary hex coordinates
		real(8) :: x, y ! local normalised cartesian coordinates, relative to the centre of hex (i, j)
		
		! normalise
		norm_coords(:) = (cart_coords(:) - ori(:)) / sca(:)
		
		! check which box the coordinate falls into
		j = ceiling(norm_coords(2) / sqrt(3.0/4.0) - (2.0/3.0))
		i = ceiling(norm_coords(1) - 0.5 * real(1 + j))
		
		! calculate local coordinates
		x = norm_coords(1) - real(i) - 0.5 * real(j)
		y = norm_coords(2) - sqrt(3.0/4.0) * real(j)
		
		! use local coordinates to check if the coordinate falls in a discrepancy between the box and hex schemes
		if (y > (1.0 - x) / sqrt(3.0)) then
		    j = j + 1
		else if (y > (1.0 + x) / sqrt(3.0)) then
		    j = j + 1
			i = i - 1
		end if
		hex_coords(1) = i
		hex_coords(2) = j
		
		! vertical coordinate
		hex_coords(3) = ceiling(norm_coords(3) - 0.5)
	end function hex_coords
	
	! note: assumes column-aligned CMFD hex grid
	function hc_cmfd_coords(p_coords, ori, sca_cmfd, n_cmfd)
	    real(8), intent(in) :: p_coords(3)
		real(8), intent(in) :: ori(3)
		real(8), intent(in) :: sca_cmfd(3)
		integer, intent(in) :: n_cmfd
		integer :: hc_cmfd_coords(3)
		
		real(8) :: adj_car(3)
		real(8) :: adj_ori(3)
		integer :: adj_hex(3)
		
		adj_car = [p_coords(2), p_coords(1), p_coords(3)]
		adj_ori = [ori(2), ori(1), ori(3)]
		adj_hex = hex_coords(adj_car, adj_ori, sca_cmfd)
		hc_cmfd_coords = [adj_hex(2), adj_hex(1), adj_hex(3)]
	end function hc_cmfd_coords
	
	! note: assumes column-aligned CMFD hex grid and row-aligned FMFD hex grid
	! note: for the moment assumes sca_cmfd(3) == sca_fmfd(3), that is to say 1 axial fine mesh per coarse mesh
    function hf_fmfd_coords(p_coords, ori, sca_cmfd, sca_fmfd, n_fmfd)
	    real(8), intent(in) :: p_coords(3) ! cartesian coordinate of the particle
		real(8), intent(in) :: ori(3) ! coordinates of the centre of the centre CMFD node
		real(8), intent(in) :: sca_cmfd(3) ! size of the CMFD cells
		real(8), intent(in) :: sca_fmfd(3) ! size of the FMFD cells
		integer, intent(in) :: n_fmfd ! number of FMFD cells along each edge of CMFD cells
		integer :: hf_fmfd_coords(9)
		
		integer :: boundaries(8) ! 0 for normal, 1 if out of boundary
		real(8) :: adj_ori(3)
		real(8) :: origin(3)
		
		! temporary variables for coordinate finding
	    real(8) :: cm_car(3) ! x and y reversed
	    integer :: cm_hex(3) ! CM coordinates, x and y reversed
	    integer :: cm_cor(3) ! corrected CM coordinates
	    real(8) :: fm_car(3) ! particle coordinates relative to centre of current coarse mesh node
	    integer :: fm_cor(3) ! local FM coordinates in a given coarse mesh node
		
		integer :: i, j, k, ii, jj, kk
		
		cm_car = [p_coords(2), p_coords(1), p_coords(3)]
		adj_ori = [ori(2), ori(1), ori(3)]
		cm_hex = hex_coords(cm_car, adj_ori, sca_cmfd)
		cm_cor = [cm_hex(2), cm_hex(1), cm_hex(3)]
		fm_car(1) = p_coords(1) - ori(1) - sqrt(3.0/4.0) * cm_cor(1) * sca_cmfd(1)
		fm_car(2) = p_coords(2) - ori(2) - (cm_cor(2) + 0.5 * real(cm_cor(1))) * sca_cmfd(2)
		fm_car(3) = p_coords(3) - ori(3) - cm_cor(3) * sca_cmfd(3)
		origin = [0.0, 0.0, 0.0]
		fm_cor = hex_coords(fm_car, origin, sca_fmfd)
		! adjust coordinates; the central fine mesh cell of the FMFD grid is (n_fmfd, n_fmfd)!
		fm_cor(1) = fm_cor(1) + n_fmfd
		fm_cor(2) = fm_cor(2) + n_fmfd
		! (TODO): implement multiple fine mesh cells across the axial height of each coarse-mesh cell
		boundaries(:) = 0
		if (fm_cor(1) + fm_cor(2) < n_fmfd + 1)    boundaries(1) = 1
		if (fm_cor(1) + fm_cor(2) .ge. 3 * n_fmfd) boundaries(8) = 1
		if (fm_cor(1) < 1)                         boundaries(2) = 1
		if (fm_cor(1) > 2 * n_fmfd - 1)            boundaries(7) = 1
		if (fm_cor(2) < 1)                         boundaries(3) = 1
		if (fm_cor(2) > 2 * n_fmfd - 1)            boundaries(6) = 1
		do while (sum(boundaries) .ne. 0)
		    if (boundaries(1) == 1) then
		        if (boundaries(2) == 1) then
		    	    fm_cor(1) = 1
		    		fm_cor(2) = n_fmfd
					exit
		    	else
		    	    fm_car(1) = fm_car(1) + (1.0 + 0.5) * sca_fmfd(1) / 3.0
					fm_car(2) = fm_car(2) + sqrt(3.0/4.0) * sca_fmfd(2) / 3.0
		    	end if
		    end if
		    if (boundaries(2) == 1) then
		        if (boundaries(6) == 1) then
		    	    fm_cor(1) = 1
		    		fm_cor(2) = 2 * n_fmfd - 1
					exit
		    	else
		    	    fm_car(1) = fm_car(1) + (2.0 - 0.5) * sca_fmfd(1) / 3.0
					fm_car(2) = fm_car(2) - sqrt(3.0/4.0) * sca_fmfd(2) / 3.0
		    	end if
		    end if
			if (boundaries(6) == 1) then
			    if (boundaries(8) == 1) then
				    fm_cor(1) = n_fmfd
					fm_cor(2) = 2 * n_fmfd - 1
					exit
				else
				    fm_car(2) = fm_car(2) - 2.0 * sqrt(3.0/4.0) * sca_fmfd(2) / 3.0
				end if
			end if
			if (boundaries(8) == 1) then
			    if (boundaries(7) == 1) then
				    fm_cor(1) = 2 * n_fmfd - 1
					fm_cor(2) = n_fmfd
					exit
				else
				    fm_car(1) = fm_car(1) - (1.0 + 0.5) * sca_fmfd(1) / 3.0
					fm_car(2) = fm_car(2) - sqrt(3.0/4.0) * sca_fmfd(2) / 3.0
				end if
			end if
			if (boundaries(7) == 1) then
			    if (boundaries(3) == 1) then
				    fm_cor(1) = 2 * n_fmfd - 1
					fm_cor(2) = 1
					exit
				else
				    fm_car(1) = fm_car(1) - (2.0 - 0.5) * sca_fmfd(1) / 3.0
					fm_car(2) = fm_car(2) + sqrt(3.0/4.0) * sca_fmfd(2) / 3.0
				end if
			end if
			if (boundaries(3) == 1) then
			    if (boundaries(1) == 1) then
				    fm_cor(1) = n_fmfd
					fm_cor(2) = 1
					exit
				else
				    fm_car(2) = fm_car(2) + 2.0 * sqrt(3.0/4.0) * sca_fmfd(2) / 3.0
				end if
			end if
			fm_cor = hex_coords(fm_car, origin, sca_fmfd)
			fm_cor(1) = fm_cor(1) + n_fmfd
			fm_cor(2) = fm_cor(2) + n_fmfd
		    boundaries(:) = 0
		    if (fm_cor(1) + fm_cor(2) < n_fmfd + 1)    boundaries(1) = 1
		    if (fm_cor(1) + fm_cor(2) .ge. 3 * n_fmfd) boundaries(8) = 1
		    if (fm_cor(1) < 1)                         boundaries(2) = 1
		    if (fm_cor(1) > 2 * n_fmfd - 1)            boundaries(7) = 1
		    if (fm_cor(2) < 1)                         boundaries(3) = 1
		    if (fm_cor(2) > 2 * n_fmfd - 1)            boundaries(6) = 1
		end do
		
		i = cm_cor(1)
		j = cm_cor(2)
		k = cm_cor(3)
		ii = fm_cor(1)
		jj = fm_cor(2)
		kk = 1 ! (TODO) implement multiple fine-mesh cells per coarse-mesh cell!
		
		hf_fmfd_coords(1) = i
		hf_fmfd_coords(2) = j
		hf_fmfd_coords(3) = k
		hf_fmfd_coords(4) = ii
		hf_fmfd_coords(5) = jj
		hf_fmfd_coords(6) = kk
		hf_fmfd_coords(7) = (fcr - 1) * (ncm(2) - 1) + ii + (i - 1) * fcr + (j - 1) * (1 - fcr)
		hf_fmfd_coords(8) = jj + (i - 1) * (fcr - 1) + (j - 1) * (2 * fcr - 1)
		hf_fmfd_coords(9) = kk + (k - 1) * fcz
	end function hf_fmfd_coords
	
	subroutine hex_surf(p_coords, dir, ori, sca_cmfd, sca_fmfd, n_fmfd, i_xyz, d_mesh, i_surf)
	    real(8), intent(in) :: p_coords(3) ! cartesian coordinate of the particle
		real(8), intent(in) :: dir(3) ! unit vector showing the direction of the particle
		real(8), intent(in) :: ori(3) ! coordinates of the centre of the centre CMFD node
		real(8), intent(in) :: sca_cmfd(3) ! size of the CMFD cells
		real(8), intent(in) :: sca_fmfd(3) ! size of the FMFD cells
		integer, intent(in) :: n_fmfd ! number of FMFD cells along each edge of CMFD cells
        integer, intent(inout) :: i_xyz(3)
        real(8), intent(inout) :: d_mesh
        integer, intent(inout) :: i_surf
		
		! temporary variables for coordinate finding
	    integer :: cm_cor(3) ! corrected CM coordinates
	    real(8) :: fm_car(3) ! particle coordinates relative to centre of current coarse mesh node
	    integer :: fm_cor(3) ! local FM coordinates in a given coarse mesh node
		
		integer :: boundaries(8) ! 0 for normal, 1 if out of boundary
		integer :: hf_global(3) ! global fine-mesh coordinates
		
		integer :: inlet(9) ! receiving hf_fmfd_coords return value
		
		real(8) :: disp(3) ! normalised displacement relative to node centre
		real(8) :: vel(3) ! normalised direction vector
		
		real(8) :: dist(8) ! normalised distance to relevant cell boundary
		real(8) :: rat ! ratio between dduct and sca_fmfd(1)
		
		integer :: i
		
		!! SAME ALGORITHM AS HF_FMFD_COORDS TO OBTAIN FM_CAR AND FM_COR
		inlet = hf_fmfd_coords(p_coords, ori, sca_cmfd, sca_fmfd, n_fmfd)
		cm_cor = inlet(1:3)
		fm_cor = inlet(4:6)
		hf_global = inlet(7:9)
		! reset fm_car to its original values
		fm_car(1) = p_coords(1) - ori(1) - sqrt(3.0/4.0) * cm_cor(1) * sca_cmfd(1)
		fm_car(2) = p_coords(2) - ori(2) - (cm_cor(2) + 0.5 * real(cm_cor(1))) * sca_cmfd(2)
		fm_car(3) = p_coords(3) - ori(3) - cm_cor(3) * sca_cmfd(3)
		! get normalised fine-mesh relative cartesian coordinates
		disp(1) = (fm_car(1) / sca_fmfd(1)) - real(fm_cor(1) - n_fmfd) - 0.5 * real(fm_cor(2) - n_fmfd)
		disp(2) = (fm_car(2) / sca_fmfd(2)) - sqrt(3.0/4.0) * real(fm_cor(2) - n_fmfd)
		disp(3) = (fm_car(3) / sca_fmfd(3)) ! (TODO): implement multiple axial fine-mesh cells per coarse-mesh cell!
		vel(:) = dir(:) / sca_fmfd(:)
		
		!! CHECK FOR NONEXISTENT OR ABNORMAL FINE MESH CELL BOUNDARIES
		boundaries(:) = 1 ! 1 internal, 0 coarse mesh boundary, 2 nonexistent
		! axial boundaries always exist and are normal
		if (fm_cor(1) + fm_cor(2) == n_fmfd + 1) then
		    boundaries(1) = boundaries(1) * 0
			boundaries(3) = boundaries(3) * 2
		end if
		if (fm_cor(2) == 1) then
		    boundaries(3) = boundaries(3) * 0
			boundaries(7) = boundaries(7) * 2
		end if
		if (fm_cor(1) == 2 * n_fmfd - 1) then
		    boundaries(7) = boundaries(7) * 0
			boundaries(8) = boundaries(8) * 2
		end if
		if (fm_cor(1) + fm_cor(2) == 3 * n_fmfd - 1) then
		    boundaries(8) = boundaries(8) * 0
			boundaries(6) = boundaries(6) * 2
		end if
		if (fm_cor(2) == 2 * n_fmfd - 1) then
		    boundaries(6) = boundaries(6) * 0
			boundaries(2) = boundaries(2) * 2
		end if
		if (fm_cor(1) == 1) then
		    boundaries(2) = boundaries(2) * 0
			boundaries(1) = boundaries(1) * 2
		end if
		
		!! TEST BOUNDARIES
		rat = dduct / sca_fmfd(1)
		dist(:) = 1.0
		if (boundaries(1) == 1) dist(1) = (1.0 + 2.0 * disp(1)) / (2.0 * vel(1))
		if (boundaries(2) == 1) dist(2) = (1.0 + disp(1) - sqrt(3.0) * disp(2)) / (vel(1) - sqrt(3.0) * vel(2))
		if (boundaries(3) == 1) dist(3) = (1.0 + disp(1) + sqrt(3.0) * disp(2)) / (vel(1) + sqrt(3.0) * vel(2))
		if (boundaries(4) == 1) dist(4) = (1.0 + 2.0 * disp(3)) / (2.0 * vel(3))
		if (boundaries(5) == 1) dist(5) = (1.0 - 2.0 * disp(3)) / (-2.0 * vel(3))
		if (boundaries(6) == 1) dist(6) = (1.0 - disp(1) - sqrt(3.0) * disp(2)) / (-vel(1) - sqrt(3.0) * vel(2))
		if (boundaries(7) == 1) dist(7) = (1.0 - disp(1) + sqrt(3.0) * disp(2)) / (-vel(1) + sqrt(3.0) * vel(2))
		if (boundaries(8) == 1) dist(8) = (1.0 - 2.0 * disp(1)) / (-2.0 * vel(1))
		! assembly boundaries
		if (boundaries(1) == 0) dist(1) = (2.0 * rat + sqrt(3.0) * disp(1) + disp(2)) / (sqrt(3.0) * vel(1) + vel(2))
		if (boundaries(2) == 0) dist(2) = (2.0 * rat + sqrt(3.0) * disp(1) - disp(2)) / (sqrt(3.0) * vel(1) - vel(2))
		if (boundaries(3) == 0) dist(3) = (2.0 * rat + 2.0 * disp(2)) / (2.0 * vel(2))
		if (boundaries(6) == 0) dist(6) = (2.0 * rat - 2.0 * disp(2)) / (-2.0 * vel(2))
		if (boundaries(7) == 0) dist(7) = (2.0 * rat - sqrt(3.0) * disp(1) + disp(2)) / (-sqrt(3.0) * vel(1) + vel(2))
		if (boundaries(8) == 0) dist(8) = (2.0 * rat - sqrt(3.0) * disp(1) - disp(2)) / (-sqrt(3.0) * vel(1) - vel(2))
		dist(:) = 0.0 - dist(:)
		
		!! OUTPUT RESULTS
		i_surf = -1
		do i = 1, 8
		    if ((dist(i) < d_mesh) .and. (dist(i) > 0.0)) then
			    d_mesh = dist(i) ! find the shortest positive distance to a surface
				i_surf = i
			end if
		end do
		if (i_surf == -1 .or. d_mesh .le. 0.0) then
		    print *, "SURFACE FINDING ERROR"
		    print *, p_coords
		    print *, ori
		    print *, cm_cor
		    print *, fm_car
		    print *, fm_cor
			print *, boundaries
		    print *, disp
		    print *, vel
		    print *, i_surf, d_mesh
			stop
		end if
		
!        if(hf_global(1)==5 .and. hf_global(2)==5 .and. (i_surf==1) .and. p_coords(2)+d_mesh*dir(2)>-5) then
!            print '(A,3F12.5)', 'POSI', p_coords(1:3)
!            print '(A,3F12.5)', 'UVW ', dir(1:3)
!            print '(A,3F12.5)', 'DISP', disp(1:3)
!            print '(A,3F12.5)', 'VELO', vel(1:3)
!            print '(A,8I2)', 'BDS', boundaries(1:8)
!            print '(A,8F12.5,I2)', 'DIST', dist(1:8), i_surf
!            print '(A,3F15.8,I2)', 'ARR', p_coords(1:3) + d_mesh * dir(1:3), i_surf
!            !print *, ' '
!        endif
		!if (i_surf /= 4 .and. i_surf /= 5) then
		!    print *, p_coords + d_mesh * dir
		!end if
		
		i_xyz = hf_global
	end subroutine hex_surf
	
	! tallying function for partial flux
	subroutine hf_surf(xyz, surf, wgt)
	    integer, intent(in) :: xyz(3)
		integer, intent(in) :: surf
		real(8), intent(in) :: wgt
		
        if ( .not. fmfdon ) return
		
		select case(surf)
		case (1, 2, 3, 4)
		    fm_thread(xyz(1), xyz(2), xyz(3))%J0(surf) = fm_thread(xyz(1), xyz(2), xyz(3))%J0(surf) + wgt
		case (5, 6, 7, 8)
		    fm_thread(xyz(1), xyz(2), xyz(3))%J1(surf) = fm_thread(xyz(1), xyz(2), xyz(3))%J1(surf) + wgt
		end select
	end subroutine hf_surf
end module hex_geometry
