module XS_header 
    use omp_lib
    
    implicit none 
    
    type MG_XS 
        character(20) :: mat_id
        real(8), allocatable :: sig_tr(:), sig_abs(:), sig_cap(:), sig_fis(:), nu(:), chi(:) 
        real(8), allocatable :: sig_scat(:,:)
    
      contains 
        !procedure :: find_idx => find_mat_idx
      
    end type 
    type(MG_XS), allocatable, target :: XS_MG(:), XS_MG_temp(:)
    type(MG_XS), pointer             :: XS_MG_ptr
    integer :: n_group     !> 0 if CE mode
    
	
	type MG_delayed
		real(8), allocatable :: beta(:), lambda(:), vel(:) 
		real(8), allocatable :: spectra(:,:) 
	end type 
    type(MG_delayed), allocatable :: MGD(:), MGD_temp(:)
	
	
	
    contains 
    
    function find_mat_idx (this, mat_id) result (idx) 
        type(MG_XS) :: this(:) 
        character(*) :: mat_id 
        integer :: i, idx
        
        do i = 1, size(this) 
            if (this(i)%mat_id == mat_id) then 
                idx = i 
                return 
            endif
        enddo 
        print *, "no such mat id : ", mat_id 
        stop 
        
    end function 
    

end module 
