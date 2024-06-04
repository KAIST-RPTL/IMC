module simulation_header 
    use bank_header, only : bank
    !use geometry_header, only: cell
    
    implicit none 

    ! computing time analysis
    real(8), allocatable:: t_MC(:,:), t_det(:,:), t_tot(:,:)

    ! fission source file reading
    logical:: source_read = .false.
    character(80):: source_file
    
    type particle_stack
        type(bank), pointer :: p
      contains 
        procedure :: reset => reset_stack
    end type
    !type(particle_stack), allocatable :: active_stack(:), temp_stack(:)

    contains
    
    subroutine reset_stack(this) 
        class(particle_stack) :: this 
        nullify(this%p)
    end subroutine 


end module
