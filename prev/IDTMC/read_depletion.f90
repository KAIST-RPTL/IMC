
subroutine read_depletion 
	use constants, only: rd_dep
	use read_functions
	use variables
	
	implicit none
	logical :: file_exists
	integer :: Open_Error, File_Error
	character(4)::Card
	character::Card_Type	

	file_exists = .false.
	inquire(file="./inputfile/depletion.inp",exist=file_exists)
	if(file_exists==.false.) then
	  do_burn = .false.
	  return
	end if 
	
    open(unit=rd_dep,file="./inputfile/depletion.inp",status='old', action='read',iostat=Open_Error)
	Read_File : do
        read(rd_dep,*,iostat=File_Error) Card
        if (File_Error/=0) exit Read_File
        if (Card=="CARD" .or. Compare_String(Card,"card")) then
            backspace(rd_dep)
            read(rd_dep,*,iostat=File_Error) Card,Card_Type
            call Small_to_Capital(Card_Type)
            if (icore==score) print *, "depletion.inp :: CARD ", Card_Type," is being read..."
            call Read_Card(rd_dep,Card_Type)
        end if
    end do Read_File
	close(rd_dep)
	
	
	
end subroutine read_depletion















