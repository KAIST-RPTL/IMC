
subroutine exception_handler( subroutine_name, stop_reason )
use variables, only : icore
implicit none
character(len=*), intent(in) :: subroutine_name
character(len=*), intent(in) :: stop_reason

print *, "subroutine name :: ", subroutine_name
print *, "stop reason :: ", stop_reason
print *, "core index :: ", icore
stop

end subroutine
