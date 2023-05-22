! =============================================================================
! PROCESS_LINE extracts the first term / parameter in the line
! =============================================================================
subroutine process_line (line, ids, output) 
    implicit none 
    character(*), intent(in) :: line 
    integer :: ip, ids, jp, ide 
    character(30) :: output 
    
    jp = ids 
    do while (.true.) 
        jp = jp + 1 
        if (line(jp:jp).eq.' ') then 
            ide = jp - 1
            exit
        endif 
    enddo 
    
    output = line(ids:ide)

    ip = 0 
    do
        if (line(jp+ip:jp+ip).eq.' ' .and. (jp+ip) < len(line)) then 
            ip = ip + 1
        else
            exit
        endif 
    enddo 
    ids = jp - 1 + ip

end subroutine 
