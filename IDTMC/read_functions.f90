module read_functions 

implicit none

contains 
    subroutine Small_to_Capital(Small)
        implicit none
        character(*),intent(inout)::Small
        integer::Length,int_01        
        Length=len(Small)
        do int_01 =1, Length
            select case(Small(int_01:int_01))
            case("a") 
                Small(int_01:int_01)="A"
            case("b")  
                Small(int_01:int_01)="B"
            case("c")  
                Small(int_01:int_01)="C"
            case("d")  
                Small(int_01:int_01)="D"
            case("e")  
                Small(int_01:int_01)="E"
            case("f")  
                Small(int_01:int_01)="F"
            case("g")  
                Small(int_01:int_01)="G"
            case("h")  
                Small(int_01:int_01)="H"
            case("i")  
                Small(int_01:int_01)="I"
            case("j")  
                Small(int_01:int_01)="J"
            case("k")  
                Small(int_01:int_01)="K"
            case("l")  
                Small(int_01:int_01)="L"
            case("m")  
                Small(int_01:int_01)="M"
            case("n")  
                Small(int_01:int_01)="N"
            case("o")  
                Small(int_01:int_01)="O"
            case("p")  
                Small(int_01:int_01)="P"
            case("q")  
                Small(int_01:int_01)="Q"
            case("r")  
                Small(int_01:int_01)="R"
            case("s")  
                Small(int_01:int_01)="S"
            case("t")  
                Small(int_01:int_01)="T"
            case("u")  
                Small(int_01:int_01)="U"
            case("v")  
                Small(int_01:int_01)="V"
            case("w")  
                Small(int_01:int_01)="W"
            case("x")  
                Small(int_01:int_01)="X"
            case("y")  
                Small(int_01:int_01)="Y"
            case("z")  
                Small(int_01:int_01)="Z"
            end select
        end do
    end subroutine
    !忙式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式忖!
    !弛           01_02. Compare_String           弛!
    !戌式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式戎!
    function Compare_String(Input_01,Input_02)
        implicit none
        logical::Compare_String
        logical,dimension(:),allocatable::Logical_String
        character(*)::Input_01,Input_02
        character,dimension(:),allocatable::S_01,S_02
        character,dimension(:),allocatable::P_01,P_02
        character::C_01,C_02
        integer::i,j,N_01,N_02,T_01,T_02
        integer::Identifier,Max
        ! 01_02_01. Initialization
        ! 01_02_01_01. Obtain Length Information
        N_01=len(Input_01)
        N_02=len(Input_02)
        ! 01_02_01_02. Remove blank space
        Identifier=0
        do i=1,N_01
            if (Input_01(i:i)/=" " .or. Input_01(i:i)/="    ") then 
                    Exit
            end if
            if (Input_01(i:i)==" " .or. Input_01(i:i)=="    ") then 
                Identifier=Identifier+1                
            end if
        end do
        T_01=Identifier
        Identifier=0
        do i=1,N_02
            if (Input_02(i:i)/=" " .or. Input_02(i:i)/="    ") then 
                    Exit
            end if
            if (Input_02(i:i)==" " .or. Input_02(i:i)=="    ") then 
                Identifier=Identifier+1                
            end if
        end do
        T_02=Identifier
        ! 01_02_01_03. Find Maximum Value
        if ( N_01-T_01>N_02-T_02) then
            Max=N_01-T_01
        else
            Max=N_02-T_02
        end if  
        ! 01_02_01_04. Allocate S_01,S_02 and Logical_String
        Identifier=0
        allocate(Logical_String(1:Max))
        Compare_String=.false.
        do j=1,Max
            Logical_String(j)=.false.
        end do     
        allocate(S_01(1:Max))
        allocate(S_02(1:Max))
        ! 01_02_02. Comparison     
        do i=1,Max
            S_01(i:i)="|"
            S_02(i:i)="|"
            if(i<=N_01-T_01) then
                if (Input_01(i+T_01:i+T_01)/=" " .or. Input_01(i+T_01:i+T_01)/="    ") then 
                    S_01(i:i)=Input_01(i+T_01:i+T_01)
                    call Small_to_Capital(S_01(i))
                end if
            end if
            if(i<=N_02-T_02) then
                if (Input_02(i+T_02:i+T_02)/=" " .or. Input_02(i+T_02:i+T_02)/="    ") then 
                    S_02(i:i)=Input_02(i+T_02:i+T_02)
                    call Small_to_Capital(S_02(i))
                end if
            end if
        end do
        ! 01_02_03. Calculate Compare_String  
        do i=1,Max            
            if(S_01(i)==S_02(i)) then
                Logical_String(i)=.true.
            end if            
        end do
         do j=1,Max
            if(Logical_String(j)==.true.) then
                Identifier=Identifier+1
            end if     
        end do
        if (Identifier==Max) then
                Compare_String=.true.
        end if       
    end function
    !忙式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式忖!
    !弛         01_03. Card Error Print           弛!
    !戌式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式式戎!
    subroutine Card_Error(Card_Type,Char_Temp)
    implicit none
    character(*),intent(in)::Card_Type,Char_Temp
        if(Card_Type==Char_Temp) then
            write(*,'(2a)') " Error in CARD ", trim(Card_Type)
        else
            write(*,'(4a)') " Error in CARD ", trim(Card_Type),", ",trim(Char_Temp)
        end if    
        stop
    end subroutine


end module
