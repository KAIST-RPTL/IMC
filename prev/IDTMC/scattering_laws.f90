module scattering_laws 

use constants, only : wt_coord
use particle_header 
use ace_header, only : EnergyDistDataForm, ace
use randoms, only : rang

implicit none

contains 

subroutine law_selector (erg, iso, mu, dist, iMT, law)
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    real(8), intent(inout) :: mu
    integer, intent(in) :: iMT
    integer, intent(in) :: law
    
    !print '(A16, I3, A8, I3, A5, I3)', 'selected law :: ', law, ' for MT ', ace(iso)%MT(iMT), ' TYR ', ace(iso)%TY(iMT)
    select case (law)
    case (1)
        call LAW1(erg, iso, mu,  dist, iMT) 
    case (3)  
        call LAW3 (erg, iso, mu, dist, iMT)
    case (4)  
        call LAW4 (erg, iso, mu, dist, iMT) 
    case (44) 
        call LAW44(erg, iso, mu, dist, iMT)  
    case (61) 
        call LAW61(erg, iso, mu, dist, iMT)  
        
    
    case (7) 
        call LAW7 (erg, iso, mu, dist, iMT)  
    case (9)                       
        call LAW9 (erg, iso, mu, dist, iMT)  
    case (11)                      
        call LAW11(erg, iso, mu, dist, iMT)  
    
    case (66)
        call LAW66(erg, iso, mu, dist, iMT)  
        print *, "WARNING :: LAW66 is called"
        
    !> Unprepared Laws 
    case (2) 
        print *, "LAW2 (Discrete photon energy) is not prepared"
        stop 
    case (5) 
        print *, "LAW5 (General evaporation spectrum) is not prepared"
        stop 
    case (22) 
        print *, "LAW22 (UK-law2; Tabular linear functions of incident energy out) is not prepared"
        stop 
    case (24) 
        print *, "LAW24 (UK-law6; Equiprobable energy multiplier) is not prepared"
        stop 
    case default  
        print *, "UNKNOWN LAW", law
        stop
    end select
    
    
end subroutine



! ========================================================= !
!  >>>>>  law 1 (from endf Law 1) 
!        -- tabular equiprobable energy bins.
! ========================================================= !
subroutine LAW1 (erg, iso, mu, dist, iMT)
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    real(8), intent(inout) :: mu
    integer, intent(in) :: iMT

    integer :: i, k
    integer :: pt1, pt2, pt3 
    integer :: NR, NE, NET 
    real(8) :: ipfac
    real(8), allocatable :: Ein(:), Eout1(:), Eout2(:) 
    real(8) :: E_1, E_K
    real(8) :: E_, Eout

    NR  = dist%LDAT(1)
    if (NR /= 0) then 
        print *, "WARNING :: NBT exists (law 1)"
        !allocate(NBT(NR)); allocate(INTP(NR)) 
        !NBT (1:NR) = dist%LDAT(2   :NR+1  ) 
        !INTP(1:NR) = dist%LDAT(NR+2:2*NR+1) 
    endif

    NE  = dist%LDAT(2+2*NR)
    NET = dist%LDAT(3+2*NR+NE) 
    allocate(Ein(1:NE)) 
    allocate(Eout1(1:NET))
    allocate(Eout2(1:NET))


    Ein(1:NE) = dist%LDAT(3+2*NR:3+2*NR+NE-1)
    pt1 = 1; pt2 = NE 
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= Ein(pt3)) then 
            pt1 = pt3 
        else 
            pt2 = pt3 
        endif 
    enddo
    ipfac = max(0.d0, min(1.d0,(erg-Ein(pt1))/(Ein(pt2)-Ein(pt1))))
    
    k = rang()*NET + 1

    Eout1(1:NET) = dist%LDAT(4+2*NR+NE+(pt1-1)*NET:4+2*NR+NE+(pt1)*NET-1) 
    Eout2(1:NET) = dist%LDAT(4+2*NR+NE+(pt2-1)*NET:4+2*NR+NE+(pt2)*NET-1) 

    E_1 = Eout1(1)+ipfac*(Eout2(1)-Eout1(1))
    E_K = Eout1(NET) + ipfac*(Eout2(NET)-Eout1(NET)) 

    !> randomly select energy i or i+1 
    if(rang() > ipfac) then 
        E_ = Eout1(k) + rang()*(Eout1(k+1)-Eout1(k))
        Eout = E_1 + (E_ - Eout1(1))*(E_K-E_1)/(Eout1(NET)-Eout1(1))
    else 
        E_ = Eout2(k) + rang()*(Eout2(k+1)-Eout2(k))
        Eout = E_1 + (E_ - Eout2(1))*(E_K-E_1)/(Eout2(NET)-Eout2(1))
    endif 

    !> Check CM or LAB 
    !if (ace(iso)%TY(iMT)<0) then
    !    Eout = Eout + (erg + 2*mu*(ace(iso)%atn+1)*sqrt(erg * Eout))/(ace(iso)%atn+1)**2
    !endif
    erg = Eout


end subroutine 


! ========================================================= !
!  >>>>>  law 3 & 33 -- level scattering.
! law 3 applies only for neutron in-neutron out scattering.
! law 33 for a more general combination of particle types.
! ========================================================= !
subroutine LAW3 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso
    real(8), intent(in) :: mu
    integer, intent(in) :: iMT
    integer :: i 
    real(8) :: E_in, E_outCM, E_outLAB
    real(8) :: temp
    
    
    E_outCM  = dist%LDAT(2)*(erg - dist%LDAT(1)) 
    if (ace(iso)%TY(iMT)<0) then 
        erg         = E_outCM
    else 
        temp      = erg + 2*mu*(ace(iso)%atn +1)*sqrt(erg * E_outCM)
        E_outLAB = E_outCM + temp/(ace(iso)%atn+1)**2
        erg = E_outLAB
    endif
    
    !print *, 'law 3 ::', ace(iso)%library, ace(iso)%TY(iMT)
    
end subroutine 



! ========================================================= !
!  >>>>>  law 4 (from endf law 1) 
!        -- continuous erg tabular distribution.
! ========================================================= !
subroutine LAW4 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    integer :: i, j
    logical :: found
    
    integer :: NR
    integer :: NE 
    integer, allocatable :: NP(:), K(:) 
    real(8), allocatable, dimension(:) :: NBT, INTP
    real(8), allocatable, dimension(:) :: E, L 
    !real(8), allocatable, dimension(:) :: EOUT, PDF, CDF 
    
    integer :: pt1, pt2, pt3 
    real(8) :: rn
    real(8) :: ipfac
    
    integer :: INTT_, INTT, ND
    real(8) :: c1, c2, c3, c4   ! c(i,k), c(i+1,k), c(i,k+1), c(i+1,k+1) 
                                ! respectively in page 2-43 ENDF Manual Vol.1
    real(8) :: p1, p2            ! p(l,k), p(l,k+1)
    real(8) :: E1, E2            ! E(l,k), E(l,k+1) 
    real(8) :: Eout, E_, temp
    
    NR = dist%LDAT(1)
    if (NR /= 0) then 
        print *, "WARNING :: NBT exists (law 4)"
        allocate(NBT(NR)); allocate(INTP(NR)) 
        NBT (1:NR) = dist%LDAT(2   :NR+1  ) 
        INTP(1:NR) = dist%LDAT(NR+2:2*NR+1) 
    endif
    
    NE = dist%LDAT(2+2*NR)
    allocate(E(1:NE)); E(1:NE) = dist%LDAT(3+2*NR   :3+2*NR+NE-1  ) 
    allocate(L(NE)); L(1:NE) = dist%LDAT(3+2*NR+NE:2+2*NR+2*NE) 
    
    !> Determine the incident Energy interval
    pt1 = 1; pt2 = NE
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= E(pt3) ) then 
            pt1 = pt3 
        else 
            pt2 = pt3
        endif 
    enddo 
    ipfac = max(0.d0, min(1.d0,(erg-E(pt1))/(E(pt1+1)-E(pt1))))
    
    !print *, 'points ',pt1, pt2, E(pt1)
    
    !> Sample the outgoing Energy     
    allocate(K(1:NE)); allocate(NP(1:NE))
    K(1)  = 3+2*NR+2*NE
    NP(1) = dist%LDAT(K(1)+1)
    do i = 2, NE 
        K(i)  = K(i-1)+2+3*NP(i-1)
        NP(i) = dist%LDAT(K(i)+1)
    enddo 
    
    INTT_ = dist%LDAT(K(pt1))
    INTT  = mod(INTT_, 10)
    ND      = (INTT_ - INTT)/10
    
    
    rn = rang(); found = .false.
    do i = 1, ND
        if(found) exit 
        c1 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i-1)
        c2 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i-1)
        c3 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i)
        c4 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i)
        
        if (rn > c1+ipfac*(c2-c1) .and. rn < c3+ipfac*(c4-c3)) then 
            Eout  = dist%LDAT(K(pt1)+2+i-1) + ipfac*(dist%LDAT(K(pt2)+2+i-1) - dist%LDAT(K(pt1)+2+i-1))
            found = .true. 
        endif 
        
    enddo 
    
    pt3 = pt1
    if (rang() < ipfac) pt3 = pt2
    do i = ND+1, NP(pt3)
        if(found) exit 
        
        !print *, dist%LDAT(K(pt1)+2+i-1),  dist%LDAT(K(pt1)+2+2*NP(pt1)+i-1)

        c1 = dist%LDAT(K(pt3)+2+2*NP(pt3)+i-1)
        c3 = dist%LDAT(K(pt3)+2+2*NP(pt3)+i)        
        
        if (rn > c1 .and. rn < c3) then 
            p1 = dist%LDAT(K(pt3)+2+NP(pt3)+i-1)
            p2 = dist%LDAT(K(pt3)+2+NP(pt3)+i)
            E1 = dist%LDAT(K(pt3)+2+i-1)
            E2 = dist%LDAT(K(pt3)+2+i)
            
            
            if (INTT == 1 .or. p1==p2) then  !> histogram
                E_ = E1 + (rn - c1)/p1
            elseif(INTT == 2) then      !> linear-linear
                temp = (p2-p1)/(E2-E1)
                E_ = E1 + (sqrt(p1**2 + 2*temp*(rn-c1)) -p1)/temp
            else 
                print *, 'ERROR :: LAW4, INTT NUMBER IS WRONG AS ', INTT
            endif 
            
            !pt3 = pt1  !> pt3 is for l
            !if (pt1 == pt2) pt1 = pt1-1
            E1 = dist%LDAT(K(pt1)+2) + ipfac*(dist%LDAT(K(pt1+1)+2)-dist%LDAT(K(pt1)+2))
            E2 = dist%LDAT(K(pt1)+2+NP(pt1)-1) + ipfac*(dist%LDAT(K(pt1+1)+2+NP(pt1+1)-1) - dist%LDAT(K(pt1)+2+NP(pt1)-1))
            !E2 = dist%LDAT(K(pt1)+2+i-1) + ipfac*(dist%LDAT(K(pt1+1)+2+i-1) - dist%LDAT(K(pt1)+2+i-1))
            
            Eout = E1 + (E_-dist%LDAT(K(pt3)+2))*(E2-E1)/(dist%LDAT(K(pt3)+2+NP(pt3)-1) - dist%LDAT(K(pt3)+2))
            !Eout = E1 + (E_-dist%LDAT(K(pt3)+2))*(E2-E1)/(dist%LDAT(K(pt3)+2+i-1) - dist%LDAT(K(pt3)+2))
            
            found = .true. 
        endif 
    enddo 
    !if(found /= .true.) then 
    !    print *, "ERROR :: now found"
    !    stop
    !endif 
    !> Check CM or LAB 
    !if (ace(iso)%TY(iMT) < 0) then
    !    Eout = Eout + (erg + 2*mu*(ace(iso)%atn+1)*sqrt(erg * Eout))/(ace(iso)%atn+1)**2
    !endif
    erg = Eout
        
end subroutine 

! ========================================================= !
!  >>>>>  law 44 (from endf law 1) 
!        -- kalbach-87 correlated formalism.
! ========================================================= !
subroutine LAW44 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    integer :: i, j
    logical :: found
    
    integer :: NR
    integer :: NE 
    integer, allocatable :: NP(:), K(:) 
    real(8), allocatable, dimension(:) :: NBT, INTP
    real(8), allocatable, dimension(:) :: E!, I
    !real(8), allocatable, dimension(:) :: EOUT, PDF, CDF 
    
    integer :: pt1, pt2, pt3 
    real(8) :: rn
    real(8) :: ipfac
    
    integer :: INTT_, INTT, ND
    real(8) :: c1, c2, c3, c4   ! c(i,k), c(i+1,k), c(i,k+1), c(i+1,k+1) 
                                ! respectively in page 2-43 ENDF Manual Vol.1
    real(8) :: p1, p2            ! p(l,k), p(l,k+1)
    real(8) :: E1, E2            ! E(l,k), E(l,k+1) 
    real(8) :: A, R, T
    real(8) :: Eout, E_, temp
    
    
    
    NR = dist%LDAT(1)
    if (NR /= 0) then 
        print *, "WARNING :: NBT exists (law 44)"
        allocate(NBT(NR)); allocate(INTP(NR)) 
        NBT (1:NR) = dist%LDAT(2   :NR+1  ) 
        INTP(1:NR) = dist%LDAT(NR+2:2*NR+1) 
    endif
    
    
    NE = dist%LDAT(2+2*NR)
    allocate(E(1:NE)); E(1:NE) = dist%LDAT(3+2*NR   :3+2*NR+NE-1  ) 
    !allocate(I(NE)); I(1:NE) = dist%LDAT(3+2*NR+NE:2+2*NR+2*NE) 
    !> Determine the incident Energy interval
    pt1 = 1; pt2 = NE
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= E(pt3) ) then 
            pt1 = pt3 
        else 
            pt2 = pt3
        endif 
    enddo 
    ipfac = max(0.d0, min(1.d0,(erg-E(pt1))/(E(pt1+1)-E(pt1))))
    
    !> Sample the outgoing Energy
    allocate(K(1:NE)); allocate(NP(1:NE))
    K(1)  = 3+2*NR+2*NE
    NP(1) = dist%LDAT(K(1)+1)
    do i = 2, NE 
        K(i)  = K(i-1)+2+5*NP(i-1)
        NP(i) = dist%LDAT(K(i)+1)
    enddo 
    
    
    INTT_ = dist%LDAT(K(pt1))
    INTT  = mod(INTT_, 10)
    ND      = (INTT_ - INTT)/10
    
    rn = rang(); found = .false.
    do i = 1, ND              !> discrete spectra 
        if(found) exit 
        c1 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i-1)
        c2 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i-1)
        c3 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i)
        c4 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i)
        
        if (rn > c1+ipfac*(c2-c1) .and. rn < c3+ipfac*(c4-c3)) then 
            Eout  = dist%LDAT(K(pt1)+2+i-1) + ipfac*(dist%LDAT(K(pt2)+2+i-1) - dist%LDAT(K(pt1)+2+i-1))
            R = dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1) &
                + ipfac*(dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1) - dist%LDAT(K(pt2)+2+3*NP(pt2)+i-1))
            A = dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1) &
                + ipfac*(dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1) - dist%LDAT(K(pt2)+2+4*NP(pt2)+i-1))
            found = .true. 
        endif 
    enddo 
    
    !print *, dist%LDAT(K(pt1)+2+3*NP(pt1)-1)
    pt3 = pt1
    if (rang() < ipfac) pt3 = pt2
    do i = ND+1, NP(pt3)-1     !> continuous spectra 
        if(found) exit 
        
        c1 = dist%LDAT(K(pt3)+2+2*NP(pt3)+i-1)
        c3 = dist%LDAT(K(pt3)+2+2*NP(pt3)+i)
        if (rn > c1 .and. rn < c3) then 
            p1 = dist%LDAT(K(pt3)+2+NP(pt3)+i-1)
            p2 = dist%LDAT(K(pt3)+2+NP(pt3)+i)
            E1 = dist%LDAT(K(pt3)+2+i-1)
            E2 = dist%LDAT(K(pt3)+2+i)
            
            if (INTT == 1 .or. p1==p2) then  !> histogram
                E_ = E1 + (rn - c1)/p1
                R  = dist%LDAT(K(pt3)+2+3*NP(pt3)+i-1)
                A  = dist%LDAT(K(pt3)+2+4*NP(pt3)+i-1)
                
                
                !if (A == 0) then 
                !    print *, A, dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1), pt1
                !    print *, dist%LDAT(K(pt1)+2+4*NP(pt1) : K(pt1)+2+5*NP(pt1)-1)
                !endif 
            elseif(INTT == 2) then      !> linear-linear
                temp = (p2-p1)/(E2-E1)
                E_ = E1 + (sqrt(p1**2 + 2*temp*(rn-c1)) -p1)/temp 
                
                R  = dist%LDAT(K(pt3)+2+3*NP(pt3)+i-1) &
                    + (dist%LDAT(K(pt3)+2+3*NP(pt3)+i)-dist%LDAT(K(pt3)+2+3*NP(pt3)+i-1))&
                    * (E_ - E1)/(E2-E1)
                
                A  = dist%LDAT(K(pt3)+2+4*NP(pt3)+i-1) &
                    + (dist%LDAT(K(pt3)+2+4*NP(pt3)+i)-dist%LDAT(K(pt3)+2+4*NP(pt3)+i-1))&
                    * (E_ - E1)/(E2-E1)
            else 
                print *, 'ERROR :: LAW44, INTT NUMBER IS WRONG AS ', INTT
            endif 
            
            !pt3 = pt1  !> pt3 is for l
            !if (pt1 == pt2) pt1 = pt1-1
            E1 = dist%LDAT(K(pt1)+2) + ipfac*(dist%LDAT(K(pt1+1)+2)-dist%LDAT(K(pt1)+2))
            E2 = dist%LDAT(K(pt1)+2+NP(pt1)-1) + ipfac*(dist%LDAT(K(pt1+1)+2+NP(pt1+1)-1) - dist%LDAT(K(pt1)+2+NP(pt1)-1))
            !E2 = dist%LDAT(K(pt1)+2+i-1) + ipfac*(dist%LDAT(K(pt1+1)+2+i-1) - dist%LDAT(K(pt1)+2+i-1))

            Eout = E1 + (E_-dist%LDAT(K(pt3)+2))*(E2-E1)/(dist%LDAT(K(pt3)+2+NP(pt3)-1) - dist%LDAT(K(pt3)+2))
            !Eout = E1 + (E_-dist%LDAT(K(pt3)+2))*(E2-E1)/(dist%LDAT(K(pt3)+2+i-1) - dist%LDAT(K(pt3)+2))
            
            found = .true. 
        endif 
    enddo 
    
    !if(found /= .true.) then 
    !    print *, "ERROR :: not found"
    !    stop
    !endif 
    !> Sample CM mu 
    rn = rang() 
    if (rang() > R) then 
        T = (2.*rn -1.)*sinh(A)
        mu = log(T+sqrt(T**2 +1.))/A
    else 
        mu = log(rn*exp(A)+(1.-rn)*exp(-A))/A
    endif 
    
    !if (mu /= mu) then 
    !    print *, mu, T, A, R 
    !    print *, INTT, pt1
    !    print *, dist%LDAT(K(pt1)+2+4*NP(pt1) : K(pt1)+2+5*NP(pt1)-1)        
    !endif 
    
    erg = Eout
    
end subroutine 


! ========================================================= !
!  >>>>>  law 61 (from endf law 1) 
!        -- correlated tab energy-angle dist.
! ========================================================= !
subroutine LAW61 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    integer :: i, j
    logical :: found
    
    integer :: NR
    integer :: NE 
    integer, allocatable :: NP(:), K(:) 
    real(8), allocatable, dimension(:) :: NBT, INTP
    real(8), allocatable, dimension(:) :: E 
    integer, allocatable, dimension(:) :: Loc
    !real(8), allocatable, dimension(:) :: EOUT, PDF, CDF 
    
    integer :: pt1, pt2, pt3 
    real(8) :: rn
    real(8) :: ipfac
    
    integer :: INTT_, INTT, ND
    real(8) :: c1, c2, c3, c4   ! c(i,k), c(i+1,k), c(i,k+1), c(i+1,k+1) 
                                ! respectively in page 2-43 ENDF Manual Vol.1
    real(8) :: p1, p2            ! p(l,k), p(l,k+1)
    real(8) :: E1, E2            ! E(l,k), E(l,k+1) 
    real(8) :: Eout, E_, temp
    
    !> angle distribution parameters 
    integer :: L, LC, JJ, NP_ang
    real(8) :: CSOUT1, CSOUT2, PDF1, PDF2, CDF1, CDF2
        
    
    NR = dist%LDAT(1)
    if (NR /= 0) then 
        print *, "WARNING :: NBT exists (law 61)"
        allocate(NBT(NR)); allocate(INTP(NR)) 
        NBT (1:NR) = dist%LDAT(2   :NR+1  ) 
        INTP(1:NR) = dist%LDAT(NR+2:2*NR+1) 
    endif
    
    NE = dist%LDAT(2+2*NR)
    allocate(E(NE)); E(1:NE) = dist%LDAT(3+2*NR   :3+2*NR+NE-1  ) 
    allocate(loc(NE)); Loc(1:NE) = dist%LDAT(3+2*NR+NE : 3+2*NR+2*NE-1)
    
    !> Determine the incident Energy interval
    pt1 = 1; pt2 = NE
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= E(pt3) ) then 
            pt1 = pt3 
        else 
            pt2 = pt3
        endif 
    enddo 
    ipfac = max(0.d0, min(1.d0,(erg-E(pt1))/(E(pt1+1)-E(pt1))))
    
    !> Sample the outgoing Energy     
    allocate(K(1:NE)); allocate(NP(1:NE))
    !K(1)  = 3+2*NR+2*NE
    !NP(1) = dist%LDAT(K(1)+1)
    !print *, 1, K(1), NP(1)
    !do i = 2, NE 
    !    K(i)  = K(i-1)+2+4*NP(i-1)
    !    NP(i) = dist%LDAT(K(i)+1)
    !    INTT_ = dist%LDAT(K(i))
    !    print *, i, K(i), NP(i), INTT_
    !    print *, dist%LDAT(K(i)-1:K(i)+1) 
    !enddo 
    
    do i = 1, NE
        K(i)    = loc(i) - dist%IDAT+1
        !INTT_    = dist%LDAT(loc(i) - dist%IDAT + 1)
        NP(i)    = dist%LDAT(K(i)+1)
    enddo 
    
    INTT_ = dist%LDAT(K(pt1))
    INTT  = mod(INTT_, 10)
    ND      = (INTT_ - INTT)/10
        
    rn = rang(); found = .false.
    do i = 1, ND              !> discrete spectra 
        if(found) exit 
        c1 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i-1)
        c2 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i-1)
        c3 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i)
        c4 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i)
        
        if (rn > c1+ipfac*(c2-c1) .and. rn < c3+ipfac*(c4-c3)) then 
            Eout  = dist%LDAT(K(pt1)+2+i-1) + ipfac*(dist%LDAT(K(pt2)+2+i-1) - dist%LDAT(K(pt1)+2+i-1))
            found = .true. 
        endif 
        
    enddo
    
    pt3 = pt1
    if (rang() < ipfac) pt3 = pt2
    do i = ND+1, NP(pt3)-1     !> continuous spectra 
        if(found) exit 
        
        c1 = dist%LDAT(K(pt3)+2+2*NP(pt3)+i-1)
        c3 = dist%LDAT(K(pt3)+2+2*NP(pt3)+i)
        
        if (rn > c1 .and. rn < c3) then 
            p1 = dist%LDAT(K(pt3)+2+NP(pt3)+i-1)
            p2 = dist%LDAT(K(pt3)+2+NP(pt3)+i)
            E1 = dist%LDAT(K(pt3)+2+i-1)
            E2 = dist%LDAT(K(pt3)+2+i)
            
            if ( INTT == 1 .or. p1 == p2 ) then  !> histogram
                E_ = E1 + (rn - c1)/p1
            elseif(INTT == 2) then      !> linear-linear
                temp = (p2-p1)/(E2-E1)
                E_ = E1 + (sqrt(p1**2 + 2*temp*(rn-c1)) -p1)/temp 
            else 
                print *, 'ERROR :: LAW61, INTT NUMBER IS WRONG AS ', INTT
            endif 
            
            !pt3 = pt1  !> pt3 is for l
            if (pt1 == pt2) pt1 = pt1-1
            E1 = dist%LDAT(K(pt1)+2) + ipfac*(dist%LDAT(K(pt1+1)+2)-dist%LDAT(K(pt1)+2))
            E2 = dist%LDAT(K(pt1)+2+NP(pt1)-1) + ipfac*(dist%LDAT(K(pt1+1)+2+NP(pt1+1)-1) - dist%LDAT(K(pt1)+2+NP(pt1)-1))
            
            Eout = E1 + (E_-dist%LDAT(K(pt3)+2))*(E2-E1)/(dist%LDAT(K(pt3)+2+NP(pt3)-1) - dist%LDAT(K(pt3)+2))
            j = i
            found = .true. 
        endif 
    enddo 
    !print *, 'check mid'
    !if (found == .false.) then 
    !    print *, "ERROR :: PROBABILITY NOT FOUND"
    !    stop
    !endif
    !> Sample CM mu ==========================================
    LC = dist%LDAT(K(pt3)+2+3*NP(pt3)+j-1)
    !L  = ace(iso)%JXS(11)+abs(LC)-1 
    !L  = LC - (5+2*dist%NR+2*dist%NE)-1 - ace(iso)%JXS(11)
    L = (K(pt3)+2+4*NP(pt3)) +(LC - dist%LDAT(K(pt3)+2+3*NP(pt3))) - 1
    
    !print *, 'check out'
    !if(L <= 0) then 
    !    print *, '                       L', L 
    !    print *, '                     LCi', LC
    !    print *, 'Outgoing energy grid (j)', j, '/', NP(pt1), NP(pt2)
    !    print *, '          starting point', (K(pt1)+2+4*NP(pt1))
    !    print *, '                     LC1', dist%LDAT(K(pt1)+2+3*NP(pt1))
    !    print *, '                  K(pt1)', K(pt1)
    !    print *, '                 NP(pt1)', NP(pt1)
    !    do i = 1, size(dist%LDAT(:))
    !        write(wt_coord,*) dist%LDAT(i)
    !    enddo 
    !    stop 
    !endif 
    
    
    
    deallocate(NP)
    if (LC == 0) then
        !> isotropic scattering (CM system)
        mu = 2.0d0*rang()-1.0d0
    else  !> Tabular probability angular distribution
        rn = rang() 
        NP_ang = dist%LDAT(L+2)
        
        pt1 = 1
        pt2 = NP_ang
        do 
            if (pt2-pt1 ==1) exit 
            pt3 = (pt1+pt2)/2 
            if(rn >= dist%LDAT(L+2+2*NP_ang+pt3-1)) then 
                pt1 = pt3
            else 
                pt2 = pt3
            endif
        enddo
        
        JJ         = dist%LDAT(L+1)
        CSOUT1    = dist%LDAT(L+2+pt1-1)
        CSOUT2    = dist%LDAT(L+2+pt1)
        PDF1    = dist%LDAT(L+2+NP_ang+pt1-1)
        PDF2    = dist%LDAT(L+2+NP_ang+pt1)
        CDF1    = dist%LDAT(L+2+2*NP_ang+pt1-1)
        CDF2    = dist%LDAT(L+2+2*NP_ang+pt1)
        
        
        
        temp = (PDF2-PDF1)/(CSOUT2-CSOUT1)
        if (JJ == 1 .or. temp == 0) then !> Histogram 
            mu = CSOUT1 + (rn-CDF1)/PDF1
        elseif (JJ == 2) then !> Linear-Linear 
            mu = CSOUT1+(sqrt(max(0.0d0,PDF1**2 + 2.0d0*temp*(rn-CDF1)))-PDF1)/temp
        else 
            print *, "ERROR LAW61 :: UNKNOWN INTERPOLATION TYPE IN SCATTERING ANGLE DISTRIBUTIOIN"
        endif    
    endif
    
    !if (abs(mu) > 1. ) then 
    !    print *, 'WARNING :: abs cosine larger than 1', mu 
    !    print *, JJ, temp
    !    print *, CSOUT1
    !    if (JJ == 1 .or. temp == 0) then !> Histogram 
    !        print *, (rn-CDF1)/PDF1
    !    elseif (JJ == 2) then !> Linear-Linear 
    !        print *, (sqrt(max(0.0d0,PDF1**2 + 2.0d0*temp*(rn-CDF1)))-PDF1)/temp
    !    endif 
    !    !print *, (rn-CDF1)/PDF1
    !    !print *, (sqrt(max(0.0d0,PDF1**2 + 2.0d0*temp*(rn-CDF1)))-PDF1)/temp
    !    !mu = sign(1.0d0, mu)
    !    
    !    stop 
    !endif
    !> Check CM or LAB 
    !if (ace(iso)%TY(iMT) < 0) then
    !    Eout = Eout + (erg + 2*mu*(ace(iso)%atn+1)*sqrt(erg * Eout))/(ace(iso)%atn+1)**2
    !endif
    erg = Eout
    
end subroutine 


! ========================================================= !
!  >>>>>  law 44_ANG (from endf law 1) 
!        -- kalbach-87 correlated formalism.
!        -- Only for angle distribution.
! ========================================================= !
subroutine LAW44_ANG (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(in) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    integer :: i, j
    logical :: found
    
    integer :: NR
    integer :: NE 
    integer, allocatable :: NP(:), K(:) 
    real(8), allocatable, dimension(:) :: NBT, INTP
    real(8), allocatable, dimension(:) :: E
    
    integer :: pt1, pt2, pt3 
    real(8) :: rn
    real(8) :: ipfac
    
    integer :: INTT_, INTT, ND
    real(8) :: c1, c2, c3, c4   ! c(i,k), c(i+1,k), c(i,k+1), c(i+1,k+1) 
                                ! respectively in page 2-43 ENDF Manual Vol.1
    real(8) :: p1, p2            ! p(l,k), p(l,k+1)
    real(8) :: E1, E2            ! E(l,k), E(l,k+1) 
    real(8) :: A, R, T
    real(8) :: E_, temp 
    
    NR = dist%LDAT(1)
    if (NR /= 0) then 
        print *, "WARNING :: NBT exists (law 44 ang)"
        allocate(NBT(NR)); allocate(INTP(NR)) 
        NBT (1:NR) = dist%LDAT(2   :NR+1  ) 
        INTP(1:NR) = dist%LDAT(NR+2:2*NR+1) 
    endif
    
    
    NE = dist%LDAT(2+2*NR)
    allocate(E(NE)); E(1:NE) = dist%LDAT(3+2*NR   :2+2*NR+NE  ) 
    !allocate(I(NE)); I(1:NE) = dist%LDAT(3+2*NR+NE:2+2*NR+2*NE) 
    !> Determine the incident Energy interval
    pt1 = 1; pt2 = NE
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= E(pt3) ) then 
            pt1 = pt3 
        else 
            pt2 = pt3
        endif 
    enddo 
    ipfac = max(0.d0, min(1.d0,(erg-E(pt1))/(E(pt1+1)-E(pt1))))
    
    !> Sample the outgoing Energy
    allocate(K(1:NE)); allocate(NP(1:NE))
    K(1)  = 3+2*NR+2*NE
    NP(1) = dist%LDAT(K(1)+1)
    do i = 2, NE 
        K(i)  = K(i-1)+2+5*NP(i-1)
        NP(i) = dist%LDAT(K(i)+1)
    enddo 
        
    INTT_ = dist%LDAT(K(pt1))
    INTT  = mod(INTT_, 10)
    ND      = (INTT_ - INTT)/10
    
    rn = rang(); found = .false.
    do i = 1, ND              !> discrete spectra 
        if(found) exit 
        c1 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i-1)
        c2 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i-1)
        c3 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i)
        c4 = dist%LDAT(K(pt2)+2+2*NP(pt2)+i)
        
        if (rn > c1+ipfac*(c2-c1) .and. rn < c3+ipfac*(c4-c3)) then 
            R = dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1) + &
                ipfac*(dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1) - dist%LDAT(K(pt2)+2+3*NP(pt2)+i-1))
            A = dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1) + &
                ipfac*(dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1) - dist%LDAT(K(pt2)+2+4*NP(pt2)+i-1))
            found = .true. 
        endif 
    enddo 
    
    !print *, dist%LDAT(K(pt1)+2+3*NP(pt1)-1)
    if (rang() < ipfac) pt1 = pt2
    do i = ND+1, NP(pt1)-1     !> continuous spectra 
        if(found) exit 
        
        c1 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i-1)
        c3 = dist%LDAT(K(pt1)+2+2*NP(pt1)+i)
        if (rn > c1 .and. rn < c3) then 
            p1 = dist%LDAT(K(pt1)+2+NP(pt1)+i-1)
            p2 = dist%LDAT(K(pt1)+2+NP(pt1)+i)
            E1 = dist%LDAT(K(pt1)+2+i-1)
            E2 = dist%LDAT(K(pt1)+2+i)
            
            if (INTT == 1 .or. p1==p2) then  !> histogram
                R  = dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1)
                A  = dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1)
                
            elseif(INTT == 2) then      !> linear-linear
                temp = (p2-p1)/(E2-E1)
                E_ = E1 + (sqrt(p1**2 + 2*temp*(rn-c1)) -p1)/temp 
                R  = dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1) &
                    + (dist%LDAT(K(pt1)+2+3*NP(pt1)+i)-dist%LDAT(K(pt1)+2+3*NP(pt1)+i-1))&
                    * (E_ - E1)/(E2-E1)
                
                A  = dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1) &
                    + (dist%LDAT(K(pt1)+2+4*NP(pt1)+i)-dist%LDAT(K(pt1)+2+4*NP(pt1)+i-1))&
                    * (E_ - E1)/(E2-E1)
            else 
                print *, 'ERROR :: LAW44-ANG, INTT NUMBER IS WRONG AS ', INTT
            endif 
                        
            found = .true. 
        endif 
    enddo 
    
    !> Sample CM mu 
    if (rang() > R) then 
        T = (2.0d0*rang() -1.0d0)*sinh(A)
        mu = log(T+sqrt(T**2 +1.))/A
    else 
        rn = rang() 
        mu = log(rn*exp(A)+(1.0d0-rn)*exp(-A))/A
    endif 
    
    
end subroutine 



! ========================================================= !
!  >>>>>  law 7 (from endf law 7) 
!        -- simple maxwell fission spectrum.
! ========================================================= !
subroutine LAW7 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    integer :: i, j
    
    integer :: NR
    integer :: NE 
    
    integer :: pt1, pt2, pt3, itpType
    real(8) :: rn1, rn2, rn3, rn4
    real(8) :: ipfac
    real(8), allocatable :: E(:), T_tbl(:)
	real(8), allocatable :: NBT(:), INTP(:) 
    real(8) :: Eout, C, U, T, temp
    
    NR = dist%LDAT(1) 
    if (NR /= 0) then 
        !print *, "WARNING :: NBT exists (law 7)"
        allocate(NBT(NR)); allocate(INTP(NR)) 
        NBT (1:NR) = dist%LDAT(2   :2+NR-1) 
        INTP(1:NR) = dist%LDAT(NR+2:2*NR+1) 
    endif

    NE = dist%LDAT(2+2*NR) 
    
    allocate(E(1:NE)) 
    allocate(T_tbl(1:NE)) 
    
    E(1:NE)     = dist%LDAT(3+2*NR:3+2*NR+NE-1)
    T_tbl(1:NE) = dist%LDAT(3+2*NR+NE:3+2*NR+2*NE-1)
    U            = dist%LDAT(3+2*NR+2*NE)
    
    pt1=1; pt2=NE
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= E(pt3)) then 
            pt1 = pt3 
        else 
            pt2 = pt3 
        endif 
    enddo
    
    itpType = 2
    do i = 1, NR 
        if (pt1 < NBT(i)) then 
            itpType = INT(i)
            exit 
        endif
    enddo
    
    select case( itpType )
    case( 1 )   !> histogram 
        T = T_tbl(pt1)
    case( 2 )   !> linear-linear 
        T = T_tbl(pt1)+(T_tbl(pt2)-T_tbl(pt1))*(erg-E(pt1))/(E(pt2)-E(pt1))
    case( 3 )   !> linear-log
        T = T_tbl(pt1)+(T_tbl(pt2)-T_tbl(pt1))*log(erg/E(pt1))/log(E(pt2)/E(pt1))
    case( 4 )   !> log-linear 
        T = T_tbl(pt1)*(T_tbl(pt2)/T_tbl(pt1))**((erg-E(pt1))/(E(pt2)-E(pt1)))
    case( 5 )   !> log-log 
        T = T_tbl(pt1)*(T_tbl(pt2)/T_tbl(pt1))**(log(erg/E(pt1))/log(E(pt2)/E(pt1)))
    end select
    
    !ipfac = max(0.d0, min(1.d0,(erg-E(pt1))/(E(pt2)-E(pt1))))
    !
    !T    = T_tbl(pt1) + ipfac*(T_tbl(pt2)-T_tbl(pt1)) 
    !temp = (erg - U)/T
    !C    = T**(-1.5d0) * (sqrt(PI)/2.0d0 * erf(sqrt(temp)) - sqrt(temp)*exp(-temp))**(-1.0d0)
    
    temp = 10
    do while (temp > 1)
        rn1 = rang() 
        rn2 = rang() 
        temp = rn1**2 + rn2**2
    enddo 
    rn3 = rang()
    rn4 = rang() 
    
    Eout = -T*(rn1**2 * log(rn3)/(temp) + log(rn4) )
    
    !> Check CM or LAB 
    !if (ace(iso)%TY(iMT)<0) then
    !    Eout = Eout + (erg + 2*mu*(ace(iso)%atn+1)*sqrt(erg * Eout))/(ace(iso)%atn+1)**2
    !endif
    erg = Eout
    
end subroutine


! ========================================================= !
!  >>>>>  law 9 (from endf law 9) 
!        -- evaporation spectrum.
! ========================================================= !
subroutine LAW9 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(in) :: mu
    
    integer :: NR
    integer :: NE 
    
    integer :: pt1, pt2, pt3, i, idx
    real(8) :: rn1, rn2
    real(8) :: ipfac
    real(8), allocatable :: E(:), T_tbl(:)
    real(8), allocatable :: NBT(:), INT(:)
    real(8) :: Eout, T, temp, U
    integer :: itpType=0 
    
    NR = dist%LDAT(1) 
    NE = dist%LDAT(2+2*NR) 
    NR = dist%LDAT(1) 
    if (NR /= 0) then 
        !print *, "WARNING :: NBT exists (law 9)"
        allocate(NBT(1:NR)); NBT(:) = dist%LDAT(2    : 2+NR-1)
        allocate(INT(1:NR)); INT(:) = dist%LDAT(2+NR : 2+2*NR-1)
    endif
    allocate(E(1:NE)) 
    allocate(T_tbl(1:NE)) 
    E(1:NE)     = dist%LDAT(3+2*NR:3+2*NR+NE-1)
    T_tbl(1:NE) = dist%LDAT(3+2*NR+NE:3+2*NR+2*NE-1)
    U            = dist%LDAT(3+2*NR+2*NE)

    
    pt1=1; pt2=NE
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= E(pt3)) then 
            pt1 = pt3 
        else 
            pt2 = pt3 
        endif 
    enddo
    
    itpType = 2
    do i = 1, NR 
        if (pt1 < NBT(i)) then 
            itpType = INT(i)
            exit 
        endif
    enddo
    
    select case( itpType )
    case( 1 )   !> histogram 
        T = T_tbl(pt1)
    case( 2 )   !> linear-linear 
        T = T_tbl(pt1)+(T_tbl(pt2)-T_tbl(pt1))*(erg-E(pt1))/(E(pt2)-E(pt1))
    case( 3 )   !> linear-log
        T = T_tbl(pt1)+(T_tbl(pt2)-T_tbl(pt1))*log(erg/E(pt1))/log(E(pt2)/E(pt1))
    case( 4 )   !> log-linear 
        T = T_tbl(pt1)*(T_tbl(pt2)/T_tbl(pt1))**((erg-E(pt1))/(E(pt2)-E(pt1)))
    case( 5 )   !> log-log 
        T = T_tbl(pt1)*(T_tbl(pt2)/T_tbl(pt1))**(log(erg/E(pt1))/log(E(pt2)/E(pt1)))
    end select    
    !ipfac = max(0.d0, min(1.d0,(erg-E(pt1))/(E(pt2)-E(pt1))))
    !T     = T_tbl(pt1) + ipfac*(T_tbl(pt2)-T_tbl(pt1)) 
    
    !print *, 'T', T, pt1,pt2, erg
    !do i = 1, NE 
    !    print *, i, E(i), T_tbl(i)
    !enddo 
    !stop
    
    do
        rn1 = rang() 
        rn2 = rang() 
        Eout = -T*log(rn1*rn2) 
        
        if (Eout <= (erg - U)) exit
    enddo 
    
    erg = Eout
    
    
    if (allocated(T_tbl)) deallocate(T_tbl)
    if (allocated(E))     deallocate(E)
    if (allocated(INT))   deallocate(INT)
    if (allocated(NBT))   deallocate(NBT)
    
end subroutine




! ========================================================= !
!  >>>>>  law 11 (from endf law 11) 
!        -- energy dependent watt spectrum.
! ========================================================= !
subroutine LAW11 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    
    integer :: NRa, NRb 
    integer :: NEa, NEb 
    integer :: L
    
    integer :: pt1, pt2, pt3
    real(8) :: rn1, rn2
    real(8) :: ipfac
    real(8), allocatable :: Ea(:), Eb(:), a_tbl(:), b_tbl(:)
    real(8) :: Eout, a, b, temp, U, g
    
    
    NRa = dist%LDAT(1)
    if(NRa /= 0) print *, "WARNING :: NBT exists (law 11)"
    NEa = dist%LDAT(2+2*NRa) 
    allocate(Ea(1:NEa)) 
    allocate(a_tbl(1:NEa)) 
    Ea(1:NEa)    = dist%LDAT(3+2*NRa     : 3+2*NRa+NEa-1)
    a_tbl(1:NEa) = dist%LDAT(3+2*NRa+NEa : 3+2*NRa+2*NEa-1)
    L = 3+2*(NRa + NEa)
    
    NRb = dist%LDAT(L) 
    NEb = dist%LDAT(L+1+2*NRb)
    if(NRb /= 0) print *, "WARNING :: NBT exists (law 11)"
    allocate(Eb(1:NEb)) 
    allocate(b_tbl(1:NEb))
    Eb(1:NEb)    = dist%LDAT(L+2+2*NRb        : L+2+2*NRb+NEb-1) 
    b_tbl(1:NEb) = dist%LDAT(L+2+2*NRb+NEb    : L+2+2*NRb+2*NEb-1) 
    U            = dist%LDAT(L+2+2*NRb+2*NEb) 
    
    
    pt1=1; pt2=NEa
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= Ea(pt3)) then 
            pt1 = pt3 
        else 
            pt2 = pt3 
        endif 
    enddo
    ipfac = max(0.d0, min(1.d0,(erg-Ea(pt1))/(Ea(pt2)-Ea(pt1))))
    a     = a_tbl(pt1) + ipfac*(a_tbl(pt2)-a_tbl(pt1)) 
    
    pt1=1; pt2=NEb
    do 
        if(pt2-pt1 == 1) exit 
        pt3 = (pt2+pt1)/2 
        if (erg >= Eb(pt3)) then 
            pt1 = pt3 
        else 
            pt2 = pt3 
        endif 
    enddo
    ipfac = max(0.d0, min(1.d0,(erg-Eb(pt1))/(Eb(pt2)-Eb(pt1))))
    b     = b_tbl(pt1) + ipfac*(b_tbl(pt2)-b_tbl(pt1)) 
    
    
    g = sqrt((1+a*b/8.0d0)**2 -1) + (1+a*b/8.0d0)
    
    loop : do
        rn1 = rang() 
        rn2 = rang() 
        Eout = -a*g*log(rn1)
        
        temp = ((1-g)*(1-log(rn1))-log(rn2))**2 
        if (temp <= b*Eout) exit loop
    enddo loop
    
    !> Check CM or LAB 
    !if (ace(iso)%TY(iMT)<0) then
    !    Eout = Eout + (erg + 2*mu*(ace(iso)%atn+1)*sqrt(erg * Eout))/(ace(iso)%atn+1)**2
    !endif
    erg = Eout
    
end subroutine

    
    
    
! ========================================================= !
!  >>>>>  law 66 (from endf law 6) 
!        -- n-body phase space distribution.
! ========================================================= !
subroutine LAW66 (erg, iso, mu, dist, iMT) 
    !type(particle), intent(inout) :: p 
    real(8), intent(inout) :: erg
    type(EnergyDistDataForm), intent(in) :: dist
    integer, intent(in) :: iso 
    integer, intent(in) :: iMT
    real(8), intent(inout) :: mu
    
    integer :: NPSX 
    real(8) :: Ap, A
    real(8) :: rn1, rn2, rn3, rn4, rn5, rn6, rn7, rn8, rn9
    real(8) :: pp, x, y, T 
    real(8) :: Eout, Emax
    
    
    NPSX = dist%LDAT(1) 
    Ap     = dist%LDAT(2)
    A     = ace(iso)%atn 
    
    Emax = (Ap-1)/Ap *(A*erg/(A+1) + ace(iso)%Q(iMT))
    
    do
        rn1 = rang() 
        rn2 = rang() 
        if (rn1**2 + rn2**2 <= 1) exit
    enddo
    
    do
        rn3 = rang() 
        rn4 = rang() 
        if (rn3**2 + rn4**2 <= 1) exit
    enddo
    
    select case (NPSX) 
    case(3) 
        rn5 = rang() 
        pp  = rn5
    case(4) 
        rn5 = rang() 
        rn6 = rang() 
        pp  = rn5*rn6 
    case(5) 
        rn5 = rang() 
        rn6 = rang() 
        rn7 = rang() 
        rn8 = rang() 
        pp  = rn5*rn6*rn7*rn8 
    case default 
        print *, "WARNING :: N-BODY PHASE SPACE DISTRIBUTION - BODY NUMBER : ", NPSX
    end select 
    
    rn9 = rang()
    x = -rn1*log(rn1**2 + rn2**2)/(rn1**2 + rn2**2) - log(rn9)
    y = -rn3*log(rn3**2 + rn4**2)/(rn3**2 + rn4**2) - log(pp)
    T = x/(x+y)
    
    Eout = T*Emax
    mu   = 2*rn2 - 1
    
    !> Check CM or LAB 
    !if (ace(iso)%TY(iMT)<0) then
    !    Eout = Eout + (erg + 2*mu*(ace(iso)%atn+1)*sqrt(erg * Eout))/(ace(iso)%atn+1)**2
    !endif
    erg = Eout
    
    
end subroutine

! ========================================================= !
!    ENDF Interpolation scheme 
! ========================================================= !
subroutine ENDF_interpolation (scheme)
    integer, intent(in) :: scheme 
    
    select case (scheme) 
    case (1) 
        
    case (2) 
        
    case (3) 
        
    case (4) 
        
    case (5)
        
    case (6) 
        print *, "ENDF INTERPOLTAION :: Gamow charged-particle penetrability law is not prepared"
    case default 
        print *, "ERROR :: UNKNOWN ENDF INTERPOLATION SCHEME -", scheme 
        stop 
    end select
    
    
end subroutine 

end module 

