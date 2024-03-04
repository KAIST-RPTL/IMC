module surface_header

    use constants!, only : INFINITY, TOOLONG, EPSILON
    use omp_lib
    
    implicit none 
    
    
        
    type Surface
        character(20) :: surf_id            !> surface id
        integer       :: surf_type            
        integer       :: bc                 !> 0 : normal
                                            !> 1 : vacuum
                                            !> 2 : reflective
        
        real(8), dimension(10) :: parmtrs    !> p   : 
                                            !> sqr : 
                                            !> cyl : 
                                            !> sph : 

        integer, allocatable :: &
             neighbor_pos(:), &               ! List of cells on positive side
             neighbor_neg(:)                  ! List of cells on negative side
            
    contains 
        !procedure :: find_idx => find_surf_idx
    
    end type
    type(Surface), allocatable, target :: Surfaces(:), Surfaces_temp(:)
    type(surface) :: surf_obj
    

    
    
    contains
    
    subroutine read_surf (surfobj, args, nargs)
        type(Surface) :: surfobj
		character(*) :: args(:)
		integer :: nargs
        character(30):: surf_type
        integer :: i
		
		read(args(2), *) surfobj%surf_id
		read(args(3), *) surf_type
        surfobj%surf_type = surf_type_converter(surf_type)
		
		do i = 1, nargs-3
			read(args(i+3), *) surfobj%parmtrs(i)
        enddo 
        surfobj%bc        = 0 ! by default 
		
		
    end subroutine 
    
	subroutine EMSG_surf(surf_type, nargs,line_idx)
		integer :: surf_type, nargs,line_idx
		integer :: n 
		
		select case (surf_type)
		case(1:3) 
			n = 4
		case(4:9)
			n = 6
		case(10) 
			n = 7 
		case(11:12) 
			n = 6
        case(13)
            n = 7
        case(14)
            n = 9
		case default 
			n = 0 
		end select
		
		if (nargs /= n) then 
			print '(a,i4,a)', "geom.inp (Line ",line_idx,") Wrong number of surface parameters "
			stop 
		endif 
		
	end subroutine 
    
    
    
    function surf_type_converter(surf_type) result(output) 
        character(*), intent(in) :: surf_type
        integer :: output
        output = 0 
        
        if (surf_type.eq.'px')   output = 1
        if (surf_type.eq.'py')   output = 2
        if (surf_type.eq.'pz')   output = 3
        if (surf_type.eq.'sqcx') output = 4
        if (surf_type.eq.'sqcy') output = 5
        if (surf_type.eq.'sqcz') output = 6
        if (surf_type.eq.'cylx') output = 7
        if (surf_type.eq.'cyly') output = 8
        if (surf_type.eq.'cylz') output = 9
        if (surf_type.eq.'sph')  output = 10
		if (surf_type.eq.'hexxc') output = 11
        if (surf_type.eq.'hexyc') output = 12
        if (surf_type.eq.'rect')  output = 13
        if (surf_type.eq.'cuboid') output = 14
        
        if (output == 0) then 
            print *, '******************************************'
            print *, 'ERROR: WRONG SURFACE TYPE - ', surf_type 
            print *, '******************************************'
            stop  
        endif 
        
        
    end function 
    
    function find_surf_idx(this, surf_id) result (idx)
        type(surface) :: this(:)
        integer ::  i, idx
        character(*):: surf_id
        
        do i = 1, size(this) 
            if (this(i)%surf_id == surf_id) then 
                idx = i 
                goto 1
            endif
        enddo 
        print *, "no such surface id : ", surf_id 
        stop 
1        continue
    end function 


    function surf_neg_or_pos(this, xyz) result(neg) 
        type(surface) :: this
        real(8) :: xyz(3), xyz_tr(3), val
        logical :: neg
        integer :: type
        real(8) :: d_temp(6) 
		real(8) :: r, tmp
		
        neg = .false.
        
        type = this % surf_type
        select case (type)
        case(1) !> px
            if (xyz(1) < this%parmtrs(1)) neg = .true.
        case(2) !> py
            if (xyz(2) < this%parmtrs(1)) neg = .true.
        case(3) !> pz
            if (xyz(3) < this%parmtrs(1)) neg = .true.
            
        case(6) !> sqcz
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            xyz_tr(3) = xyz(3)
            if ((abs(xyz_tr(1))) < this%parmtrs(3)&
                .and.abs(xyz_tr(2)) < this%parmtrs(3))  &
                neg = .true.
            
        case(9) !> cylz
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            xyz_tr(3) = xyz(3)
            
            val = sqrt((xyz_tr(1))**2 + (xyz_tr(2))**2 ) 
            if (val < this%parmtrs(3)) neg = .true.
            
        case(10) !> sph
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            xyz_tr(3) = xyz(3) - this%parmtrs(3)
            
            val = sqrt((xyz_tr(1))**2 + (xyz_tr(2))**2 +(xyz_tr(3))**2 ) 
            if (val < this%parmtrs(4)) neg = .true.
            
		!case(11) !> hexz
		!	! index  1 
		!	!       ---
		!	!	 4 /   \ 2
		!	!    5 \   / 3
		!	!       ---
		!	!		 6
		!	
        !    xyz_tr(1) = xyz(1) - this%parmtrs(1)
        !    xyz_tr(2) = xyz(2) - this%parmtrs(2)
		!	r = this%parmtrs(3)
		!	
		!	! idx = 1
		!	if (xyz_tr(2) > r) return
		!	
		!	! idx = 2
		!	val = -2.0*xyz_tr(1) + 2.0*r 
		!	if (xyz_tr(2) > val) return
		!	
		!	! idx = 3
		!	val = 2.0*xyz_tr(1) - 2.0*r 
		!	if (xyz_tr(2) < val) return
		!	
		!	! idx = 4
		!	val = 2.0*xyz_tr(1) + 2.0*r 
		!	if (xyz_tr(2) < val) return
		!	
		!	! idx = 5
		!	val = -2.0*xyz_tr(1) - 2.0*r 
		!	if (xyz_tr(2) > val) return
		!	
		!	! idx = 6
		!	if (xyz_tr(2) < -r) return
		!	
		!	neg = .true.
		!	
        !========== TEST FOR HEXAGONAL ===========\
        case(11) !> hexxc
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            r = this%parmtrs(3)
            tmp = sqrt(3.0)*0.5d0
            ! X-type hex. prism
            ! PLANE INDICATOR
            !
            !1  /\ 3
            !2 |  | 2
            !3  \/ 1
            neg = .false.
            ! SLAB 1
            val = xyz_tr(1)
            if(val>=r .or. val<=-r) return
            ! SLAB 2
            val = xyz_tr(2)*tmp-xyz_tr(1)*0.5d0
            if(val>=r .or. val<=-r) return
            ! SLAB 3
            val = xyz_tr(2)*tmp+xyz_tr(1)*0.5d0
            if(val>=r .or. val<=-r) return
            neg = .true.
        case(12) !> hexyc
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            r = this%parmtrs(3)
            tmp = sqrt(3.0)*0.5d0
            ! Y-type hex. prism
            ! PLANE INDICATOR
            !    1
            !    __
            !2  /  \ 3
            !3  \__/ 2
            !  
            !    1
            ! SLAB #1
            neg = .false.
            val = xyz_tr(2)
            if(val>=r .or. val<=-r) return
            ! SLAB #2
            val = xyz_tr(1)*tmp-xyz_tr(2)*0.5d0
            if(val>=r .or. val<=-r) return
            ! SLAB #3
            val = xyz_tr(1)*tmp+xyz_tr(2)*0.5d0
            if(val>=r .or. val<=-r) return
            neg = .true.
        case(13) !> Rect
            if( xyz(1)>this%parmtrs(1) .and. &
                xyz(1)<this%parmtrs(2) .and. &
                xyz(2)>this%parmtrs(3) .and. &
                xyz(2)<this%parmtrs(4) ) neg = .true.
        case(14) !> cuboid
            if( xyz(1)>this%parmtrs(1) .and. &
                xyz(1)<this%parmtrs(2) .and. &
                xyz(2)>this%parmtrs(3) .and. &
                xyz(2)<this%parmtrs(4) .and. &
                xyz(3)>this%parmtrs(5) .and. &
                xyz(3)<this%parmtrs(6)) neg = .true.
        end select
        
        
    end function 

    !> distance to surface boundary
    function surf_gp(A,B,C,D,xyz,uvw) result(dist) 
        real(8),intent(in) :: A,B,C,D, xyz(3), uvw(3)
        real(8) :: dist
        
        dist = (D-(A*xyz(1)+B*xyz(2)+C*xyz(3)))/(A*uvw(1)+B*uvw(2)+C*uvw(3))
        
    end function  
    
    function surf_px(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        
        if (surf%surf_type /= px) print *, "ERROR : WRONG SURFACE" 
        
        if (uvw(1) == 0) then 
            dist = INFINITY
        else 
            dist = (surf%parmtrs(1)-xyz(1))/uvw(1)
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    function surf_py(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        
        if (surf%surf_type /= py) print *, "ERROR : WRONG SURFACE" 
        
        if (uvw(2) == 0) then 
            dist = INFINITY
        else 
            dist = (surf%parmtrs(1)-xyz(2))/uvw(2)
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    function surf_pz(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        
        if (surf%surf_type /= pz) print *, "ERROR : WRONG SURFACE" 
        
        !if (uvw(3) == 0) then 
        !    dist = INFINITY
        !else 
            dist = (surf%parmtrs(1)-xyz(3))/uvw(3)
        !endif
        
        if (dist<0) dist = INFINITY
        
    end function
    
    
    
    function surf_cylz(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        real(8) :: a, k, c, xyz_(3), val 
        integer :: i 
        
        if (surf%surf_type /= cylz) print *, "ERROR : WRONG SURFACE" 
		dist = INFINITY
		
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2)
        
        a = uvw(1)**2 + uvw(2)**2
        c = xyz_(1)**2 + xyz_(2)**2 - surf%parmtrs(3)**2
        k = xyz_(1)*uvw(1) + xyz_(2)*uvw(2)
        val = k**2 - a*c
        
        if ((a == 0).or.(val < 0)) then 
            dist = INFINITY
        elseif (c < 0) then
            dist = (-k + sqrt(val))/a
            
        elseif (c > 0) then
            dist = (-k - sqrt(val))/a
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    
    function surf_sqcz(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: d(4), xyz_(3), temp, r
        real(8) :: dist 
        integer :: i
        
        if (surf%surf_type /= sqcz) print *, "ERROR : WRONG SURFACE" 
        
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2) 
        xyz_(3)   = xyz(3)  

        r = surf%parmtrs(3) 
        d(1) = (-r-xyz_(2))/uvw(2)
        d(2) = (-r-xyz_(1))/uvw(1)
        d(3) = ( r-xyz_(2))/uvw(2)
        d(4) = ( r-xyz_(1))/uvw(1)
        
        temp = xyz_(1)+d(1)*uvw(1)
        if ((temp < -r).or.(temp > r)) d(1) = INFINITY
        
        temp = xyz_(2)+d(2)*uvw(2)
        if ((temp < -r).or.(temp > r)) d(2) = INFINITY
        
        temp = xyz_(1)+d(3)*uvw(1)
        if ((temp < -r).or.(temp > r)) d(3) = INFINITY
        
        temp = xyz_(2)+d(4)*uvw(2)
        if ((temp < -r).or.(temp > r)) d(4) = INFINITY
        
        do i = 1, 4 
            if (d(i) < 0) d(i) = INFINITY 
        enddo 
        
        
        dist = minval(d(:))
    end function
    
    function surf_sph(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(3), k, c
        real(8) :: dist, temp
        integer :: i 
        
        if (surf%surf_type /= sph) print *, "ERROR : WRONG SURFACE" 
        
        k = 0; c = 0 
        do i = 1, 3 
            xyz_(i) = xyz(i) - surf%parmtrs(i) 
            k = k + xyz_(i) * uvw(i)
            c = c + xyz_(i)**2 
        enddo 
        c = c - surf%parmtrs(4)**2
        temp = k**2 - c
        if (temp < 0) then
            dist = INFINITY
        elseif (c < 0) then 
            dist = -k + sqrt(temp)
        elseif (c > 0) then 
            dist = -k - sqrt(temp)
		elseif (c == 0 ) then 
			!print *, 'WTF? particle on the sphere'
			dist = 0
        endif 
        
        if (dist < 0) dist = INFINITY 
        
        
    end function
    
    function surf_hexz(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(3), temp, r, val 
        real(8) :: dist, d_temp(6)
        integer :: i
        real(8) :: a, b, d 
		
        !if (surf%surf_type /= hexz) print *, "ERROR : WRONG SURFACE" 
        
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2) 
        xyz_(3)   = xyz(3)  
		r		  = surf%parmtrs(3)
		
		dist = INFINITY
		
		! idx = 1
		a = 0; b = 1; d = r
		temp = (d-a*xyz_(1)-b*xyz_(2))/(a*uvw(1)+b*uvw(2))
		val  = xyz_(1)+uvw(1)*temp
		if (temp < 0 .or. val < -r/sqrt(3.0) .or. val > r/sqrt(3.0)) temp = INFINITY 
		d_temp(1)=temp
		
		! idx = 2
		a = 2; b = 1; d = 2.0*r
		temp = (d-a*xyz_(1)-b*xyz_(2))/(a*uvw(1)+b*uvw(2))
		val  = xyz_(2)+uvw(2)*temp
		if (temp < 0 .or. val < 0 .or. val > r) temp = INFINITY 
		d_temp(2)=temp
		
		
		! idx = 3
		a = -2; b = 1; d = -2.0*r
		temp = (d-a*xyz_(1)-b*xyz_(2))/(a*uvw(1)+b*uvw(2))
		val  = xyz_(2)+uvw(2)*temp
		if (temp < 0 .or. val < -r .or. val > 0) temp = INFINITY 
		d_temp(3)=temp
		
		
		! idx = 4
		a = -2; b = 1; d = 2.0*r
		temp = (d-a*xyz_(1)-b*xyz_(2))/(a*uvw(1)+b*uvw(2))
		val  = xyz_(2)+uvw(2)*temp
		if (temp < 0 .or. val < 0 .or. val > r) temp = INFINITY 
		d_temp(4)=temp
		
		
		! idx = 5
		a = 2; b = 1; d = -2.0*r
		temp = (d-a*xyz_(1)-b*xyz_(2))/(a*uvw(1)+b*uvw(2))
		val  = xyz_(2)+uvw(2)*temp
		if (temp < 0 .or. val < -r .or. val > 0) temp = INFINITY 
		d_temp(5)=temp
		
		! idx = 6
		a = 0; b = 1; d = -r
		temp = (d-a*xyz_(1)-b*xyz_(2))/(a*uvw(1)+b*uvw(2))
		val  = xyz_(1)+uvw(1)*temp
		if (temp < 0 .or. val < -r/sqrt(3.0) .or. val > r/sqrt(3.0)) temp = INFINITY 
		d_temp(5)=temp
		
		
		
        dist = minval(d_temp(:))
        
    end function	
	
	
	
    
    function surf_hexxc(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(3), ixyz(3), iuvw(3), temp, r, val 
        real(8) :: dist!, d_temp(3)
        integer :: i
        real(8) :: a, b, d(6), temp1, temp2
		
        if (surf%surf_type /= hexxc) print *, "ERROR : WRONG SURFACE" 
        
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2) 
        xyz_(3)   = xyz(3)  
		r		  = surf%parmtrs(3)
        val = sqrt(3.d0)*0.5d0
		
        ixyz(1) = xyz_(1) !X, perp. to surf.2
        ixyz(2) = xyz_(2)*val-xyz_(1)*0.5d0 !Perpendicular to surf. 1
        ixyz(3) = xyz_(2)*val+xyz_(1)*0.5d0 !Perpendicular to surf. 3
        !ixyz(4) = xyz_(3) !Z
        
        iuvw(1) = uvw(1) !U, perp. to surf.2
        iuvw(2) = uvw(2)*val-uvw(1)*0.5d0 !Perp. to surf. 1
        iuvw(3) = uvw(2)*val+uvw(1)*0.5d0 !Perp. to surf. 2
        !iuvw(4) = uvw(3) !W

        d(1) = ( r-ixyz(1))/iuvw(1)
        d(2) = (-r-ixyz(1))/iuvw(1)
        d(3) = ( r-ixyz(2))/iuvw(2)
        d(4) = (-r-ixyz(2))/iuvw(2)
        d(5) = ( r-ixyz(3))/iuvw(3)
        d(6) = (-r-ixyz(3))/iuvw(3)

        
        temp1 = ixyz(2)+d(1)*iuvw(2); temp2 = ixyz(3)+d(1)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(1)<0 .or. iuvw(1)==0) d(1) = INFINITY

        temp1 = ixyz(2)+d(2)*iuvw(2); temp2 = ixyz(3)+d(2)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(2)<0 .or. iuvw(1)==0) d(2) = INFINITY

        temp1 = ixyz(1)+d(3)*iuvw(1); temp2 = ixyz(3)+d(3)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(3)<0 .or. iuvw(2)==0) d(3) = INFINITY
    
        temp1 = ixyz(1)+d(4)*iuvw(1); temp2 = ixyz(3)+d(4)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(4)<0 .or. iuvw(2)==0) d(4) = INFINITY

        temp1 = ixyz(1)+d(5)*iuvw(1); temp2 = ixyz(2)+d(5)*iuvw(2)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(5)<0 .or. iuvw(3)==0) d(5) = INFINITY

        temp1 = ixyz(1)+d(6)*iuvw(1); temp2 = ixyz(2)+d(6)*iuvw(2)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(6)<0 .or. iuvw(3)==0) d(6) = INFINITY

        dist = minval(d(:))
    end function

    function surf_hexyc(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(3), ixyz(3), iuvw(3), temp, r, val 
        real(8) :: dist, d_temp(6)
        integer :: i
        real(8) :: a, b, d (6), temp1, temp2
		
        if (surf%surf_type /= hexyc) print *, "ERROR : WRONG SURFACE" 
        
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2) 
        xyz_(3)   = xyz(3)  
		r		  = surf%parmtrs(3)
        val = sqrt(3.0)*0.5d0
		
        ixyz(1) = xyz_(1)*val-xyz_(2)*0.5d0
        ixyz(2) = xyz_(1)*val+xyz_(2)*0.5d0
        ixyz(3) = xyz_(2)

        iuvw(1) = uvw(1)*val-uvw(2)*0.5d0
        iuvw(2) = uvw(1)*val+uvw(2)*0.5d0
        iuvw(3) = uvw(2)
            
        d(1) = ( r-ixyz(1))/iuvw(1)
        d(2) = (-r-ixyz(1))/iuvw(1)
        d(3) = ( r-ixyz(2))/iuvw(2)
        d(4) = (-r-ixyz(2))/iuvw(2)
        d(5) = ( r-ixyz(3))/iuvw(3)
        d(6) = (-r-ixyz(3))/iuvw(3)

        temp1 = ixyz(2)+d(1)*iuvw(2); temp2 = ixyz(3)+d(1)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(1)<0 .or. iuvw(1)==0) d(1) = INFINITY

        temp1 = ixyz(2)+d(2)*iuvw(2); temp2 = ixyz(3)+d(2)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(2)<0 .or. iuvw(1)==0) d(2) = INFINITY

        temp1 = ixyz(1)+d(3)*iuvw(1); temp2 = ixyz(3)+d(3)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(3)<0 .or. iuvw(2)==0) d(3) = INFINITY

        temp1 = ixyz(1)+d(4)*iuvw(1); temp2 = ixyz(3)+d(4)*iuvw(3)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(4)<0 .or. iuvw(2)==0) d(4) = INFINITY

        temp1 = ixyz(1)+d(5)*iuvw(1); temp2 = ixyz(2)+d(5)*iuvw(2)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(5)<0 .or. iuvw(3)==0) d(5) = INFINITY

        temp1 = ixyz(1)+d(6)*iuvw(1); temp2 = ixyz(2)+d(6)*iuvw(2)
        if((abs(temp1)>r).or.(abs(temp2)>r) .or. d(6)<0 .or. iuvw(3)==0) d(6) = INFINITY

        dist = minval(d(:))
    end function

    function surf_rect(surf, xyz, uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(2), x0, x1, y0, y1
        integer :: i
        real(8) :: dist, d(4), temp

        if(surf % surf_type /= rect) print *, 'ERROR: WRONG SURF.'

        x0 = surf % parmtrs(1); x1 = surf % parmtrs(2)
        y0 = surf % parmtrs(3); y1 = surf % parmtrs(4)

        xyz_(:) = xyz(1:2)

        d(1) = (x0 - xyz_(1))/uvw(1)
        d(2) = (x1 - xyz_(1))/uvw(1)
        d(3) = (y0 - xyz_(2))/uvw(2)
        d(4) = (y1 - xyz_(2))/uvw(2)

        temp = xyz_(2) + d(1) * uvw(2)
        if((temp < -y0) .or. (temp > y1)) d(1) = INFINITY
        temp = xyz_(2) + d(2) * uvw(2)
        if((temp < -y0) .or. (temp > y1)) d(2) = INFINITY
        temp = xyz_(1) + d(3) * uvw(1)
        if((temp < -x0) .or. (temp > x1)) d(3) = INFINITY
        temp = xyz_(1) + d(4) * uvw(1)
        if((temp < -x0) .or. (temp > x1)) d(4) = INFINITY

        do i = 1, 4
            if(d(i) < 0) d(i) = INFINITY
        enddo

        dist = minval(d(:))
    end function

    function surf_cuboid(surf, xyz, uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(3), x0, x1, y0, y1, z0, z1
        integer :: i
        real(8) :: dist, d(6), temp1, temp2

        if(surf % surf_type /= CUBOID) print *, 'ERROR: WRONG SURF.'

        x0 = surf % parmtrs(1); x1 = surf % parmtrs(2)
        y0 = surf % parmtrs(3); y1 = surf % parmtrs(4)
        z0 = surf % parmtrs(5); z1 = surf % parmtrs(6)

        xyz_(:) = xyz(1:3)

        d(1) = (x0 - xyz_(1))/uvw(1)
        d(2) = (x1 - xyz_(1))/uvw(1)
        d(3) = (y0 - xyz_(2))/uvw(2)
        d(4) = (y1 - xyz_(2))/uvw(2)
        d(5) = (z0 - xyz_(3))/uvw(3)
        d(6) = (z1 - xyz_(3))/uvw(3)

        temp1 = xyz_(2) + d(1) * uvw(2)
        temp2 = xyz_(3) + d(1) * uvw(3)
        if((temp1 < -y0) .or. (temp1 > y1) .or. &
           (temp2 < -z0) .or. (temp2 > z1)) d(1) = INFINITY
        temp1 = xyz_(2) + d(2) * uvw(2)
        temp2 = xyz_(3) + d(2) * uvw(3)
        if((temp1 < -y0) .or. (temp1 > y1) .or. &
           (temp2 < -z0) .or. (temp2 > z1)) d(2) = INFINITY
        temp1 = xyz_(1) + d(3) * uvw(1)
        temp2 = xyz_(3) + d(3) * uvw(3)
        if((temp1 < -x0) .or. (temp1 > x1) .or. &
           (temp2 < -z0) .or. (temp2 > z1)) d(3) = INFINITY
        temp1 = xyz_(1) + d(4) * uvw(1)
        temp2 = xyz_(3) + d(4) * uvw(3)
        if((temp1 < -x0) .or. (temp1 > x1) .or. &
           (temp2 < -z0) .or. (temp2 > z1)) d(4) = INFINITY
        temp1 = xyz_(1) + d(5) * uvw(1)
        temp2 = xyz_(2) + d(5) * uvw(2)
        if((temp1 < -x0) .or. (temp1 > x1) .or. &
           (temp2 < -y0) .or. (temp2 > y1)) d(5) = INFINITY
        temp1 = xyz_(1) + d(6) * uvw(1)
        temp2 = xyz_(2) + d(6) * uvw(2)
        if((temp1 < -x0) .or. (temp1 > x1) .or. &
           (temp2 < -y0) .or. (temp2 > y1)) d(6) = INFINITY
        

        do i = 1, 6
            if(d(i) < 0) d(i) = INFINITY
        enddo

        dist = minval(d(:))
    end function


    subroutine surf_select(surf,xyz,uvw, dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist
		
		dist = 0
        if (surf%surf_type == px)   dist = surf_px(surf,xyz,uvw)
        if (surf%surf_type == py)   dist = surf_py(surf,xyz,uvw)
        if (surf%surf_type == pz)   dist = surf_pz(surf,xyz,uvw)
        if (surf%surf_type == cylz) dist = surf_cylz(surf,xyz,uvw)
        if (surf%surf_type == sqcz) dist = surf_sqcz(surf,xyz,uvw)
        if (surf%surf_type == sph)  dist = surf_sph(surf,xyz,uvw)
        if (surf%surf_type == hexxc) dist = surf_hexxc(surf,xyz,uvw)
        if (surf%surf_type == hexyc) dist = surf_hexyc(surf,xyz,uvw)
        if (surf%surf_type == rect)  dist = surf_rect(surf, xyz,uvw)
        if (surf%surf_type == cuboid) dist = surf_cuboid(surf,xyz,uvw)
        !if (dist<TINY_BIT*10) print *, 'STUCK?',xyz(1:2),dist,surf%surf_type
    end subroutine
    
    
end module 



