module multiphysics 
	use tetrahedral
	use strings,  only : parse, readline, removesp, uppercase, insertstr, value
	use evaluate, only : defparam, evalexpr

	implicit none
	
	logical :: do_T=.false., do_TM=.false.
	integer, parameter :: rp = dp
	
	
	real(kind=rp), allocatable, dimension(:) :: T_node
	real(kind=rp) :: T_tet
	integer, allocatable :: T_partlist(:)   ! list of parts for temperature analysis
	integer, allocatable :: TM_partlist(:)  ! list of parts for thermo-mechanical analysis
	integer, parameter :: maxstring = 100
	
	type typemat 
		character(10) :: name
		logical :: coolant=.false.
		real(kind=rp) :: Tf 
		real(kind=rp) :: rho 
		real(kind=rp) :: c_p
		
		character(maxstring)  :: k
		character(maxstring)  :: strain
		character(maxstring)  :: E
		character(maxstring)  :: G
		character(maxstring)  :: hf
	endtype
	type(typemat), allocatable :: mat(:) 
	integer :: n_mat=0, i, n
	character(maxstring) :: args_temp(10)
	integer :: curr_line=0
	integer :: io_input = 0
	integer :: rd_mp = 1
	
	
	! Target Properties
	
	
	
	! Transient Parameter 
	real(8) :: t_tot, dt
	integer :: n_timestep, i_timestep
	
	real(kind=rp), allocatable, dimension(:) :: T_1, T_0
	
	

	contains 
	
	
	
	
	subroutine HeatTransfer
		implicit none
		
		! 1. Heat transfer 
		integer :: i, j, k, ii, jj, i_tet
		real(kind=rp), allocatable, dimension(:) :: Q, Q_std, Rh, Rh_temp, Rq
		real(kind=rp), allocatable, dimension(:,:) :: Kc, Kh, Kc_temp, Kh_temp 
		real(kind=rp), allocatable, dimension(:,:) :: B_T, B
		real(kind=rp) :: area, val, vol
		integer :: idx(3)
		real(kind=rp) :: Tf
		real(kind=rp) :: h=0, k_c
		real(kind=rp) :: T_tet
		logical :: perimeter
		
		! 1. Solve for heat transfer ==============================================	
		! Allocate matrices for heat transfer
		allocate(Q(1:num_tet)) 
		allocate(Q_std(1:num_tet)) 
		allocate(Kc(1:num_node,1:num_node),Kh(1:num_node,1:num_node))
		allocate(Rq(1:num_node), Rh(1:num_node))
		allocate(Kc_temp(1:4,1:4),Kh_temp(1:4,1:4),Rh_temp(1:3)) 
		allocate(B(1:3,1:4), B_T(1:4,1:3))
		allocate(T_node(1:num_node))
		
		! Initialize 
		Kc(:,:) = 0; Kh(:,:) = 0; Rq(:) = 0; Rh(:) = 0 
		
		
		
		! Read power density data 
		open(1, file="power.out",action="read", status="old")
		! real power
		do i_tet = 1, num_tet 
			read(1,*) Q(i_tet), Q_std(i_tet)
		enddo 
		close(1)
		
		
		
		
		! Tetrahedron loop
		do i_tet = 1, num_tet 
		
			if ( .not. in_the_list (T_partlist,Tet(i_tet)%part)) cycle 
			
			do i = 1, 3 
				do j = 1, 4 
					B  (i,j) = Tet(i_tet)%abcd(i+1,j)
					B_T(j,i) = B(i,j)
				enddo 
			enddo 
			
			! get thermal conductivity (k_c) value
			call defparam('T',T_tet)
			call evalexpr(mat(Tet(i_tet)%mat)%k,k_c)
			
			Kc_temp(:,:) = 0.0d0
			do i = 1, 4 
				do j = 1, 4
					do k = 1, 3 
						Kc_temp(i,j) = Kc_temp(i,j) + k_c*B_T(i,k)*B(k,j)/(36.0*Tet(i_tet)%vol)
					enddo 
				enddo 
			enddo 
			
			area = 0
			Kh_temp(:,:) = 0.0d0
			face: do i = 1, 4 
				
				
				!! ¾ß¸Å (TODO) 
				!perimeter = .true. 
				!do j = 0, 2
				!	ii = mod(i+j, 4) + 1
				!	val = sqrt(node(tet(i_tet)%node(ii))%xyz(1)**2 + node(tet(i_tet)%node(ii))%xyz(2)**2)
				!	!if (Tet(i_tet)%mat == 3) print *, val 
				!	if ( abs(val - 0.47600) > 1.0d-3) then 
				!		perimeter = .false.
				!	endif
				!enddo 
			
				! calculate Kh if neighboring coolant
				if (Tet(i_tet)%neighbor(i) < 0 .and. perimeter == .true.) then 
				!if (Tet(Tet(i_tet)%neighbor(i))%mat == 5 .and. perimeter == .true.) then 
				!if (mat(Tet(Tet(i_tet)%neighbor(i))%mat)%coolant == .true. .and. mat(Tet(i_tet)%mat)%coolant == .false.) then 
				!if (Tet(i_tet)%neighbor(i) < 0) then 
					! calculate convective heat transfer coefficient (h)
					h = 0.96
					Tf = 300
					
					Kh_temp(:,:) = 1.0d0 
					do j = 1, 4 
						Kh_temp(j,j) = 2.0d0 
					enddo 
					
					k = 1
					do j = 1, 4 
						if (i == j) cycle 
						idx(k) = j
						k = k+1
					enddo 
					
					area = area_tri(node(Tet(i_tet)%node(idx(1)))%xyz,&
									node(Tet(i_tet)%node(idx(2)))%xyz,&
									node(Tet(i_tet)%node(idx(3)))%xyz)
									
					do j = 1, 3
						Rh(Tet(i_tet)%node(idx(j))) = Rh(Tet(i_tet)%node(idx(j))) + h*Tf*area*(1.0/3.0);
					enddo 
					
					Kh_temp(i,:) = 0 
					Kh_temp(:,i) = 0 
					
					exit face 
				endif
			
			
			
				if (Tet(i_tet)%neighbor(i) < 0) cycle face

				!! calculate Kh if neighboring coolant
				!if (mat(Tet(Tet(i_tet)%neighbor(i))%mat)%coolant) then 
				!	! calculate convective heat transfer coefficient (h)
				!	call defparam('T',T_tet)
				!	call evalexpr(mat(Tet(Tet(i_tet)%neighbor(i))%mat)%hf,h)
				!	
				!	Tf = mat(Tet(Tet(i_tet)%neighbor(i))%mat)%Tf
				!	
				!	Kh_temp(:,:) = 1.0d0 
				!	do j = 1, 4 
				!		Kh_temp(j,j) = 2.0d0 
				!	enddo 
				!	
				!	k = 1
				!	do j = 1, 4 
				!		if (i == j) cycle 
				!		idx(k) = j
				!		k = k+1
				!	enddo 
				!	
				!	area = area_tri(node(Tet(i_tet)%node(idx(1)))%xyz,&
				!					node(Tet(i_tet)%node(idx(2)))%xyz,&
				!					node(Tet(i_tet)%node(idx(3)))%xyz)
				!					
				!	do j = 1, 3
				!		Rh(Tet(i_tet)%node(idx(j))) = Rh(Tet(i_tet)%node(idx(j))) + h*Tf*area*(1.0/3.0);
				!	enddo 
				!	
				!	Kh_temp(i,:) = 0 
				!	Kh_temp(:,i) = 0 
				!	
				!	exit face 
				!endif
			enddo face
			
			val = (1.0/12.0)*h*area 
			Kh_temp(:,:) = Kh_temp(:,:) * val
			
			do i = 1, 4 
				ii = Tet(i_tet)%node(i)
				do j = 1, 4 
					jj = Tet(i_tet)%node(j)
					Kc(ii,jj) = Kc(ii,jj) + Kc_temp(i,j)
					Kh(ii,jj) = Kh(ii,jj) + Kh_temp(i,j)
				enddo 
				! Make Rq matrix
				Rq(ii) = Rq(ii) + Q(i_tet)*Tet(i_tet)%vol * 0.25d0 
			enddo 
			
			
		enddo
		
		! Union matrix 
		Kc(:,:) = Kc(:,:) + Kh(:,:)
		Rq(:) = Rq(:) + Rh(:)
		
		
		
		! Solve the system matrix
		call matrixSolverForHeatTransfer(Kc,Rq,T_node,num_node)
		
		
		! Print T_node 
		open(1, file="Temperature.out",action="write", status="replace")
		do i = 1, num_node
			!if (Kc(i,i) == 0) cycle 
			!if (T_node(i) == 0) cycle 
			!write(1,'(4es16.6)') Node(i)%xyz(:), T_node(i)
			
			if (T_node(i) == 0) then 
				write(1,'(4es16.6,I)') Node(i)%xyz(:), 0.0/0.0
			else 
				write(1,'(4es16.6,I)') Node(i)%xyz(:), T_node(i)
			endif 
		enddo 
		
		close(1)
		
		!open(1, file="Temperature_tri.out",action="write", status="replace")
		!do i = 1, num_tri 
		!	val = T_node(Tri(i)%node(1))*(1.0d0/Tri(i)%r(1)) &
		!		+ T_node(Tri(i)%node(2))*(1.0d0/Tri(i)%r(2)) &
		!		+ T_node(Tri(i)%node(3))*(1.0d0/Tri(i)%r(3))
		!		
		!	vol = (1.0d0/Tri(i)%r(1))+(1.0d0/Tri(i)%r(2))+(1.0d0/Tri(i)%r(3))
		!	
		!	write(1,'(4es16.6,I)') i, val / vol
		!enddo 
		!close(1)
		
		
		
		
		
		! deallocate matrices 
		
		deallocate(Q,Q_std)
		deallocate(Kc,Kh)
		deallocate(Rq, Rh)
		deallocate(Kc_temp,Kh_temp,Rh_temp) 
		deallocate(B, B_T)		
		
		
		
	end subroutine 
	
	
	
	subroutine matrixSolverForHeatTransfer(A,B,x,n) 
		real(kind=rp), intent(in) :: A(:,:), B(:)
		real(kind=rp), intent(inout) :: x(:) 
		integer, intent(in) :: n 
		
		integer :: i, j 
		integer :: cnt, i_cnt, j_cnt
		real(kind=rp), allocatable, dimension(:,:) :: A_act
		real(kind=rp), allocatable, dimension(:)	 :: B_act, x_act
		integer :: nrhs, lda, ldb, ierr
		integer, allocatable :: ipiv(:)
		
		cnt = 0 
		do i = 1, n 
			if (A(i,i) /= 0) cnt = cnt + 1
		enddo 
		
		allocate(A_act(1:cnt,1:cnt))
		allocate(B_act(1:cnt))
		allocate(x_act(1:cnt))
		allocate(ipiv(1:cnt))
		
		i_cnt = 0; 
		do i = 1, n 
			if (A(i,i) == 0) cycle
			i_cnt = i_cnt + 1
			j_cnt = 0; 
			do j = 1, n 
				if (A(j,j) == 0) cycle
				j_cnt = j_cnt + 1 
				A_act(i_cnt,j_cnt) = A(i,j)
			enddo 
			B_act(i_cnt) = B(i)
		enddo 
		
		nrhs = 1
		lda = cnt
		ldb = cnt
		
		if (rp == sp) then 
			call sgesv( cnt, nrhs, A_act, lda, ipiv, B_act, ldb, ierr );
		else 
			call dgesv( cnt, nrhs, A_act, lda, ipiv, B_act, ldb, ierr );
		endif 
		
		
		if( ierr > 0 ) then
			print *, "The diagonal element of the triangular factor of A," 
			print *, "U(x,x) is zero, so that A is singular", ierr, ierr 
			print *, "the solution could not be computed."
			stop
		endif 
		!call CG(A_act,B_act,x_act,cnt)
		!call gauss_2(A_act,B_act,x_act,cnt)
		
		
		
		i_cnt = 0
		do i = 1, n 
			if (A(i,i) == 0) then 
				x(i) = 0 				! node out of interest 
			else 
				i_cnt = i_cnt + 1	
				!x(i) = x_act(i_cnt)	! fill the output
				x(i) = B_act(i_cnt)	! fill the output
			endif 
		enddo 
		
		deallocate(A_act, B_act, ipiv)
		deallocate(x_act)
	end subroutine 	
	
	
	
	

	subroutine finalize_multiphysics() 
		if(allocated(Tet)) deallocate(Tet)
		if(allocated(T_node)) deallocate(T_node)
	end subroutine
	
	

	
	function in_the_list (list,item) result(in) 
		integer, intent(in) :: list(:)
		integer, intent(in) :: item
		logical :: in
		integer :: i, n
		
		in =.false. 
		n = size(list)
		do i = 1, n 
			if (item == list(i)) then 
				in = .true.
				exit 
			endif
		enddo 
		
	end function 
	
	
	
	
	
	subroutine read_multiphysics_input
		integer :: nargs, nmatargs,ios
		character(maxstring) :: args(100), matargs(100)
		integer :: i, imat
		integer :: n
		real(kind=rp) :: val 
		
		call readandparse(rd_mp, args, nargs,io_input) 
		
		if (io_input /= 0) return 
				
		select case(uppercase(args(1)))
		case ('MAT')
			call value(args(2),imat,ios)
			!call removesp(matargs(1))
			n = 0 
			READMAT: do  
				n = n + 1
				call readandparse(rd_mp, matargs, nmatargs,io_input)
				select case (uppercase(matargs(1)))
				case('NAME') 
					mat(imat)%name = trim(matargs(2))
				case('K')
					if (trim(matargs(2)) == '=') then 
						mat(imat)%k = '' 
						do i = 3, nmatargs
							call insertstr(mat(imat)%k,trim(matargs(i)),index(mat(imat)%k,' '))
						enddo 
					else 
						mat(imat)%k    = trim(matargs(2))
					endif 
				case('STRAIN') 
					if (trim(matargs(2)) == '=') then 
						mat(imat)%strain = '' 
						do i = 3, nmatargs
							call insertstr(mat(imat)%strain,trim(matargs(i)),index(mat(imat)%strain,' '))
						enddo 
					else 
						mat(imat)%strain = trim(matargs(2))
					endif 
				case('E') 
					if (trim(matargs(2)) == '=') then 
						mat(imat)%E = '' 
						do i = 3, nmatargs
							call insertstr(mat(imat)%E,trim(matargs(i)),index(mat(imat)%E,' '))
						enddo 
					else 
						mat(imat)%E    = trim(matargs(2))
					endif 
				
				case('G')
					if (trim(matargs(2)) == '=') then 
						mat(imat)%G = '' 
						do i = 3, nmatargs
							call insertstr(mat(imat)%G,trim(matargs(i)),index(mat(imat)%G,' '))
						enddo 
					else 
						mat(imat)%G    = trim(matargs(2))
					endif 
				case('HF')
					if (trim(matargs(2)) == '=') then 
						mat(imat)%hf = '' 
						do i = 3, nmatargs
							call insertstr(mat(imat)%hf,trim(matargs(i)),index(mat(imat)%hf,' '))
						enddo 
					else 
						mat(imat)%hf    = trim(matargs(2))
					endif 
				case('COOLANT')
					read(matargs(2),'(l)') mat(imat)%coolant
				case('TF')
					if (trim(matargs(2))=='=') then 
						read(matargs(3),*) mat(imat)%Tf
					else 
						read(matargs(2),*) mat(imat)%Tf
					endif 
					
					
				case('DENSITY')
					if (trim(matargs(2))=='=') then 
						read(matargs(3),*) mat(imat)%rho
					else 
						read(matargs(2),*) mat(imat)%rho
					endif 
				case('C_P')
					if (trim(matargs(2))=='=') then 
						read(matargs(3),*) mat(imat)%c_p
					else 
						read(matargs(2),*) mat(imat)%c_p
					endif 
					
					
				case('ENDMAT')
					exit READMAT
				case('')
				case default 
					read(matargs(3),*) val
					call defparam(matargs(1),val)
				end select 
			enddo READMAT
			
		case('TEMPERATURE') 
			read(args(2),*) do_T 
			if (do_T) then 
				call readandparse(rd_mp, args, nargs,io_input)
				allocate(T_partlist(1:nargs)) 
				do i = 1, nargs
					read(args(i),*) T_partlist(i)
				enddo 
			endif
			
		case('THERMOMECHANIC')
			read(args(2),*) do_TM
			if (do_TM) then 
				call readandparse(rd_mp, args, nargs,io_input)
				allocate(TM_partlist(1:nargs)) 
				do i = 1, nargs
					read(args(i),*) TM_partlist(i)
				enddo 
			endif
			
			
		case('TRANSIENT') 
			read(args(2),*) do_transient
			if (do_transient) then 
				call readandparse(rd_mp, args, nargs,io_input)
				read(args(1),*) t_tot
				read(args(2),*) dt
				n_timestep = t_tot / dt
			endif
			
			
		case default
			write(*,'(a7,I2,a19,a)'), "(Line ", curr_line,") NO SUCH OPTION - ", args(1)
			stop
		end select
		
		
		
		
	end subroutine read_multiphysics_input
	
	subroutine readandparse(rdf, args, nargs, ierr) 
		integer, intent(in) :: rdf
		character(len=*),dimension(:),intent(inout) :: args
		integer,intent(inout) :: nargs
		integer,intent(inout) :: ierr
		character(len=132) :: str, line
		character(1) :: delims=' '
		integer :: idx, lenline
		
		str(:) = ''
		do 
			idx = len_trim(str)
			call readline(rdf,line,ierr,curr_line)
			lenline=len_trim(line)
			if (lenline==0) exit
			str(idx+1:idx+lenline) = line(1:lenline)
			if (line(lenline:lenline) /= '\') exit
		enddo 
		call parse(str,delims,args,nargs)
			
	end subroutine
	
	
	integer function findloc (array, element) 
		implicit none 
		
		integer, dimension(:), intent(in) :: array 
		integer, intent(in) :: element
		integer :: i, n 
		
		n = size (array) 
		
		findloc = 0 
		do i = 1, n 
			if (array(i) == element) then 
				findloc = i 
				return 
			endif 
		enddo 
		
	end function 
	
	
	
	
end module multiphysics