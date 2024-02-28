module ace_module

use variables
use constants, only: K_B
use ace_header
use strings, only : parse

implicit none 

contains

subroutine set_ace_iso(iso, lib)

!==============================================================================
implicit none
integer, intent(in) :: iso
character(*), intent(in) :: lib
integer :: i, j, k
integer :: line
integer :: iso0K
integer :: pt1, pt2
integer :: anum, mnum
integer :: ix, lmt 
integer, parameter :: ace_read_handler = 20171116
integer :: ierr

integer :: min_egrid, max_egrid, num_egrid, loc1, loc2, loc3, max_len
character(100) :: tmp1, args(100), tmparg
integer :: nargs, length, liblen
 
if(E_mode==0) return

!if(icore==score) print *, "   Start reading ace-format nuclear data"

!Ace-format data
  !Open ace format library
  if(icore==score) print *, iso, trim(ace(iso)%library)

  open(ace_read_handler, file=trim(ace(iso)%library), action="read")
  
  ! excited state if mnum > 300
  !read(ace(iso)%library(1:5), '(f)') val
  !i = mod(int(val), 1000)
!  write(*,*) iso, (ace(iso)%library)
!  if(ace(iso)%library(5:5)=='.') then
!      read(ace(iso)%library(1:4),'(i)') val
!  else
!      read(ace(iso)%library(1:5),'(i)') val
!  endif
!  i = mod(val,1000)
!  if (i > 300) then 
!    ace(iso)%excited = .true.
!  else 
!    ace(iso)%excited = .false. 
!  endif 
  liblen = len(trim(lib))
  do
    read(ace_read_handler,*) tmp1
    call parse(tmp1,' ',args, nargs)
    tmparg = adjustl(trim(args(1)))
    length = len(trim(tmparg))
    
    if(tmparg(length-2:length) == lib(liblen-2:liblen)) exit
  enddo
  backspace(ace_read_handler)
  !1st line
  read(ace_read_handler,'(i6, 4X, f12.6, es12.4)') ace(iso)%ZAID, ace(iso)%atn, ace(iso)%temp
  i = mod(ace(iso)%ZAID,1000)
  if (i > 300) then
      ace(iso)%excited = .true.
  else
      ace(iso)%excited = .false.
  endif

  !2~6 line
  Do line=2,6
    read(ace_read_handler,*)  
  Enddo
  
  
  !7~12 line
  read(ace_read_handler,10) (ace(iso)%NXS(i), i=1,16), (ace(iso)%JXS(i),i=1,32)
  !Allocate XSS array 
  if( allocated ( XSS ) ) deallocate ( XSS )
  allocate( XSS( 1 : ace(iso)%NXS(1)+4 ) )
  
  
  !Read XSS array 
  pt1 = 1
  pt2 = 4
  do
    if (pt2 >= ace(iso)%NXS(1)) then
      read(ace_read_handler,*) XSS(pt1:ace(iso)%NXS(1))
      exit
    end if
    read(ace_read_handler,*) XSS(pt1:pt2)
    pt1 = pt1 + 4
    pt2 = pt2 + 4
  end do
  !print *,iso, trim(ace(iso)%library), ace(iso)%NXS(3), ace(iso)%JXS(21)+1, xss(ace(iso)%jxs(21))
  

  call set_ESZ( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_NU(  iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_MTR( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_LQR( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_TYR( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_SIG( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_AND( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_DLW( iso, ace(iso)%NXS, ace(iso)%JXS, 1 )
  call set_DLW( iso, ace(iso)%NXS, ace(iso)%JXS, 2 )
  call set_DLW( iso, ace(iso)%NXS, ace(iso)%JXS, 3 )
  call set_NYD( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_FIS( iso, ace(iso)%NXS, ace(iso)%JXS )
  if(do_ures) call set_UNR(iso, ace(iso)%NXS, ace(iso)%JXS )

  !Maximum energy in XSS table
  if( Emax < XSS(ace(iso)%NXS(3)) ) Emax = XSS(ace(iso)%NXS(3))
  
  !Set recoverable energy per fission [MeV]
  !Ref. Eq. (7.19) in Nuclear Engineering Fundamentals: A Practical Perspective
  !First proposed by Unik and Ginlder [1970]
  
  allocate( ace(iso) % sigqf ( 1 : ace(iso)%NXS(3) ) ); ace(iso) % sigqf = 0d0
  !if(ace(iso)%JXS(21)/=0) then
      do i = 1, ace(iso) % NXS(4)
        if(ace(iso) % MT(i) == 18 .or. ace(iso) % TY(i) == 19) then
            ace(iso) % sigqf(:) = ace(iso) % sigqf(:) + ace(iso) % sig_MT(i) % cx(:) * ace(iso) % Q(i) * 202.27d0 / 193.4054d0
        endif
      enddo
  !endif

  if(ace(iso)%JXS(21)/=0) then
        ace(iso) % qval = sum( ace(iso) % sigqf ) / sum( ace(iso) % sigf )
  else
        ace(iso)%qval = 0.d0
  end if

  !Find (n,g) reaction cross section location (MT=102 or ENDF_NG) for burnup
  !JXS(28) is used to store location onum)%frac
  ace(iso)%JXS(28) = 0 !initialization of (n,g) pointer
  lmt = ace(iso)%JXS(3) - 1  !location of MTR Block
  do ix = 1, ace(iso)%NXS(4) !find MT=102 among reactions excluding elastic rx
    if( XSS( lmt + ix ) == ENDF_NG ) then
      ace(iso)%JXS(28) = ace(iso)%JXS(7) + &
&                 nint( XSS( ace(iso)%JXS(6) + ix - 1 ) ) 
      !XSS( JXS(28)-1 ) := first energy grid index for MT=102 (IE)
      !XSS( JXS(28) )   := number of consecutive energies for MT=102 (NE)
      !XSS( JXS(28)+IC ) := cross section for MT=102 at relative energy index IC
      !Since IE for (n,g) is 1, absolute energy grid index is used for IC.
      exit
    end if
  end do
  if(ace(iso)%JXS(28) == 0 ) then
    if(icore==score) print *, ace(iso)%ZAID, "(n,g) cross section is not found" 
  end if 

  !Deallocate XSS array for current isotope
  deallocate( XSS )
  !Close file handler for current isotope
  close(ace_read_handler)
10 format(8i9/8i9/8i9/8i9/8i9/8i9)

end subroutine set_ace_iso

subroutine set_sab_iso(iso)
  implicit none
  integer, intent(in) :: iso
  integer :: i, j, k
  integer :: line
  integer :: pt1, pt2
  integer :: anum, mnum
  integer :: ix, lmt 
  integer, parameter :: ace_read_handler = 20171116
  !real :: val
  integer :: val
  integer :: ierr
  
  integer :: min_egrid, max_egrid, num_egrid, loc1, loc2, loc3, max_len

  !Open ace format library
  if(icore==score) print *, '    S(a,b) lib: ', iso, trim(sab(iso)%library)

  open(ace_read_handler, file=trim(sab(iso)%library), action="read", iostat= ierr)
  
  !1~6 line
  Do line=1,6
    read(ace_read_handler,*)  
  Enddo
  
  
  !7~12 line
  read(ace_read_handler,10) sab(iso)%NXS(1:16), sab(iso)%JXS(1:32)
  !Allocate XSS array 
  if( allocated ( XSS ) ) deallocate ( XSS )
  allocate( XSS(1:sab(iso)%NXS(1)) )
  
  
  !Read XSS array 
  pt1 = 1
  pt2 = 4
  do
    if (pt2 >= sab(iso)%NXS(1)) then
      read(ace_read_handler,*) XSS(pt1:sab(iso)%NXS(1))
      exit
    end if
    read(ace_read_handler,*) XSS(pt1:pt2)
    pt1 = pt1 + 4
    pt2 = pt2 + 4
  end do
  
  call set_ITIE( iso, sab(iso)%NXS, sab(iso)%JXS )
  call set_ITCE( iso, sab(iso)%NXS, sab(iso)%JXS )
  call set_ITXE( iso, sab(iso)%NXS, sab(iso)%JXS )
  call set_ITCA( iso, sab(iso)%NXS, sab(iso)%JXS )

  deallocate(XSS)
  close(ace_read_handler)
10 format(8i9/8i9/8i9/8i9/8i9/8i9)

end subroutine set_sab_iso

subroutine set_ace0K_iso(iso0K)
!==============================================================================
implicit none
integer, intent(in) :: iso0K
integer :: i, j, k
integer :: line
integer :: pt1, pt2
integer :: anum, mnum
integer :: ix, lmt 
integer, parameter :: ace_read_handler = 20171116
!real :: val
integer :: val
integer :: ierr

integer :: min_egrid, max_egrid, num_egrid, loc1, loc2, loc3, max_len
 
if(E_mode==0) return

  !Open ace format library
  if(icore==score) print *, n_iso0K, trim(ace0K(iso0K)%library)
  open(ace_read_handler, file=trim(ace0K(iso0K)%library), action="read")
  
  ! 1st line
  read(ace_read_handler,11) ace0K(iso0K)%ZAID, ace0K(iso0K)%atn, ace0K(iso0K)%temp
  11 format(i6,4x,f12.6,es12.4)

  ! 2~6 line
  do line = 2, 6
    read(ace_read_handler,*)  
  enddo
  
  ! 7~12 line
  read(ace_read_handler,10) ace0K(iso0K)%NXS(1:16), ace0K(iso0K)%JXS(1:32)
  !Allocate XSS array 
  if( allocated ( XSS ) ) deallocate ( XSS )
  allocate( XSS(1:ace0K(iso0K)%NXS(1)) )

  
  !Set energy and elastic scattering cross section
  pt1 = 1
  pt2 = 4
  do
    if( pt2 > ace0K(iso0K)%NXS(3)*4 ) then
      read(ace_read_handler,*) XSS(pt1:ace0K(iso0K)%NXS(3)*4)
      exit
    end if
    read(ace_read_handler,*) XSS( pt1 : pt2 )
    pt1 = pt1 + 4
    pt2 = pt2 + 4
  end do

  call set_ESZ0K( iso0K, ace0K(iso0K)%NXS, ace0K(iso0K)%JXS )

  deallocate( XSS )
  close(ace_read_handler)

10 format(8i9/8i9/8i9/8i9/8i9/8i9)


end subroutine set_ace0K_iso

subroutine set_ace

!==============================================================================
implicit none
integer :: i, j, k
integer :: line
integer :: iso, iso0K
integer :: pt1, pt2
integer :: anum, mnum
integer :: ix, lmt 
integer, parameter :: ace_read_handler = 20171116
!real :: val
integer :: val
integer :: ierr

integer :: min_egrid, max_egrid, num_egrid, loc1, loc2, loc3, max_len
 
if(E_mode==0) return

if(icore==score) print *, "   Start reading ace-format nuclear data"

!Ace-format data
Emax = 0.d0
READ_ACE_ISO:Do iso = 1, num_iso
  !Open ace format library
  if(icore==score) print *, iso, trim(ace(iso)%library)

  open(ace_read_handler, file=trim(library_path)//trim(ace(iso)%library), action="read")
  
  ! excited state if mnum > 300
  !read(ace(iso)%library(1:5), '(f)') val
  !i = mod(int(val), 1000)
!  if(ace(iso)%library(5:5)=='.') then
!      read(ace(iso)%library(1:4),'(i)') val
!  else
!      read(ace(iso)%library(1:5),'(i)') val
!  endif
!  i = mod(val,1000)
!  if (i > 300) then 
!    ace(iso)%excited = .true.
!  else 
!    ace(iso)%excited = .false. 
!  endif 
  
  !1st line
  read(ace_read_handler,'(i6, 4X, f12.6, es12.4)') ace(iso)%ZAID, ace(iso)%atn, ace(iso)%temp

  !2~6 line
  Do line=2,6
    read(ace_read_handler,*)  
  Enddo
  
  
  !7~12 line
  read(ace_read_handler,10) (ace(iso)%NXS(i), i=1,16), (ace(iso)%JXS(i),i=1,32)
  !Allocate XSS array 
  if( allocated ( XSS ) ) deallocate ( XSS )
  allocate( XSS( 1 : ace(iso)%NXS(1)+4 ) )
  
  
  !Read XSS array 
  pt1 = 1
  pt2 = 4
  do
    if (pt2 >= ace(iso)%NXS(1)) then
      read(ace_read_handler,*) XSS(pt1:ace(iso)%NXS(1))
      exit
    end if
    read(ace_read_handler,*) XSS(pt1:pt2)
    pt1 = pt1 + 4
    pt2 = pt2 + 4
  end do
  !print *,iso, trim(ace(iso)%library), ace(iso)%NXS(3), ace(iso)%JXS(21)+1, xss(ace(iso)%jxs(21))
  

  call set_ESZ( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_NU(  iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_MTR( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_LQR( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_TYR( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_SIG( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_AND( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_DLW( iso, ace(iso)%NXS, ace(iso)%JXS, 1 )
  call set_DLW( iso, ace(iso)%NXS, ace(iso)%JXS, 2 )
  call set_DLW( iso, ace(iso)%NXS, ace(iso)%JXS, 3 )
  call set_NYD( iso, ace(iso)%NXS, ace(iso)%JXS )
  call set_FIS( iso, ace(iso)%NXS, ace(iso)%JXS )
  if(do_ures) call set_UNR(iso, ace(iso)%NXS, ace(iso)%JXS )

  
  !Maximum energy in XSS table
  if( Emax < XSS(ace(iso)%NXS(3)) ) Emax = XSS(ace(iso)%NXS(3))
  
  !Set recoverable energy per fission [MeV]
  !Ref. Eq. (7.19) in Nuclear Engineering Fundamentals: A Practical Perspective
  !First proposed by Unik and Ginlder [1970]
  
  if(ace(iso)%JXS(21)/=0) then
    anum = ace(iso)%NXS(2)/1000
    mnum = ace(iso)%NXS(2) - anum*1000
    !ace(iso)%qval = 1.29927d-3*(anum**2)*sqrt(dble(mnum))+33.12d0
    do i=1,ace(iso)%NXS(4)
        if(ace(iso)%MT(i)==18) then
        ace(iso)%qval = ace(iso)%Q(i)*1.0458343d0
        ! ADJUST Qval with Reference U-235 Heating (Serpent method)
!        if(icore==score) print *, 'KAPPA', anum*1000+mnum, ace(iso)%qval
        exit
        endif
    enddo
  else
    ace(iso)%qval = 0.d0
  end if

  !Find (n,g) reaction cross section location (MT=102 or ENDF_NG) for burnup
  !JXS(28) is used to store location onum)%frac
  ace(iso)%JXS(28) = 0 !initialization of (n,g) pointer
  lmt = ace(iso)%JXS(3) - 1  !location of MTR Block
  do ix = 1, ace(iso)%NXS(4) !find MT=102 among reactions excluding elastic rx
    if( XSS( lmt + ix ) == ENDF_NG ) then
      ace(iso)%JXS(28) = ace(iso)%JXS(7) + &
&                 nint( XSS( ace(iso)%JXS(6) + ix - 1 ) ) 
      !XSS( JXS(28)-1 ) := first energy grid index for MT=102 (IE)
      !XSS( JXS(28) )   := number of consecutive energies for MT=102 (NE)
      !XSS( JXS(28)+IC ) := cross section for MT=102 at relative energy index IC
      !Since IE for (n,g) is 1, absolute energy grid index is used for IC.
      exit
    end if
  end do
  if(ace(iso)%JXS(28) == 0 ) then
    if(icore==score) print *, ace(iso)%ZAID, "(n,g) cross section is not found" 
  end if 

  ! Create new absorption XS
  ace(iso) % siga(:) = ace(iso) % sigd(:)
  if(allocated(ace(iso) % sigf)) ace(iso) % siga(:) = ace(iso) % siga(:) +  ace(iso) % sigf(:)

  !Deallocate XSS array for current isotope
  deallocate( XSS )
  !Close file handler for current isotope
  close(ace_read_handler)

  !open(ace_read_handler, file="../ACE_293K/94236.70c_0293", action="read")
  !close(ace_read_handler)
  
End do READ_ACE_ISO

if ( sab_iso /= 0 ) then
if ( icore == score ) print *, "   Read S(a,b) scattering law tables"
READ_SAB_ISO : do iso = 1, sab_iso

  !Open ace format library
  if(icore==score) print *, iso, trim(sab(iso)%library)

  open(ace_read_handler, file=trim(library_path)//trim(sab(iso)%library), action="read", iostat= ierr)
  
  !1~6 line
  Do line=1,6
    read(ace_read_handler,*)  
  Enddo
  
  
  !7~12 line
  read(ace_read_handler,10) sab(iso)%NXS(1:16), sab(iso)%JXS(1:32)
  !Allocate XSS array 
  if( allocated ( XSS ) ) deallocate ( XSS )
  allocate( XSS(1:sab(iso)%NXS(1)) )
  
  
  !Read XSS array 
  pt1 = 1
  pt2 = 4
  do
    if (pt2 >= sab(iso)%NXS(1)) then
      read(ace_read_handler,*) XSS(pt1:sab(iso)%NXS(1))
      exit
    end if
    read(ace_read_handler,*) XSS(pt1:pt2)
    pt1 = pt1 + 4
    pt2 = pt2 + 4
  end do
  
  call set_ITIE( iso, sab(iso)%NXS, sab(iso)%JXS )
  call set_ITCE( iso, sab(iso)%NXS, sab(iso)%JXS )
  call set_ITXE( iso, sab(iso)%NXS, sab(iso)%JXS )
  call set_ITCA( iso, sab(iso)%NXS, sab(iso)%JXS )

  deallocate(XSS)
  close(ace_read_handler)

end do READ_SAB_ISO
end if



if ( n_iso0K /= 0 ) then
if ( icore == score ) print *, "   Read ace-format nuclear data at 0K"
READ_ACE0K_ISO : Do iso0K = 1, n_iso0K

  !Open ace format library
  if(icore==score) print *, iso0K, trim(ace0K(iso0K)%library)
  open(ace_read_handler, file=trim(library_path)//trim(ace0K(iso0K)%library), action="read")
  
  ! 1st line
  read(ace_read_handler,11) ace0K(iso0K)%ZAID, ace0K(iso0K)%atn, ace0K(iso0K)%temp
  11 format(i6,4x,f12.6,es12.4)

  ! 2~6 line
  do line = 2, 6
    read(ace_read_handler,*)  
  enddo
  
  ! 7~12 line
  read(ace_read_handler,10) ace0K(iso0K)%NXS(1:16), ace0K(iso0K)%JXS(1:32)
  !Allocate XSS array 
  if( allocated ( XSS ) ) deallocate ( XSS )
  allocate( XSS(1:ace0K(iso0K)%NXS(1)) )

  
  !Set energy and elastic scattering cross section
  pt1 = 1
  pt2 = 4
  do
    if( pt2 > ace0K(iso0K)%NXS(3)*4 ) then
      read(ace_read_handler,*) XSS(pt1:ace0K(iso0K)%NXS(3)*4)
      exit
    end if
    read(ace_read_handler,*) XSS( pt1 : pt2 )
    pt1 = pt1 + 4
    pt2 = pt2 + 4
  end do

  call set_ESZ0K( iso0K, ace0K(iso0K)%NXS, ace0K(iso0K)%JXS )

  deallocate( XSS )
  close(ace_read_handler)

end do READ_ACE0K_ISO
end if


10 format(8i9/8i9/8i9/8i9/8i9/8i9)


end subroutine set_ace



!==============================================================================

subroutine set_ESZ( iso, NXS, JXS )

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: pt1, pt2
type (AceFormat), pointer :: ac 

!Set pointer
ac => ace(iso)

!Allocate arrays 
allocate( ac % E( 1 : NXS(3) ) )
allocate( ac % sigt( 1 : NXS(3) ) )
allocate( ac % sigd( 1 : NXS(3) ) )
allocate( ac % sigel( 1 : NXS(3) ) )
allocate( ac % H( 1 : NXS(3) ) )
allocate( ac % siga( 1 : NXS(3) ) )

!Energies
pt1 = JXS(1)
pt2 = JXS(1) + NXS(3) - 1
ac % E( 1 : NXS(3) ) = XSS( pt1 : pt2 )

!Total cross sections
pt1 = pt2 + 1
pt2 = pt1 + NXS(3) - 1
ac % sigt( 1 : NXS(3) ) = XSS( pt1 : pt2 )

!Absorption cross sections
pt1 = pt2 + 1
pt2 = pt1 + NXS(3) - 1
ac % sigd( 1 : NXS(3) ) = XSS( pt1 : pt2 )

!Elastic scattering cross sections
pt1 = pt2 + 1
pt2 = pt1 + NXS(3) - 1
ac % sigel( 1 : NXS(3) ) = XSS( pt1 : pt2 )

!Average heating numbers
pt1 = pt2 + 1
pt2 = pt1 + NXS(3) - 1
ac % H( 1 : NXS(3) ) = XSS( pt1 : pt2 )


end subroutine set_ESZ



!==============================================================================
! SET_ESZ0K
!==============================================================================
subroutine set_ESZ0K( iso, NXS, JXS )
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: pt1, pt2
type (AceFormat0K), pointer :: ac0 
integer:: imin, imax


! Set pointer
ac0 => ace0K(iso)

! Allocate arrays 
allocate( ac0 % erg( 1 : NXS(3) ) )
allocate( ac0 % xs0( 1 : NXS(3) ) )
! Energies
pt1 = JXS(1)
pt2 = pt1 + NXS(3) - 1
ac0 % erg( 1 : NXS(3) ) = XSS( pt1 : pt2 )
! Elastic scattering cross sections
pt1 = JXS(1) + 3*NXS(3)
pt2 = pt1 + NXS(3) - 1
ac0 % xs0( 1 : NXS(3) ) = XSS( pt1 : pt2 )


!! Epithermal energy range
!call GET_IERG_DBRC(iso,imin,DBRC_E_min)
!call GET_IERG_DBRC(iso,imax,DBRC_E_max)
!
!
!! Reallocate
!deallocate(ac0%erg,ac0%xs0)
!allocate(ac0%erg(1:(imax-imin+1)))
!allocate(ac0%xs0(1:(imax-imin+1)))
!!   energy
!pt1 = jxs(1)+imin-1
!pt2 = pt1+(imax-imin)
!ac0%erg(1:) = XSS(pt1:pt2)
!!   elastic scattering cross sections
!pt1 = jxs(1)+3*nxs(3)+imin-1
!pt2 = pt1+(imax-imin)
!ac0%xs0(1:) = XSS(pt1:pt2)


end subroutine set_ESZ0K


! =============================================================================
! GET_IERG_DBRC
! =============================================================================
subroutine GET_IERG_DBRC(iso_,ierg_,erg)
    integer, intent(in)::  iso_
    integer, intent(out):: ierg_
    real(8), intent(in)::  erg
    type(AceFormat0K), pointer:: ac
    integer:: low, mid, high
    integer:: ne

    ac => ace0K(iso_)

!    ! beyond the energy boundary
!    ne = size(ac%erg(:))
!    if ( erg < ac%erg(1) ) then
!        ierg_ = 1
!        return
!    elseif ( erg > ac%erg(ne) ) then    
!        ierg_ = ne-1
!        return
!    end if

    ! binary search
    ne = size(ac%erg(:))
    low = 1
    high = ne
    if ( erg > ac%erg(ne) ) then
        ierg_ = ne-1
    else
    do while ( low+1 /= high ) 
        mid = (low+high)/2
        if ( erg < ac%erg(mid) ) then
            high = mid
        else
            low = mid
        end if
    end do
    ierg_ = low
    end if

    if ( associated(ac) ) nullify(ac)

end subroutine

!==============================================================================

subroutine set_NU( iso, NXS, JXS ) 

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: KNU
integer :: NC
integer :: NR, NE, IE
integer :: pcid
integer :: idx
type (AceFormat), pointer :: ac

!Set pointer
ac => ace(iso)

!check existence of NU block
if( JXS(2) == 0 ) then
  ac % nu_block_exist = .false.
  ac % nu_del_block_exist = .false.
  return
end if

!Set array structure flag for total nu
if( XSS( JXS(2) ) > 0 ) then
  KNU = JXS(2)
else
  KNU = JXS(2) + abs(XSS(JXS(2))) + 1 
end if
ac % nu_tot_flag = XSS( KNU )

select case( ac % nu_tot_flag )
case(1) !> polynomial function form
  ac % nu_tot % NC = XSS( KNU+1 )
  NC = ac % nu_tot % NC
  allocate( ac % nu_tot % C( 1 : NC ) )
  ac % nu_tot % C( 1 : NC ) = XSS( KNU+2 : KNU+NC+1 )
case(2) !> tabular data form
  ac % nu_tot % NR = XSS( KNU+1 )
  NR = ac % nu_tot % NR
  allocate( ac % nu_tot % NBT( 1 : NR ) )
  allocate( ac % nu_tot % INT( 1 : NR ) )
  ac % nu_tot % NBT( 1 : NR ) = XSS( KNU+2    : KNU+2+NR-1 )
  ac % nu_tot % INT( 1 : NR ) = XSS( KNU+2+NR : KNU+2+2*NR-1 )

  ac % nu_tot % NE = XSS( KNU+2+2*NR )
  NE = ac % nu_tot % NE
  allocate( ac % nu_tot % E( 1 : NE ) )
  allocate( ac % nu_tot % F( 1 : NE ) )
  ac % nu_tot % E( 1: NE ) = XSS( KNU+3+2*NR      : KNU+3+2*NR+NE-1 )
  ac % nu_tot % F( 1: NE ) = XSS( KNU+3+2*NR+NE   : KNU+3+2*NR+2*NE-1 )

    !do IE = 1, NE
    !    if (ac%library(1:5) == '92235') print *, IE, ac % nu_tot % E(IE), ac % nu_tot % F(IE)
    !enddo 
    !if (ac%library(1:5) == '92235') stop
  
end select

!set array structure flag for delayed nu
if( JXS(24) > 0 ) then
  KNU = JXS(24)
  ac % nu_del % NR = XSS( KNU+1 )
  NR = ac % nu_del % NR
  !> if NR = 0 :: linear-linear interpolation
  !> if NR > 0 :: ENDF interpolation parameters (NBT, INT) are used
  if( NR > 0 ) then
    allocate( ac % nu_del % NBT( 1 : NR ) )
    allocate( ac % nu_del % INT( 1 : NR ) )
    ac % nu_del % NBT( 1 : NR ) = XSS( KNU+2    : KNU+2+NR-1 )
    ac % nu_del % INT( 1 : NR ) = XSS( KNU+2+NR : KNU+2+2*NR-1)
  end if

  ac % nu_del % NE = XSS( KNU+2+2*NR )
  NE = ac % nu_del % NE
  allocate( ac % nu_del % E( 1 : NE ) )
  allocate( ac % nu_del % F( 1 : NE ) )
  ac % nu_del % E( 1 : NE ) = XSS( KNU+3+2*NR      : KNU+3+2*NR+NE-1 )
  ac % nu_del % F( 1 : NE ) = XSS( KNU+3+2*NR+NE   : KNU+3+2*NR+2*NE-1 )

  !Precursor Data :: decay constant and its probability
  allocate( ac % prcr( 1 : NXS(8) ) )
  idx = JXS(25)
  do pcid = 1, NXS(8)
    ac % prcr( pcid ) % decay_const = XSS( idx ) * 1.0d8
    ac % prcr( pcid ) % NR = XSS( idx+1 )
    NR = ac % prcr( pcid ) % NR
    if( NR > 0 ) then
      allocate( ac % prcr( pcid ) % NBT(1:NR) )
      allocate( ac % prcr( pcid ) % INT(1:NR) )
      ac % prcr( pcid ) % NBT(1:NR) = XSS( idx+2 : idx+2+NR-1 )
      ac % prcr( pcid ) % INT(1:NR) = XSS( idx+2+NR : idx+2+2*NR-1 )
    end if
    ac % prcr( pcid ) % NE = XSS( idx+2+2*NR )
    NE = ac % prcr( pcid ) % NE
    allocate( ac % prcr( pcid ) % E( 1 : NE ) )
    allocate( ac % prcr( pcid ) % F( 1 : NE ) )
    ac % prcr(pcid) % E( 1 : NE ) = XSS( idx+3+2*NR : idx+3+2*NR+NE-1 )
    ac % prcr(pcid) % F( 1 : NE ) = XSS( idx+3+2*NR+NE : idx+3+2*NR+2*NE-1 )
	
	idx = idx + (3+2*NR+2*NE)
	
  end do
  
else 
  ac % nu_del_block_exist = .false. 
end if
end subroutine set_NU


!==============================================================================

subroutine set_MTR( iso, NXS, JXS ) 

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)


!NXS(4) := number of reactions excluding elastic
if( NXS(4) == 0 ) return
allocate( ace(iso) % MT( 1 : NXS(4) ) )
ace(iso) % MT( 1 : NXS(4) ) = XSS( JXS(3) : JXS(3)+NXS(4)-1 )


end subroutine set_MTR


!==============================================================================
! List of kinetic Q-values for Reactions
subroutine set_LQR( iso, NXS, JXS ) 

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)


!NXS(4) := number of reactions excluding elastic
if( NXS(4) == 0 ) return
allocate( ace(iso) % Q( 1 : NXS(4) ) )
ace(iso) % Q( 1 : NXS(4) ) = XSS( JXS(4) : JXS(4)+NXS(4)-1 )


end subroutine set_LQR


!==============================================================================
! TYpe of Reactions
subroutine set_TYR( iso, NXS, JXS ) 

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: ii

!NXS(4) := number of reactions excluding elastic
if( NXS(4) == 0 ) return
allocate( ace(iso) % TY( 1 : NXS(4) ) )
ace(iso) % TY( 1 : NXS(4) ) = XSS( JXS(5) : JXS(5)+NXS(4)-1 )


end subroutine set_TYR


!==============================================================================

subroutine set_SIG( iso, NXS, JXS ) 

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
type (CrossSectionDataForm), pointer :: sigmt
integer :: ii   !> loop index for reaction MT
integer :: IE, NE
integer :: LOCA


!NXS(4) := number of reactions excluding elastic
if( NXS(4) == 0 ) return


allocate( ace(iso) % sig_MT( 1 : NXS(4) ) )
do ii = 1, NXS(4)
  !Set pointer
  sigmt => ace(iso) % sig_MT(ii)

  LOCA = nint( XSS( JXS(6)+ii-1 ) )
  
  sigmt % IE = nint( XSS( JXS(7)+LOCA-1 ) )
  IE = sigmt % IE
  sigmt % NE = nint( XSS( JXS(7)+LOCA ) )
  NE = sigmt % NE
  
  allocate ( sigmt % cx( 1 : IE+NE-1 ) )
  sigmt % cx(:) = 0 
  sigmt % cx( IE : IE+NE-1 ) = XSS( JXS(7)+LOCA+1 : JXS(7)+LOCA+NE )
end do


end subroutine set_SIG


!==============================================================================

subroutine set_AND( iso, NXS, JXS ) 

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
type (AngularDist), pointer :: an
integer :: ii   !> loop index for reaction MT
integer :: LOCB !> relative location of location of angualr distribution data for reactio MT(ii)
integer :: LC   !> location of tables associated with energies E(IE)
integer :: NE   !> number of energies at which anuglar distributions are tabulated
integer :: IE   !> energy grid index
integer :: loc  !> temporary locator
integer :: NP   !> number of points in the distribution

allocate( ace(iso) % ang( 0 : NXS(5) ) )
allocate( ace(iso) % ang_flag( 0 : NXS(5) ) )

do ii = 0, NXS(5)
  !> ii = 0 :: elastic scattering
  !> ii > 0 :: reaction MT(ii)
  LOCB = XSS( JXS(8) + ii )
  ace(iso) % ang_flag(ii) = LOCB

  !Set pointer
  an => ace(iso) % ang(ii)

  if( LOCB == 0 ) then
    !> LOCB == 0 :: no angular distributions are given, use isotropic distribution
    cycle
  else if( LOCB > 0 ) then
    NE = nint( XSS( JXS(9)+LOCB-1 ) )
    an % NE = NE
    allocate( an % dist_flag( 1: NE ) )  !> LC
    allocate( an % E( 1: NE ) )             !> E 
    allocate( an % dist(1 : NE ) )

    do IE = 1, NE
      an % E(IE) = XSS( JXS(9)+LOCB+IE-1 )
      LC = nint( XSS( JXS(9)+LOCB+NE+IE-1 ) )
      an % dist_flag(IE) = LC
      
      if( LC == 0 ) then
        !isotropic case
        cycle
      else if( LC > 0 ) then
        !32 equiprobable bin distribution
        allocate( an % dist(IE) % LDAT(1:33))
        loc = JXS(9) + LC - 1
        an % dist(IE) % LDAT(1:33) = XSS( loc : loc + 32)
        
      else 
        !tabular angular distribution given by PDF(1:NP) and CDF(1:NP)
        !loc = JXS(9) - LC - 1
        !NP = XSS(loc+1)
        !an % dist(IE) % NP = NP 
        !an % dist(IE) % JJ = XSS(loc)
        !
        !allocate( an % dist(IE) % CSOUT( 1 : NP ) )
        !allocate( an % dist(IE) % PDF( 1 : NP ) )
        !allocate( an % dist(IE) % CDF( 1 : NP ) )
        !
        !
        !!an % dist(IE) % LDAT( 1 : 3+3*NP-1 ) = XSS( loc : loc+3*3*NP-1 )
        !loc = loc + 2;  an % dist(IE) % CSOUT( 1 : NP ) = XSS( loc : loc+NP-1 )
        !loc = loc + NP; an % dist(IE) % PDF( 1 : NP )   = XSS( loc : loc+NP-1 )
        !loc = loc + NP; an % dist(IE) % CDF( 1 : NP )   = XSS( loc : loc+NP-1 )
        
        loc = JXS(9) - LC - 1
        NP = XSS(loc+1)
        allocate( an % dist(IE) % LDAT( 1 : 3+3*NP-1 ) )
        an % dist(IE) % LDAT( 1 : 3+3*NP-1 ) = XSS( loc : loc+3+3*NP-1-1 )
        !an % dist(IE) % LDAT(3:NP+2) = XSS(loc+2 : loc+2+NP-1)
      end if
    end do

  else if( LOCB == -1 ) then
    !> LOCB == -1 :: no angular distribution data are given for this reaction in the AND Block.
    !>               angular distribution data are specified through LAW(i)=44 in DLW Block.
    cycle
  else
    !> Fault
    call exception_handler( "set_AND", "LOCB flag is wrong" )
  end if
end do


end subroutine set_AND


!==============================================================================

subroutine set_DLW( iso, NXS, JXS, opt )

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer, intent(in) :: opt  !> option variables
                            !> opt=1 :: DLW block
                            !> opt=2 :: delayed neutron energy distribution 
                            !> opt=3 :: DLWP block 
type (EnergyDist), pointer :: eg
integer :: ii      !> loop index for reaction MT

integer :: LED     !> location of LDLW block
integer :: LDIS    !> location of DLW Block
integer :: NMT     !> number of reaction MT

integer :: LOCC    !> relative location of energy distribution array in DLW Block
integer :: nlaw    !> counter for number of law
integer :: ilaw    !> current law index
integer :: loc_law !> location of current law
integer :: LNW     !> location of next law
integer :: IDAT    !> location of law data array

integer :: NR      !> number of interpolation regions to define law applicability regime
integer :: NE      !> number of energies

select case( opt )
case( 1 ) !> prompt neutron energy distribution
  NMT = NXS(5)
  LED = JXS(10)
  LDIS = JXS(11)
  allocate( ace(iso) % pneg( 1 : NMT ) )
case( 2 ) !> delayed neutron energy distribution
  NMT = NXS(8)
  LED = JXS(26)
  LDIS = JXS(27)
  allocate( ace(iso) % dneg( 1 : NMT ) )
case( 3 ) !> prompt photon energy distribution 
  NMT = NXS(6)
  LED = JXS(18)
  LDIS = JXS(19)
  allocate( ace(iso) % ppeg( 1 : NMT ) )
end select


!NXS(6) !NXS(8) !NXS(5)
!LED = JXS(26)
!LDIS = JXS(27)

!LED = JXS(18)
!LDIS = JXS(19)


do ii = 1, NMT 
  nlaw = 1
  LOCC =  XSS( LED+ii-1 )

  !location of current law
  loc_law = LDIS+LOCC
  LNW = XSS( loc_law-1 )
  
  !set pointer
  select case( opt ) 
  case( 1 )
    eg => ace(iso) % pneg(ii)
  case( 2 ) 
    eg => ace(iso) % dneg(ii)
  case( 3 )
    eg => ace(iso) % ppeg(ii)
  end select 

  !count the number of laws for reaction MT(ii)
  do while( 1 )
    if( LNW  == 0 ) exit 
    nlaw = nlaw + 1
    loc_law = LDIS + LNW 
    LNW = XSS( loc_law-1 )
  enddo
  eg % nlaw = nlaw

  !allocate energy distribution array
  allocate( eg % dist( 1 : nlaw ) )

  !set enegy distribution array
  ilaw = 1
  loc_law = LDIS+LOCC
  LNW = XSS( loc_law-1 )
  do while( 1 )
    eg % dist(ilaw) % law           = XSS( loc_law )
    eg % dist(ilaw) % IDAT          = XSS( loc_law + 1)
    eg % dist(ilaw) % NR            = XSS( loc_law + 2)
    NR = eg % dist(ilaw) % NR !abbreviation
    if( NR > 0 ) then
      allocate( eg % dist(ilaw) % NBT( 1: NR ) )
      allocate( eg % dist(ilaw) % INT( 1: NR ) )
      eg % dist(ilaw) % NBT( 1 : NR ) = XSS( loc_law+3    : loc_law+3+NR-1   )
      eg % dist(ilaw) % INT( 1 : NR ) = XSS( loc_law+3+NR : loc_law+3+2*NR-1 )
    end if
    eg % dist(ilaw) % NE            = XSS( loc_law+3+2*NR)
    NE = eg % dist(ilaw) % NE !abbreviation
    allocate( eg % dist(ilaw) % E( 1: NE ) )
    allocate( eg % dist(ilaw) % F( 1: NE ) )
    eg % dist(ilaw) % E( 1 : NE )   = XSS( loc_law+4+2*NR    : loc_law+4+2*NR+NE-1 )
    eg % dist(ilaw) % F( 1 : NE )   = XSS( loc_law+4+2*NR+NE : loc_law+4+2*NR+2*NE-1 )
    
    !set energy distribution array for this law
    IDAT = XSS( loc_law+1 )
    call set_LAW( LDIS, LDIS+IDAT-1, eg % dist(ilaw) )

    !check last law   
    if( LNW == 0 ) exit

    !move to next law
    ilaw = ilaw + 1
    loc_law = LDIS+LNW
    LNW = XSS( loc_law-1 )
  end do
end do


end subroutine set_DLW


!==============================================================================

subroutine set_LAW( LDIS, loc_data, dist ) !> confirmed

!==============================================================================
implicit none
integer, intent(in) :: LDIS       !> offset of XSS array
integer, intent(in) :: loc_data   !> location of law data in XSS array
type (EnergyDistDataForm), intent(inout) :: dist !> selected energy distribution array to be set 
integer :: len  !> array length, LDAT(1:len)
integer :: K
integer :: IL
integer :: L
integer :: NR   
integer :: NE
integer :: NET
integer :: IP
integer :: NP, NP_ang
integer :: IE
integer :: ptr1, ptr2, ptr3
integer :: NRa, NRb, NEa, NEb
integer :: NPE, NPA
integer :: IPE, IPA
integer :: LC


select case( dist % law )

case( 1 )  !> Tabular Equiprobable Energy Bins (never called)

  write(err_msg, '(A)') 'Law 1 is not verified'
  call exception_handler( "set_LAW", err_msg )


  NR  = XSS( loc_data )
  NE  = XSS( loc_data + 1 + 2*NR )
  NET = XSS( loc_data + 2 + 2*NR + NE )
  len = 3 + 2*NR + NE + NET*NE
  allocate( dist % LDAT( 1 : len ) )
  dist % LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1 )


  !check LDAT array
!  print *, dist % LDAT(1), NR
!  print *, dist % LDAT(2+2*NR), NE
!  print *, dist % LDAT(3+2*NR+NE), NET
!  do IE = 1, NET
!    print *, IE, dist % LDAT(4+2*NR+NE+IE-1), XSS(loc_data + 2*NR+NE+IE-2 )
!  end do


case( 2 )  !> Discrete Photon Energy (confimred)
  allocate( dist % LDAT( 1 : 2 ) )
  dist % LDAT( 1 : 2 ) = XSS( loc_data : loc_data+1 )
  !check LDAT array
!  print *, dist % LDAT(1), dist % LDAT(2)
!  print *, XSS(loc_data), XSS(loc_data+1) 
!  stop


            
case( 3 )  !> Level Scattering (confimred)
  allocate( dist % LDAT( 1 : 2 ) )
  dist % LDAT( 1 : 2 ) = XSS( loc_data : loc_data+1 )

  !check LDAT array
!  print *, dist % LDAT(1), dist % LDAT(2)
!  print *, XSS(loc_data), XSS(loc_data+1) 
!  stop

case( 4 )  !> Continuous Tabular Distribution (confirmed)
  NR = XSS( loc_data )
  NE = XSS( loc_data + 1 + 2*NR )

  !get length of LDAT array
  len  = 2 + 2*NR + 2*NE
  do IE = 1, NE
    IL = loc_data + 1 + 2*NR + NE + IE
    K = LDIS + XSS(IL) - 1
    NP = XSS( K+1 )
    len = len + 2 + 3*NP
  end do
  allocate( dist % LDAT( 1 : len ) )
  
  !> My modification 
  dist%LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1)
  
  !> Jo's mistake ===================================================
!  !set LDAT array
!  ptr1 = 2 + 2*NR + NE
!  dist % LDAT( 1 : ptr1 ) = XSS( loc_data : loc_data+ptr1-1 )
!  ptr2 = 3 + 2*NR + 2*NE 
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NP = XSS( K+1 )
!    !set location of distribution and distribution data
!    dist % LDAT( ptr1+IE )  = ptr2
!    dist % LDAT( ptr2 : ptr2+2+3*NP-1 ) = XSS( K : K+2+3*NP-1 )
!
!    !update cumulative pointer
!    ptr2 = ptr2 + 2+3*NP
!  end do
  !> =================================================================

  !check LDAT array
!  print *, "NR", NR, dist % LDAT(1) 
!  print *, "NE", NE, dist % LDAT(2+2*NR) 
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NP = XSS( K+1 )
!    ptr1 = dist % LDAT(3 + 2*NR + NE + IE - 1)
!    print *, IE, dist % LDAT( ptr1 + 1), NP
!    do IP = 1, NP
!      print *, IE, IP, dist % LDAT( ptr1 + 2 + 2*NP + IP -1 ), XSS( K + 2 + 2*NP + IP - 1)
!    end do
!  end do
!  stop


case( 5 )  !> General Evaporation Spectrum (never called)

  write(err_msg, '(A)') 'Law 5 is not verified'
  call exception_handler( "set_LAW", err_msg )

  NR = XSS( loc_data )
  NE = XSS( loc_data + 1 + 2*NR )
  NET = XSS( loc_data + 2 + 2*NR + 2*NE )
  len = 3 + 2*NR + 2*NE + NET
  allocate( dist % LDAT( 1 : len ) )
  dist % LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1 )

! check LDAT array
!  print *, dist % LDAT(len), XSS(loc_data+len-1)
!  do IE = 1, NET
!    print *, IE, dist % LDAT( 4 + 2*NR + 2*NE + IE-1), XSS( loc_data + 3 + 2*NR + 2*NE + IE -1 )
!  end do
!  stop


case( 7 )  !> Simple Maxwell Fission Spectrum (confirmed)
  NR = XSS( loc_data )
  NE = XSS( loc_data + 1 + 2*NR )
  len = 3 + 2*NR + 2*NE
  allocate( dist % LDAT( 1 : len ) )
  dist % LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1 )
! check LDAT array
!  print *, dist % LDAT(len), XSS(loc_data+len-1)
!  do IE = 1, NE
!    print *, IE, dist % LDAT( 3 + 2*NR + NE + IE-1), XSS( loc_data + 2 + 2*NR + NE + IE-1 )
!  end do
!  stop


case( 9 )  !> Evaporation Spectrum (confirmed)
  NR = XSS( loc_data )
  NE = XSS( loc_data + 1 + 2*NR )
  len = 3 + 2*NR + 2*NE
  allocate( dist % LDAT( 1 : len ) )
  dist % LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1 )

! check LDAT array
!  print *, dist % LDAT(len), XSS(loc_data+len-1)
!  do IE = 1, NE
!    print *, IE, dist % LDAT( 3 + 2*NR + NE + IE-1), XSS( loc_data + 2 + 2*NR + NE + IE-1 )
!  end do
!  stop

case( 11 ) !> Energy Dependent Watt Spectrum (confirmed)
  NRa = XSS( loc_data )
  NEa = XSS( loc_data + 1 + 2*NRa )

  L = loc_data + 2 + 2*(NRa + NEa)
  NRb = XSS( L )
  NEb = XSS( L + 1 + 2*NRb )
  len = 5 + 2*(NRa + NEa + NRb + NEb) 
  allocate( dist % LDAT( 1 : len ) )
  dist % LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1 )

! check LDAT array
!  do IE = 1, NEb
!    print *, IE, dist % LDAT( 3 + 2*(NRa+NEa) + 2 + 2*NRb + NEb + IE-1 ), XSS( L + 2 + 2*NRb+NEb + IE - 1) 
!  end do
!  print *, dist % LDAT(len), XSS( loc_data+len-1 ) 
!  stop 


case( 22 ) !> Tabular Linear Functions 
  write(err_msg, '(A)') 'Law 22 is not prepared'
  call exception_handler( "set_LAW", err_msg )

 
case( 24 ) !> UK Law 6
  write(err_msg, '(A)') 'Law 24 is not prepared'
  call exception_handler( "set_LAW", err_msg )


case( 44 )  !> Kalbach-87 Formalism (confirmed)
  NR = XSS( loc_data )
  NE = XSS( loc_data + 1 + 2*NR )

  !get length of LDAT array
  len  = 2 + 2*NR + 2*NE
  do IE = 1, NE
    IL = loc_data + 1 + 2*NR + NE + IE
    K = LDIS + XSS(IL) - 1
    NP = XSS( K+1 )
    len = len + 2 + 5*NP
  end do
  allocate( dist % LDAT( 1 : len ) )
  
  
  !> My modification 
  dist%LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1)
  
  !> Jo's mistake ===================================================
!  !set LDAT array
!  ptr1 = 2 + 2*NR + NE
!  dist % LDAT( 1 : ptr1 ) = XSS( loc_data : loc_data+ptr1-1 )
!  
!  
!  print *, ptr1
!  print *, dist % LDAT( 1 : ptr1 )
!  
!  ptr2 = ptr1 + NE + 1
!  print *, ptr1, ptr2, NE
!  
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NP = XSS( K+1 )
!    !set location of distribution and distribution data
!    dist % LDAT( ptr1+IE )  = ptr2
!    dist % LDAT( ptr2 : ptr2+2+5*NP-1 ) = XSS( K : K+2+5*NP-1 )
!
!    !update cumulative pointer
!    ptr2 = ptr2 + 2+5*NP
!    
!   ! print '(4F10.4)', dist%ldat(3:100)
!    print *, dist % LDAT( ptr1+IE )
!    print *, dist % LDAT( ptr2 : ptr2+2+5*NP-1 )    
!    print *, ''
!    
!    print *, dist%ldat(1:3+2*NR+2*NE+2+4*NP)
!    stop 
!  end do
  
  ! ==================================================================

!  !check LDAT array
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NP = XSS( K+1 )
!
!    ptr1 = nint( dist % LDAT( 2 + 2*NR + NE + IE ) )
!    print *, NP
!    do IP = 1, NP
!      print *, "R", IP,  ( dist % LDAT( ptr1 + 1 + 3*NP + IP) ), ( XSS( K + 1 + 3*NP + IP  ) )
!    end do
!
!  end do


case( 61 ) !>  Tabular Angular Distribution
  NR = XSS( loc_data )
  NE = XSS( loc_data + 1 + 2*NR )

  !> My modification 
  K  = XSS(loc_data-1 + 3+2*NR+2*NE-1) - dist%IDAT+1
  NP = XSS(loc_data-1+K+1)
  LC = XSS(loc_data-1+K+2+4*NP-1)
  L = LDIS + abs(LC) -1 -1
  NP_ang = XSS(L+2)
  len = L+2+3*NP_ang - loc_data + 1

  allocate( dist % LDAT( 1 : len ) )
  dist%LDAT( 1 : len ) = XSS( loc_data : loc_data+len-1)
  
  
  !print *, len 
  !stop
  
  !> Jo's mistake ===================================================
  !get length of LDAT array
!  len = 2 + 2*NR + 2*NE
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NPE = XSS( K+1 )
!    len = len + 2 + 4*NPE
!    do IPE = 1, NPE
!      LC = XSS( K+1 + 3*NPE + IPE )
!      if( LC == 0 ) then
!        !isotropic distribution
!      else if( LC > 0 ) then
!        L = LDIS + abs(LC) - 1
!        NPA = XSS(L+1)
!        len = len + 2 + 3*NPA 
!      else 
!        write(err_msg, '(A)') 'Invalid angle distribution format in Law 61' 
!        call exception_handler( "set_LAW", err_msg )
!      end if
!    end do
!  end do
!  allocate( dist % LDAT( 1 : len ) )
!
!  !set LDAT array
!  ptr1 = 2 + 2*NR + NE
!  dist % LDAT( 1 : ptr1 ) = XSS( loc_data : loc_data+ptr1-1 )
!  ptr2 = ptr1 + NE + 1 
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NPE = XSS( K+1 )
!
!    !set location of energy distribution and energy distribution data
!    dist % LDAT( ptr1+IE ) = ptr2
!    dist % LDAT( ptr2 : ptr2+2 + 3*NPE -1 ) = XSS( K : K+2 + 3*NPE -1 )
!    ptr3 = ptr2 + 2 + 4*NPE 
!
!    !set location of angle distribution and angular distribution data
!    do IPE = 1, NPE
!      LC = XSS( K+1 + 3*NPE + IPE )
!
!      if( LC == 0 ) then
!        !isotropic distribution
!        dist % LDAT( ptr2 + 1 + 3*NPA + IPE ) = 0
!      else if( LC > 0 ) then
!        L = LDIS + abs(LC) - 1
!        NPA = XSS(L+1)
!        !set location of angle distribution and angular distribution data
!        dist % LDAT( ptr2 + 1 + 3*NPE + IPE ) = ptr3
!        dist % LDAT( ptr3 : ptr3+2 + 3*NPA -1) = XSS( L : L+2 + 3*NPA -1 )
!
!        !update cumulative pointer for angular distribution data
!        ptr3 = ptr3 + 2 + 3*NPA
!      else 
!        write(err_msg, '(A)') 'Invalid angle distribution format in Law 61' 
!        call exception_handler( "set_LAW", err_msg )
!      end if
!    end do
!
!    !update cumulative pointer for energy distribution data
!    ptr2 = ptr3 
!
!  end do
  !> =======================================================================
!  !check LDAT array
!  print *, "check1", nint(dist % LDAT(1)), NR
!  print *, "check2", nint(dist % LDAT(2+2*NR)), NE
!  do IE = 1, NE
!    IL = loc_data + 1 + 2*NR + NE + IE
!    K = LDIS + XSS(IL) - 1
!    NPE = XSS( K+1 )
!
!    print *, "check3", nint( dist % LDAT(  dist % LDAT(2+2*NR+NE + IE) ) ), nint( XSS(K) )
!    print *, "check4", nint( dist % LDAT(  dist % LDAT(2+2*NR+NE + IE) + 1 )), NPE 
!
!    !set location of energy distribution and energy distribution data
!    ptr2 = nint( dist % LDAT( 2 + 2*NR + NE + IE ) )
!
!    !set location of angle distribution and angular distribution data
!    do IPE = 1, NPE
!      LC = XSS( K+1 + 3*NPE + IPE )
!!      print '(4f15.2)', XSS(K+1 + IP), XSS(K+1 + NPE + IP), XSS( K+1 + 2*NPE + IP), XSS( K+1 + 3*NPE + IP)
!      if( LC == 0 ) then
!        !isotropic distribution
!      else if( LC > 0 ) then
!        L = LDIS + abs(LC) - 1
!        NPA = XSS(L+1)
!
!        ptr3 = dist % LDAT(ptr2+1 + 3*NPE + IPE)
!        print *, "check5", nint( dist % LDAT(ptr3+1) ), nint( XSS(L+1) )
!        do IPA = 1, NPA
!          print *, "check999", IPA, dist % LDAT(ptr3+1+2*NPA+IPA), XSS( L+1+2*NPA+IPA ) 
!        end do 
!      else 
!        write(err_msg, '(A)') 'Invalid angle distribution format in Law 61' 
!        call exception_handler( "set_LAW", err_msg )
!      end if
!    end do
!
!  end do


case( 66 ) !> N-body Phase Space Distribution (confirmed)
  allocate( dist % LDAT( 1 : 2 ) )
  dist % LDAT( 1 : 2 ) = XSS( loc_data : loc_data+1 )

  !check LDAT array
!  print *, dist % LDAT(1:2)
!  print *, XSS( loc_data : loc_data+1 )
!  stop

case( 67 ) !> Laboratory Angle-Energy Law
  write(err_msg, '(A)') 'Law 67 is not prepared'
  call exception_handler( "set_LAW", err_msg )


case default
  write(err_msg, '(i,A)') dist % law, ' Law can not be set' 
  call exception_handler( "set_LAW", err_msg )
end select

end subroutine set_LAW


!==============================================================================

subroutine set_NYD( iso, NXS, JXS ) !> confirmed

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: ii      !> loop index for reaction MT
integer :: JED     !> location of DLW block in XSS array
integer :: KY      !> location of neutron yield data in XSS array
integer :: NR      !> number of interpolation regions
integer :: NE      !> number of energies
type(TabularDataForm), pointer :: ny
integer :: IE

!NXS(4) := number of reactions excluding elastic
if( NXS(4) == 0 ) return

JED = JXS(11)

allocate( ace(iso) % nyd( 1 : NXS(4) ) )


do ii = 1, NXS(4)
  if( abs( ace(iso) % TY(ii) ) > 100 ) then  !> neutron yields Y(E) provieded as function of neutron energy
    !set pointer
    if ( icore == score ) print *, '         WARNING :: NYD block exists'
    
    ny => ace(iso) % nyd(ii)

    !find location of neutron yield data in XSS array
    KY = JED + abs( ace(iso) % TY(ii) ) - 101
    ny % NR = XSS( KY )
    NR = ny % NR !> abbreviation

    allocate( ny % NBT( 1 : NR ) )
    allocate( ny % INT( 1 : NR ) )

    ny % NBT(1:NR) = XSS( KY+1    : KY+NR   )
    ny % INT(1:NR) = XSS( KY+1+NR : KY+2*NR )

    ny % NE = XSS( KY+1+2*NR )
    NE = ny % NE !> abbreviation

    allocate( ny % E( 1 : NE ) )
    allocate( ny % F( 1 : NE ) )
  
    ny % E(1:NE) = XSS( KY+2+2*NR    : KY+2+2*NR+NE-1   )
    ny % F(1:NE) = XSS( KY+2+2*NR+NE : KY+2+2*NR+2*NE-1 )

    !do IE = 1, NE
    !  print *, IE, ny % E(IE), ny % F(IE)
    !end do
    
    !check neutron yield data
!    print *, ny % NE, NE
!    do IE = 1, NE
!      print *, IE, ny % E(IE), XSS(KY+1+2*NR+IE)
!    end do
!    do IE = 1, NE
!      print *, ny % F(IE), XSS(KY+1+2*NR+NE+IE)
!    end do


  else
    ace(iso) % nyd(ii) % NR = 0
    ace(iso) % nyd(ii) % NE = 0
  end if
end do


end subroutine set_NYD


!==============================================================================
!> Fis Block is a redundant quantity 

subroutine set_FIS( iso, NXS, JXS )

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)

integer :: IE, NE
type (AceFormat), pointer :: ac 
integer :: iMT, LOCA, nfis
real(8), allocatable :: temp_XS(:,:)

!Set pointer
ac => ace(iso)

!Allocate arrays 
if ( JXS(21) /= 0 ) then
    IE = XSS(JXS(21)); NE = XSS(JXS(21)+1);
    allocate( ac % sigf( 1 : NXS(3) ) )
    ac % sigf( : ) = 0 
    ac % sigf( IE : IE+NE-1 ) = XSS( JXS(21)+2 : JXS(21)+2+NE-1 )
    
else 
    nfis = 0 
    do iMT = 1, nxs(4)
        if (ac % ty(iMT) /= 19) cycle
        nfis = nfis+1
    enddo
    if (nfis == 0) return
    
    allocate(temp_XS(1:nfis, 1:NXS(3)))
    temp_XS(:,:) = 0 
    nfis = 0 
    do iMT = 1, nxs(4)
        if (ac % ty(iMT) /= 19) cycle
        nfis = nfis+1
        LOCA = XSS(JXS(6)+iMT-1)
        IE = XSS(JXS(7)+LOCA-1)
        NE = XSS(JXS(7)+LOCA)
        temp_XS(nfis,IE:IE+NE-1) = XSS( JXS(7)+LOCA+1 : JXS(7)+LOCA+NE )
    enddo 
    
    allocate( ac % sigf( 1 : NXS(3) ) )
    ac % sigf( : ) = 0 
    do IE = 1, NXS(3)
        ac % sigf (IE) = sum(temp_XS(:,IE)) 
        !if (ac%library(1:5) == '92235') print *, IE, ac%E(IE), ac%sigf(IE)
    enddo 
    !if (ac%library(1:5) == '92235') stop
endif




end subroutine set_FIS

!==============================================================================

subroutine set_UNR( iso, NXS, JXS )

!==============================================================================
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)

integer :: I, J, K, ii
integer :: pt1, pt2, N, M
type (AceFormat), pointer :: ac 

!Set pointer
ac => ace(iso)

!Set parameters

if(JXS(23)==0) return

ac % UNR % N   = XSS(JXS(23)  )
ac % UNR % M   = XSS(JXS(23)+1)
ac % UNR % INT = XSS(JXS(23)+2)
ac % UNR % ILF = XSS(JXS(23)+3)
ac % UNR % IOA = XSS(JXS(23)+4)
ac % UNR % IFF = XSS(JXS(23)+5)

if(ac % UNR % ILF > 0) then
    do ii = 1, ac % NXS(4)
        if(ac % MT(ii) == ac % UNR % ILF) then
        !if(ac % MT(ii) == 4) then
            ac % UNR % ILFidx = ii
            exit
        endif
    enddo
endif

if(ac % UNR % IOA > 0) then
    do ii = 1, ac % NXS(4)
        if(ac % MT(ii) == ac % UNR % IOA) then
            ac % UNR % IOAidx = ii
            exit
        endif
    enddo
endif

!Allocate arrays 
N = ac % UNR % N
M = ac % UNR % M

if(M<1D0 .or. N<1D0) return
!if(.not. ac % ZAID == 92238) return

n_unr = n_unr + 1
ac % UNR % unridx = n_unr

if(icore==score) write(*,*) '        UNR DATA FOUND for', ac % ZAID, n_unr

ac % UNR % URES = .true.
allocate( ac % UNR % E(1:N) )
allocate( ac % UNR % P(1:N, 1:6, 1:M) )

!Fill E 
ac % UNR % E (1:N) = XSS(JXS(23)+6:JXS(23)+5+N)

ac % UNR % Emin    = ac % UNR % E(1)
ac % UNR % Emax    = ac % UNR % E(N)

if(EUmin > ac % UNR % Emin) EUmin = ac % UNR % Emin
if(EUmax < ac % UNR % Emax) EUmax = ac % UNR % Emax

!Fill P table
pt1 = JXS(23)+6+N;
pt2 = pt1 + M - 1
do I = 1, N
    !write(*,*) ac%zaid, I, ac%UNR%E(I)
    do J = 1, 6 
        ac % UNR % P (I,J,1:M) = XSS(pt1:pt2)
        pt1 = pt2 + 1
        pt2 = pt1 + M - 1
        !write(*,'(I2,I2,<M>F10.4)') I, J, (ac%UNR%P(I,J,ii), ii=1,M)
    enddo 
enddo

if(icore==score) print *, ac % UNR % Emin, ac % UNR % Emax

end subroutine set_UNR



!==============================================================================
! set_ITIE reads energy-dependent inelastic scattering cross section
!==============================================================================
subroutine set_ITIE( iso, NXS, JXS )
implicit none
type (SAB_INEL_XS), pointer :: ab
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: pt1, pt2
integer:: ne

!Set pointer
ab => sab(iso)%itie

!Allocate arrays 
ab % ne = xss(jxs(1))
ne = ab % ne
allocate( ab % erg(1:ne) , ab % xs(1:ne) )

!Energies
pt1 = JXS(1)+1
pt2 = JXS(1)+ne  
ab % erg (1:ne) = xss(pt1:pt2)

! Inelastic cross section
pt1 = pt2+1
pt2 = pt2+ne
ab % xs (1:ne) = xss(pt1:pt2)

if ( associated(ab) ) nullify(ab)

end subroutine set_ITIE


!==============================================================================
! set_ITCE reads energy-dependent elastic scattering cross section
!==============================================================================
subroutine set_ITCE( iso, NXS, JXS )
implicit none
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: pt1, pt2
integer:: ne
type (SAB_EL_XS), pointer :: ab

if ( jxs(4) == 0 ) return

! set pointer
ab => sab(iso)%itce
print *, 'ACCESSED: ', iso

! allocate arrays 
ab % ne = xss(jxs(4))
ne = ab % ne
allocate( ab % erg(1:ne) , ab % xs(1:ne) )

! energies
pt1 = jxs(4)+1
pt2 = jxs(4)+ne  
ab % erg (1:ne) = xss(pt1:pt2)

! elastic cross section
pt1 = pt2+1
pt2 = pt2+ne
ab % xs (1:ne) = xss(pt1:pt2)

if ( associated(ab) ) nullify(ab)

end subroutine


!==============================================================================
! set_ITXE reads energy/angle distributions for inelastic scattering
!==============================================================================
subroutine set_ITXE( iso, NXS, JXS )
implicit none
type (SAB_INEL_E),  pointer:: ab
type (SAB_INEL_XS), pointer:: ie
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: pt1, pt2, pt3
integer:: ne1, ne2, na
integer:: i, j

! set pointer
ab => sab(iso)%itxe
ie => sab(iso)%itie

! allocate arrays 
ne1 = ie % ne    ! incoming energy
ne2 = nxs(4)     ! outgoing energy
na  = nxs(3)+1   ! outgoing angle
allocate( ab % erg(ne1,ne2) )
allocate( ab % ang(ne1,ne2,na) ) 

! energies & angles
do i = 1, ne1
do j = 1, ne2
    pt1 = jxs(3)+(nxs(3)+2)*(j-1)+(nxs(3)+2)*ne2*(i-1)
    pt2 = pt1 + nxs(3)+1

    ab%erg(i,j) = xss(pt1)
    ab%ang(i,j,:) = xss(pt1+1:pt2)
end do
end do

if ( nxs(2) /= 3 ) then
if ( icore == score ) then
    print*, "no equally-likely cosine for thermal scattering"
    stop
end if
end if

if ( associated(ab) ) nullify(ab)
if ( associated(ie) ) nullify(ie)

end subroutine


!==============================================================================
! set_ITCA reads anglar distributions for elastic scattering
!==============================================================================
subroutine set_ITCA( iso, NXS, JXS )
implicit none
type (SAB_EL_ANG),  pointer :: ab
type (SAB_INEL_XS), pointer :: ie
integer, intent(in) :: iso
integer, intent(in) :: NXS(1:16)
integer, intent(in) :: JXS(1:32)
integer :: pt1, pt2, pt3
integer:: ne, na
integer:: i, j


if ( jxs(4) == 0 .or. nxs(6) == -1 ) return

! set pointer
ab => sab(iso)%itca
ie => sab(iso)%itie

! allocate arrays 
ne = ie % ne    ! incoming energy
na = nxs(3)+1   ! outgoing angle
allocate( ab % ang(ne,na) ) 

! energies & angles
do i = 1, ne
    pt1 = jxs(6)+(nxs(6)+1)*(i-1)
    pt2 = pt1 + nxs(6)
    ab%ang(i,:) = xss(pt1:pt2)
end do

end subroutine

!  call set_LAND(iso, ace(iso)%NXS, ace(iso)%JXS )
!  call set_LDLW(iso, ace(iso)%NXS, ace(iso)%JXS )
!  call set_DLW(iso, ace(iso)%NXS, ace(iso)%JXS )
!  call set_YP(iso, ace(iso)%NXS, ace(iso)%JXS )


subroutine find_ACE(line, iso_, issab, dep)
    use variables, only: libpath, libname, libtemp
    character(*), intent(in) :: line
    logical, intent(in)  :: dep
    integer, intent(out) :: iso_
    logical, intent(out) :: issab
    integer :: i , j, length
    iso_ = 0
    issab= .false.
    do i = 1, num_iso
        if(trim(line) == trim(ace(i) % xslib)) then
            iso_ = i
            return
        endif
    enddo

    do i = 1, sab_iso
        if(trim(line) == trim(sab(i) % xslib)) then
            iso_ = i
            issab= .true.
        endif
    enddo

    ! libname: list of .80c, .81c, ...
    ! libpath: list of .711nc, .710nc, ...
    if(iso_ == 0) then
        do j = 1, size(libname)
            if(trim(line) == trim(libname(j))) then
                length = len_trim(line)
                if(libtemp(j) == 0d0) then
                    n_iso0K = n_iso0K + 1
                    ace0K(n_iso0K) % xslib   = trim(line)
                    ace0K(n_iso0K) % library = trim(libpath(j))
                    iso_ = n_iso0K
                    call set_ace0K_iso(iso_)
                else
                    select case(line(length:length))
                    case('c') ! XS
                        num_iso = num_iso + 1
                        ace(num_iso) % xslib   = trim(line)
                        ace(num_iso) % library = trim(libpath(j))
                        iso_ = num_iso
                        call set_ace_iso(iso_, trim(line))
                        acerecord(j) = .false.
                        if ( dep ) ace(num_iso) % depletable = .true.
                    case('t') ! SAB
                        sab_iso = sab_iso + 1
                        sab(sab_iso) % xslib   = trim(line)
                        sab(sab_iso) % library = trim(libpath(j))
                        sab(sab_iso) % temp    = libtemp(j) * K_B
                        iso_ = sab_iso
                        issab= .true.
                        call set_sab_iso(iso_)
                    end select
                endif
                return
            endif
        enddo
    endif

    if(iso_ == 0) write(*,*) 'Unknown Nuclide: ', line

end subroutine

end module
