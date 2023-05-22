module PCMFD
    use FMFD_HEADER
    !use VARIABLES, only: icore, score
    use VARIABLES   ! ***
    implicit none
    integer:: ix0, ix1, iy0, iy1, iz0, iz1
    contains

! =============================================================================
! OUT_OF_ZZ determines if a region is in or out of the zigzag boundary region
! =============================================================================
function OUT_OF_ZZ0(io,jo,ko)
    logical:: OUT_OF_ZZ0
    integer, intent(in):: io, jo, ko
    integer:: mo, no

    if  ( io == 1 .or. io == ncm(1) ) then
        OUT_OF_ZZ0 = .true.
        return
    end if
    if  ( jo == 1 .or. jo == ncm(2) ) then
        OUT_OF_ZZ0 = .true.
        return
    end if
    if  ( ko == 1 .or. ko == ncm(3) ) then
        OUT_OF_ZZ0 = .true.
        return
    end if
    
    if ( .not. zigzagon ) then
        OUT_OF_ZZ0 = .false.
        return
    end if

    do mo = 1, zz_div
    if ( zza0(mo) < io .and. io <= zza0(mo+1) ) then
        no = mo
        exit
    end if
    end do

    if ( zza1(no) < jo .and. jo <= zza2(no) ) then
        OUT_OF_ZZ0 = .false.
    else
        OUT_OF_ZZ0 = .true.
    end if

end function

! =============================================================================
! OUT_OF_ZZ determines if a region is in or out of the zigzag boundary region
! =============================================================================
function OUT_OF_ZZ(io,jo)
    logical:: OUT_OF_ZZ
    integer, intent(in):: io, jo
    integer:: mo, no
    
    if ( .not. zigzagon ) then
        OUT_OF_ZZ = .false.
        return
    end if
    OUT_OF_ZZ = .false.
    do mo = 1, zz_div
    if ( zzc0(mo) < io .and. io <= zzc0(mo+1) ) then
        no = mo
        exit
    end if
    end do

    if ( zzc1(no) < jo .and. jo <= zzc2(no) ) then
        OUT_OF_ZZ = .false.
    else
        OUT_OF_ZZ = .true.
    end if

end function

! =============================================================================
! OUT_OF_ZZ determines if a region is in or out of the zigzag boundary region
! plus one more region
! =============================================================================
function OUT_OF_ZZ1(io,jo)
    logical:: OUT_OF_ZZ1
    integer, intent(in):: io, jo
    integer:: mo, no

    if ( .not. zigzagon ) then
        OUT_OF_ZZ1 = .false.
        return
    end if
    
    do mo = 1, zz_div
    if ( zzc0(mo) < io .and. io <= zzc0(mo+1) ) then
        no = mo
        exit
    end if
    end do

    if ( zzc1(no)-1 < jo .and. jo <= zzc2(no)+1 ) then
        OUT_OF_ZZ1 = .false.
    else
        OUT_OF_ZZ1 = .true.
    end if

end function

! =============================================================================
! OUT_OF_ZZ determines if a region is in or out of the zigzag boundary region
! plus one more region
! =============================================================================
subroutine IN_ZZ(io,jo,ko,inz1,inz2)
    integer, intent(in):: io, jo, ko
    logical:: inz1, inz2
    integer:: mo, no

    if ( ko == 1 .or. ko == ncm(3) ) then
        inz1 = .true.; inz2 = .true.; return
    end if
    
    if ( io == 1 ) then
    if ( zzc1(1) < jo .and. jo <= zzc2(1) ) then
        inz1 = .true.; inz2 = .true.; return
    else if ( zzc1(1) == jo .or. jo == zzc2(1)+1 ) then
        inz1 = .false.; inz2 = .true.; return
    else
        inz1 = .false.; inz2 = .false.; return
    end if
    end if

    if ( io == ncm(1) ) then
    if ( zzc1(zz_div) < jo .and. jo <= zzc2(zz_div) ) then
        inz1 = .true.; inz2 = .true.; return
    else if ( zzc1(zz_div) == jo .or. jo == zzc2(zz_div)+1 ) then
        inz1 = .false.; inz2 = .true.; return
    else
        inz1 = .false.; inz2 = .false.; return
    end if
    end if

    do mo = 1, zz_div
    if ( zzc0(mo) < io .and. io <= zzc0(mo+1) ) then
        no = mo
        exit
    end if
    end do

    if ( zzc1(no)+1 < jo .and. jo <= zzc2(no)-1 ) then
        inz1 = .true.; inz2 = .false.
    else if ( zzc1(no)+1 == jo .or. jo == zzc2(no) ) then
        inz1 = .true.; inz2 = .true.
    else if ( zzc1(no) == jo .or. jo == zzc2(no)+1 ) then
        inz1 = .false.; inz2 = .true.
    else
        inz1 = .false.; inz2 = .false.
    end if

end subroutine


! =============================================================================
! L_PDHAT
! =============================================================================
subroutine L_PDHAT
    use VARIABLES, only: curr_cyc
    implicit none

    do ii = 2, nfm(1)
        fmDh(ii,:,:,1) = (fmJ0(ii,:,:,1)-5D-1 &
            *fmDt(ii,:,:,1)*(fphi1(ii,:,:)-fphi1(ii-1,:,:)))/fphi1(ii,:,:)
    end do
    do ii = 1, nfm(1)-1
        fmDh(ii,:,:,2) = (fmJ1(ii,:,:,2)+5D-1 &
            *fmDt(ii,:,:,2)*(fphi1(ii+1,:,:)-fphi1(ii,:,:)))/fphi1(ii,:,:)
    end do
    do jj = 2, nfm(2)
        fmDh(:,jj,:,3) = (fmJ0(:,jj,:,3)-5D-1 &
            *fmDt(:,jj,:,3)*(fphi1(:,jj,:)-fphi1(:,jj-1,:)))/fphi1(:,jj,:)
    end do
    do jj = 1, nfm(2)-1
        fmDh(:,jj,:,4) = (fmJ1(:,jj,:,4)+5D-1 &
            *fmDt(:,jj,:,4)*(fphi1(:,jj+1,:)-fphi1(:,jj,:)))/fphi1(:,jj,:)
    end do
    do kk = 2, nfm(3)
        fmDh(:,:,kk,5) = (fmJ0(:,:,kk,5)-5D-1 &
            *fmDt(:,:,kk,5)*(fphi1(:,:,kk)-fphi1(:,:,kk-1)))/fphi1(:,:,kk)
    end do
    do kk = 1, nfm(3)-1
        fmDh(:,:,kk,6) = (fmJ1(:,:,kk,6)+5D-1 &
            *fmDt(:,:,kk,6)*(fphi1(:,:,kk+1)-fphi1(:,:,kk)))/fphi1(:,:,kk)
    end do

    ! Boundary condition
    if ( .not. zigzagon ) then
    ii = 1;      fmDh(ii,:,:,1) = -fmJn(ii,:,:,1)/fphi1(ii,:,:)
    ii = nfm(1); fmDh(ii,:,:,2) = +fmJn(ii,:,:,2)/fphi1(ii,:,:)
    jj = 1;      fmDh(:,jj,:,3) = -fmJn(:,jj,:,3)/fphi1(:,jj,:)
    jj = nfm(2); fmDh(:,jj,:,4) = +fmJn(:,jj,:,4)/fphi1(:,jj,:)
    else
    do ii = 1, zz_div
        fmDh(zzf1(ii)+1,zzf0(ii)+1:zzf0(ii+1),:,1) = &
            -fmJn(zzf1(ii)+1,zzf0(ii)+1:zzf0(ii+1),:,1) &
            /fphi1(zzf1(ii)+1,zzf0(ii)+1:zzf0(ii+1),:)
        fmDh(zzf2(ii),zzf0(ii)+1:zzf0(ii+1),:,2) = &
            +fmJn(zzf2(ii),zzf0(ii)+1:zzf0(ii+1),:,2) &
            /fphi1(zzf2(ii),zzf0(ii)+1:zzf0(ii+1),:)
        fmDh(zzf0(ii)+1:zzf0(ii+1),zzf1(ii)+1,:,3) = &
            -fmJn(zzf0(ii)+1:zzf0(ii+1),zzf1(ii)+1,:,3) &
            /fphi1(zzf0(ii)+1:zzf0(ii+1),zzf1(ii)+1,:)
        fmDh(zzf0(ii)+1:zzf0(ii+1),zzf2(ii),:,4) = &
            +fmJn(zzf0(ii)+1:zzf0(ii+1),zzf2(ii),:,4) &
            /fphi1(zzf0(ii)+1:zzf0(ii+1),zzf2(ii),:)

        if ( zzf1(ii) /= 0 ) then
        fmDh(zzf1(ii),zzf0(ii)+1:zzf0(ii+1),:,2) = 0
        fmDh(zzf0(ii)+1:zzf0(ii+1),zzf1(ii),:,4) = 0
        end if
        if ( zzf2(ii) /= nfm(1) ) then
        fmDh(zzf2(ii)+1,zzf0(ii)+1:zzf0(ii+1),:,1) = 0
        fmDh(zzf0(ii)+1:zzf0(ii+1),zzf2(ii)+1,:,3) = 0
        end if

    end do
    end if
    kk = 1;      fmDh(:,:,kk,5) = -fmJn(:,:,kk,5)/fphi1(:,:,kk)
    kk = nfm(3); fmDh(:,:,kk,6) = +fmJn(:,:,kk,6)/fphi1(:,:,kk)

end subroutine

! =============================================================================
! L_PBC
! =============================================================================
subroutine L_PBC
    implicit none

    deltf1 = fmDh

    ! interface boundary
    do ii = 2, ncm(1); ix0 = 1+(ii-1)*fcr
        deltf1(ix0,:,:,1) = (fmDt(ix0,:,:,1)+fmDh(ix0,:,:,1)) &
                            / (1D0+2D0*fmDt(ix0,:,:,1))
    end do
    do ii = 1, ncm(1)-1; ix1 = ii*fcr
        deltf1(ix1,:,:,2) = (fmDt(ix1,:,:,2)+fmDh(ix1,:,:,2)) &
                            / (1D0+2D0*fmDt(ix1,:,:,2))
    end do
    do jj = 2, ncm(2); iy0 = 1+(jj-1)*fcr
        deltf1(:,iy0,:,3) = (fmDt(:,iy0,:,3)+fmDh(:,iy0,:,3)) &
                            / (1D0+2D0*fmDt(:,iy0,:,3))
    end do
    do jj = 1, ncm(2)-1; iy1 = jj*fcr
        deltf1(:,iy1,:,4) = (fmDt(:,iy1,:,4)+fmDh(:,iy1,:,4)) &
                            / (1D0+2D0*fmDt(:,iy1,:,4))
    end do
    do kk = 2, ncm(3); iz0 = 1+(kk-1)*fcz
        deltf1(:,:,iz0,5) = (fmDt(:,:,iz0,5)+fmDh(:,:,iz0,5)) &
                            / (1D0+2D0*fmDt(:,:,iz0,5))
    end do
    do kk = 1, ncm(3)-1; iz1 = kk*fcz
        deltf1(:,:,iz1,6) = (fmDt(:,:,iz1,6)+fmDh(:,:,iz1,6)) &
                            / (1D0+2D0*fmDt(:,:,iz1,6))
    end do

end subroutine

! =============================================================================
! L_MATRIX
! =============================================================================
subroutine L_PMATRIX
    implicit none
    real(8):: deno(nfm(1),nfm(2),nfm(3))  ! denominator

    !$omp parallel default(shared) private(ii,jj,kk)
    !$omp do
    do ii = 1, nfm(1)
        if ( ii /= 1 ) &
        Mfm(ii,:,:,3) = -(deltf0(ii,:,:,1)+deltf1(ii-1,:,:,2))/dfm(1)
        if ( ii /= nfm(1) ) &
        Mfm(ii,:,:,5) = -(deltf0(ii,:,:,2)+deltf1(ii+1,:,:,1))/dfm(1)
    end do
    !$omp end do
    !$omp do
    do jj = 1, nfm(2)
        if ( jj /= 1 ) &
        Mfm(:,jj,:,2) = -(deltf0(:,jj,:,3)+deltf1(:,jj-1,:,4))/dfm(2)
        if ( jj /= nfm(2) ) &
        Mfm(:,jj,:,6) = -(deltf0(:,jj,:,4)+deltf1(:,jj+1,:,3))/dfm(2)
    end do
    !$omp end do
    !$omp do
    do kk = 1, nfm(3)
        if ( kk /= 1 ) &
        Mfm(:,:,kk,1) = -(deltf0(:,:,kk,5)+deltf1(:,:,kk-1,6))/dfm(3)
        if ( kk /= nfm(3) ) &
        Mfm(:,:,kk,7) = -(deltf0(:,:,kk,6)+deltf1(:,:,kk+1,5))/dfm(3)
    end do
    !$omp end do
    !$omp end parallel

    Mfm(:,:,:,4) = &
    +(deltf0(:,:,:,1)+deltf1(:,:,:,1)+deltf0(:,:,:,2)+deltf1(:,:,:,2))/dfm(1) &
    +(deltf0(:,:,:,3)+deltf1(:,:,:,3)+deltf0(:,:,:,4)+deltf1(:,:,:,4))/dfm(2) &
    +(deltf0(:,:,:,5)+deltf1(:,:,:,5)+deltf0(:,:,:,6)+deltf1(:,:,:,6))/dfm(3) &
    +fm_a(:,:,:)


    if ( zigzagon ) then
    where( fphi1 == 0 )
        Mfm(:,:,:,1) = 0
        Mfm(:,:,:,2) = 0
        Mfm(:,:,:,3) = 0
        Mfm(:,:,:,4) = 0
        Mfm(:,:,:,5) = 0
        Mfm(:,:,:,6) = 0
        Mfm(:,:,:,7) = 0
    end where
    end if

    ! -------------------------------------------------------------------------
    !   source term
    !$omp parallel default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    !$omp do
    do ii = 2, ncm(1); ix0 = 1+(ii-1)*fcr
        deno(ix0,:,:) = 1D0+2D0*fmDt(ix0,:,:,1)
        jsrc(ix0,:,:,1) = 4D0*fmDt(ix0,:,:,1)/deno(ix0,:,:)
        fsrc(ix0,:,:,1) = fmDh(ix0-1,:,:,2)/deno(ix0,:,:)
    end do
    !$omp end do
    !$omp do
    do ii = 1, ncm(1)-1; ix1 = ii*fcr
        deno(ix1,:,:) = 1D0+2D0*fmDt(ix1,:,:,2)
        jsrc(ix1,:,:,2) = 4D0*fmDt(ix1,:,:,2)/deno(ix1,:,:)
        fsrc(ix1,:,:,2) = fmDh(ix1+1,:,:,1)/deno(ix1,:,:)
    end do
    !$omp end do
    !$omp do
    do jj = 2, ncm(2); iy0 = 1+(jj-1)*fcr
        deno(:,iy0,:) = 1D0+2D0*fmDt(:,iy0,:,3)
        jsrc(:,iy0,:,3) = 4D0*fmDt(:,iy0,:,3)/deno(:,iy0,:)
        fsrc(:,iy0,:,3) = fmDh(:,iy0-1,:,4)/deno(:,iy0,:)
    end do
    !$omp end do
    !$omp do
    do jj = 1, ncm(2)-1; iy1 = jj*fcr
        deno(:,iy1,:) = 1D0+2D0*fmDt(:,iy1,:,4)
        jsrc(:,iy1,:,4) = 4D0*fmDt(:,iy1,:,4)/deno(:,iy1,:)
        fsrc(:,iy1,:,4) = fmDh(:,iy1+1,:,3)/deno(:,iy1,:)
    end do
    !$omp end do
    !$omp do
    do kk = 2, ncm(3); iz0 = 1+(kk-1)*fcz
        deno(:,:,iz0) = 1D0+2D0*fmDt(:,:,iz0,5)
        jsrc(:,:,iz0,5) = 4D0*fmDt(:,:,iz0,5)/deno(:,:,iz0)
        fsrc(:,:,iz0,5) = fmDh(:,:,iz0-1,6)/deno(:,:,iz0)
    end do
    !$omp end do
    !$omp do
    do kk = 1, ncm(3)-1; iz1 = kk*fcz
        deno(:,:,iz1) = 1D0+2D0*fmDt(:,:,iz1,6)
        jsrc(:,:,iz1,6) = 4D0*fmDt(:,:,iz1,6)/deno(:,:,iz1)
        fsrc(:,:,iz1,6) = fmDh(:,:,iz1+1,5)/deno(:,:,iz1)
    end do
    !$omp end do
    !$omp end parallel

end subroutine

! =============================================================================
! L_PSOURCE
! =============================================================================
subroutine L_PSOURCE(keff)
    implicit none
    real(8), intent(in):: keff

    ! neutron fission source
    fm_s(:,:,:) = fm_nf(:,:,:)*fphi1(:,:,:)/keff

    ! interface BC
    !$omp parallel default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    !$omp do
    do ii = 2, ncm(1); ix0 = 1+(ii-1)*fcr
        fm_s(ix0,:,:) = fm_s(ix0,:,:) &
            +(jsrc(ix0,:,:,1)*fmJ1(ix0,:,:,1) &
            +fsrc(ix0,:,:,1)*fphi1(ix0-1,:,:))/dfm(1)
    end do
    !$omp end do
    !$omp do
    do ii = 1,ncm(1)-1; ix1 = ii*fcr
        fm_s(ix1,:,:) = fm_s(ix1,:,:) & 
            +(jsrc(ix1,:,:,2)*fmJ0(ix1,:,:,2) &
            +fsrc(ix1,:,:,2)*fphi1(ix1+1,:,:))/dfm(1)
    end do
    !$omp end do
    !$omp do
    do jj = 2, ncm(2); iy0 = 1+(jj-1)*fcr
        fm_s(:,iy0,:) = fm_s(:,iy0,:) &
            +(jsrc(:,iy0,:,3)*fmJ1(:,iy0,:,3) &
            +fsrc(:,iy0,:,3)*fphi1(:,iy0-1,:))/dfm(2)
    end do
    !$omp end do
    !$omp do
    do jj = 1, ncm(2)-1; iy1 = jj*fcr
        fm_s(:,iy1,:) = fm_s(:,iy1,:) &
            +(jsrc(:,iy1,:,4)*fmJ0(:,iy1,:,4) &
            +fsrc(:,iy1,:,4)*fphi1(:,iy1+1,:))/dfm(2)
    end do
    !$omp end do
    !$omp do
    do kk = 2, ncm(3); iz0 = 1+(kk-1)*fcz
        fm_s(:,:,iz0) = fm_s(:,:,iz0) &
            +(jsrc(:,:,iz0,5)*fmJ1(:,:,iz0,5) &
            +fsrc(:,:,iz0,5)*fphi1(:,:,iz0-1))/dfm(3)
    end do
    !$omp end do
    !$omp do
    do kk = 1, ncm(3)-1; iz1 = kk*fcz
        fm_s(:,:,iz1) = fm_s(:,:,iz1) &
            +(jsrc(:,:,iz1,6)*fmJ0(:,:,iz1,6) &
            +fsrc(:,:,iz1,6)*fphi1(:,:,iz1+1))/dfm(3)
    end do
    !$omp end do
    !$omp end parallel

    ! net current
!    if ( wholecore .and. zigzagon ) then
!    do ii = 1, zz_div
!        fm_s(zzb1(ii)+1,zzb0(ii)+1:zzb0(ii+1),:) = &
!            fm_s(zzb1(ii)+1,zzb0(ii)+1:zzb0(ii+1),:) + &
!            fmJ3(zzb1(ii)+1,zzb0(ii)+1:zzb0(ii+1),:,1)/dfm(1)
!        fm_s(zzb2(ii),zzb0(ii)+1:zzb0(ii+1),:) = &
!            fm_s(zzb2(ii),zzb0(ii)+1:zzb0(ii+1),:) - &
!            fmJ3(zzb2(ii),zzb0(ii)+1:zzb0(ii+1),:,2)/dfm(1)
!        fm_s(zzb0(ii)+1:zzb0(ii+1),zzb1(ii)+1,:) = &
!            fm_s(zzb0(ii)+1:zzb0(ii+1),zzb1(ii)+1,:) + &
!            fmJ3(zzb0(ii)+1:zzb0(ii+1),zzb1(ii)+1,:,3)/dfm(2)
!        fm_s(zzb0(ii)+1:zzb0(ii+1),zzb2(ii),:) = &
!            fm_s(zzb0(ii)+1:zzb0(ii+1),zzb2(ii),:) - &
!            fmJ3(zzb0(ii)+1:zzb0(ii+1),zzb2(ii),:,4)/dfm(2)
!    end do
!    fm_s(:,:,1+fcz) = fm_s(:,:,1+fcz) + fmJ3(:,:,1+fcz,5)/dfm(3)
!    fm_s(:,:,nfm(3)-fcz) = fm_s(:,:,nfm(3)-fcz) - fmJ3(:,:,nfm(3)-fcz,6)/dfm(3)
!    end if
    
end subroutine

! =============================================================================
! L_OUTJ
! =============================================================================
subroutine L_POUTJ
    implicit none

    ! outgoing partial current
    !$omp parallel default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    !$omp do
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
        if ( ii /= 1 ) then
        fmJn(ix0,:,:,1) = +jsrc(ix0,:,:,1)*fmJ1(ix0,:,:,1) &
            -deltf1(ix0,:,:,1)*fphi1(ix0,:,:) &
            +fsrc(ix0,:,:,1)*fphi0(ix0-1,:,:)
        fmJ0(ix0,:,:,1) = fmJ1(ix0,:,:,1)-fmJn(ix0,:,:,1)
        end if
        if ( ii /= ncm(1) ) then
        fmJn(ix1,:,:,2) = -jsrc(ix1,:,:,2)*fmJ0(ix1,:,:,2) &
            +deltf1(ix1,:,:,2)*fphi1(ix1,:,:) &
            -fsrc(ix1,:,:,2)*fphi0(ix1+1,:,:)
        fmJ1(ix1,:,:,2) = fmJ0(ix1,:,:,2)+fmJn(ix1,:,:,2)
        end if
    end do
    !$omp end do
    !$omp do
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
        if ( jj /= 1 ) then
        fmJn(:,iy0,:,3) = +jsrc(:,iy0,:,3)*fmJ1(:,iy0,:,3) &
            -deltf1(:,iy0,:,3)*fphi1(:,iy0,:) &
            +fsrc(:,iy0,:,3)*fphi0(:,iy0-1,:)
        fmJ0(:,iy0,:,3) = fmJ1(:,iy0,:,3)-fmJn(:,iy0,:,3)
        end if
        if ( jj /= ncm(2) ) then
        fmJn(:,iy1,:,4) = -jsrc(:,iy1,:,4)*fmJ0(:,iy1,:,4) &
            +deltf1(:,iy1,:,4)*fphi1(:,iy1,:) &
            -fsrc(:,iy1,:,4)*fphi0(:,iy1+1,:)
        fmJ1(:,iy1,:,4) = fmJ0(:,iy1,:,4)+fmJn(:,iy1,:,4)
        end if
    end do
    !$omp end do
    !$omp do
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        if ( kk /= 1 ) then
        fmJn(:,:,iz0,5) = +jsrc(:,:,iz0,5)*fmJ1(:,:,iz0,5) &
            -deltf1(:,:,iz0,5)*fphi1(:,:,iz0) &
            +fsrc(:,:,iz0,5)*fphi0(:,:,iz0-1)
        fmJ0(:,:,iz0,5) = fmJ1(:,:,iz0,5)-fmJn(:,:,iz0,5)
        end if
        if ( kk /= ncm(3) ) then
        fmJn(:,:,iz1,6) = -jsrc(:,:,iz1,6)*fmJ0(:,:,iz1,6) &
            +deltf1(:,:,iz1,6)*fphi1(:,:,iz1) &
            -fsrc(:,:,iz1,6)*fphi0(:,:,iz1+1)
        fmJ1(:,:,iz1,6) = fmJ0(:,:,iz1,6)+fmJn(:,:,iz1,6)
        end if
    end do
    !$omp end do

    ! data swapping for updating next iteration incoming partial current
    !$omp do
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
        if ( ii /= 1 )      fmJ0(ix0-1,:,:,2) = fmJ0(ix0,:,:,1)
        if ( ii /= ncm(1) ) fmJ1(ix1+1,:,:,1) = fmJ1(ix1,:,:,2)
    end do
    !$omp end do
    !$omp do
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
        if ( jj /= 1 )      fmJ0(:,iy0-1,:,4) = fmJ0(:,iy0,:,3)
        if ( jj /= ncm(2) ) fmJ1(:,iy1+1,:,3) = fmJ1(:,iy1,:,4)
    end do
    !$omp end do
    !$omp do
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        if ( kk /= 1 )      fmJ0(:,:,iz0-1,6) = fmJ0(:,:,iz0,5)
        if ( kk /= ncm(3) ) fmJ1(:,:,iz1+1,5) = fmJ1(:,:,iz1,6)
    end do
    !$omp end do
    !$omp end parallel

end subroutine


! =============================================================================
! L_REFJ
! =============================================================================
subroutine L_PREFJ
    implicit none

    ! zigzag boundary
    if ( zigzagon ) then
    do ii = 1, zz_div
    if ( zzf1(ii) /= 0 ) then
        fmJn(zzf1(ii),zzf0(ii)+1:zzf0(ii+1),:,2) = &
            fmJn(zzf1(ii)+1,zzf0(ii)+1:zzf0(ii+1),:,1)
        fmJn(zzf0(ii)+1:zzf0(ii+1),zzf1(ii),:,4) = &
            fmJn(zzf0(ii)+1:zzf0(ii+1),zzf1(ii)+1,:,3)
    end if
    if ( zzf2(ii) /= nfm(1) ) then
        fmJn(zzf2(ii)+1,zzf0(ii)+1:zzf0(ii+1),:,1) = &
            fmJn(zzf2(ii),zzf0(ii)+1:zzf0(ii+1),:,2)
        fmJn(zzf0(ii)+1:zzf0(ii+1),zzf2(ii)+1,:,3) = &
            fmJn(zzf0(ii)+1:zzf0(ii+1),zzf2(ii),:,4)
    end if
    end do
    end if

    ! boundary surface
    ii = 1;      fmJn(ii,:,:,1) = -deltf1(ii,:,:,1)*fphi1(ii,:,:)
    ii = nfm(1); fmJn(ii,:,:,2) = +deltf1(ii,:,:,2)*fphi1(ii,:,:)
    jj = 1;      fmJn(:,jj,:,3) = -deltf1(:,jj,:,3)*fphi1(:,jj,:)
    jj = nfm(2); fmJn(:,jj,:,4) = +deltf1(:,jj,:,4)*fphi1(:,jj,:)
    kk = 1;      fmJn(:,:,kk,5) = -deltf1(:,:,kk,5)*fphi1(:,:,kk)
    kk = nfm(3); fmJn(:,:,kk,6) = +deltf1(:,:,kk,6)*fphi1(:,:,kk)

    ! net current & surface average
    if ( zigzagon ) then
    !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        
        ! inner surfaces
        cmJ0(ii,jj,kk,1) = sum(fmJ0(ix0,iy0:iy1,iz0:iz1,1))/fc1
        cmJ1(ii,jj,kk,1) = sum(fmJ1(ix0,iy0:iy1,iz0:iz1,1))/fc1

        cmJ0(ii,jj,kk,2) = sum(fmJ0(ix1,iy0:iy1,iz0:iz1,2))/fc1
        cmJ1(ii,jj,kk,2) = sum(fmJ1(ix1,iy0:iy1,iz0:iz1,2))/fc1

        cmJ0(ii,jj,kk,3) = sum(fmJ0(ix0:ix1,iy0,iz0:iz1,3))/fc1
        cmJ1(ii,jj,kk,3) = sum(fmJ1(ix0:ix1,iy0,iz0:iz1,3))/fc1

        cmJ0(ii,jj,kk,4) = sum(fmJ0(ix0:ix1,iy1,iz0:iz1,4))/fc1
        cmJ1(ii,jj,kk,4) = sum(fmJ1(ix0:ix1,iy1,iz0:iz1,4))/fc1

        cmJ0(ii,jj,kk,5) = sum(fmJ0(ix0:ix1,iy0:iy1,iz0,5))/fc2
        cmJ1(ii,jj,kk,5) = sum(fmJ1(ix0:ix1,iy0:iy1,iz0,5))/fc2

        cmJ0(ii,jj,kk,6) = sum(fmJ0(ix0:ix1,iy0:iy1,iz1,6))/fc2
        cmJ1(ii,jj,kk,6) = sum(fmJ1(ix0:ix1,iy0:iy1,iz1,6))/fc2
    
        ! boundary surfaces
        cmJn(ii,jj,kk,1) = sum(fmJn(ix0,iy0:iy1,iz0:iz1,1))/fc1
        cmJn(ii,jj,kk,2) = sum(fmJn(ix1,iy0:iy1,iz0:iz1,2))/fc1
        cmJn(ii,jj,kk,3) = sum(fmJn(ix0:ix1,iy0,iz0:iz1,3))/fc1
        cmJn(ii,jj,kk,4) = sum(fmJn(ix0:ix1,iy1,iz0:iz1,4))/fc1
        cmJn(ii,jj,kk,5) = sum(fmJn(ix0:ix1,iy0:iy1,iz0,5))/fc2
        cmJn(ii,jj,kk,6) = sum(fmJn(ix0:ix1,iy0:iy1,iz1,6))/fc2

    end do
    end do
    end do
    !$omp end parallel do
    else
    !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        ! inner surfaces
        cmJ0(ii,jj,kk,1) = sum(fmJ0(ix0,iy0:iy1,iz0:iz1,1))/fc1
        cmJ1(ii,jj,kk,1) = sum(fmJ1(ix0,iy0:iy1,iz0:iz1,1))/fc1

        cmJ0(ii,jj,kk,2) = sum(fmJ0(ix1,iy0:iy1,iz0:iz1,2))/fc1
        cmJ1(ii,jj,kk,2) = sum(fmJ1(ix1,iy0:iy1,iz0:iz1,2))/fc1

        cmJ0(ii,jj,kk,3) = sum(fmJ0(ix0:ix1,iy0,iz0:iz1,3))/fc1
        cmJ1(ii,jj,kk,3) = sum(fmJ1(ix0:ix1,iy0,iz0:iz1,3))/fc1

        cmJ0(ii,jj,kk,4) = sum(fmJ0(ix0:ix1,iy1,iz0:iz1,4))/fc1
        cmJ1(ii,jj,kk,4) = sum(fmJ1(ix0:ix1,iy1,iz0:iz1,4))/fc1

        cmJ0(ii,jj,kk,5) = sum(fmJ0(ix0:ix1,iy0:iy1,iz0,5))/fc2
        cmJ1(ii,jj,kk,5) = sum(fmJ1(ix0:ix1,iy0:iy1,iz0,5))/fc2

        cmJ0(ii,jj,kk,6) = sum(fmJ0(ix0:ix1,iy0:iy1,iz1,6))/fc2
        cmJ1(ii,jj,kk,6) = sum(fmJ1(ix0:ix1,iy0:iy1,iz1,6))/fc2

        ! boundary surfaces
        if ( ii == 1 ) &
        cmJn(ii,jj,kk,1) = sum(fmJn(ix0,iy0:iy1,iz0:iz1,1))/fc1
        if ( ii == ncm(1) ) &
        cmJn(ii,jj,kk,2) = sum(fmJn(ix1,iy0:iy1,iz0:iz1,2))/fc1
        if ( jj == 1 ) &
        cmJn(ii,jj,kk,3) = sum(fmJn(ix0:ix1,iy0,iz0:iz1,3))/fc1
        if ( jj == ncm(2) ) &
        cmJn(ii,jj,kk,4) = sum(fmJn(ix0:ix1,iy1,iz0:iz1,4))/fc1
        if ( kk == 1 ) &
        cmJn(ii,jj,kk,5) = sum(fmJn(ix0:ix1,iy0:iy1,iz0,5))/fc2
        if ( kk == ncm(3) ) &
        cmJn(ii,jj,kk,6) = sum(fmJn(ix0:ix1,iy0:iy1,iz1,6))/fc2
    end do
    end do
    end do
    !$omp end parallel do
    end if

end subroutine

end module

