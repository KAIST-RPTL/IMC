module transient 
    use ace_reactions
    use XS_header
    use evaluate,         only : evalexpr, defparam
    use surface_header,   only : Surfaces, Surfaces_temp, find_surf_idx
    
    
    implicit none
    
    
    ! transient parameters 
    integer :: curr_timestep
    integer :: n_timestep
    real(8) :: del_t            ! [second]
    real(8) :: time_lag
    real(8) :: curr_time
    real(8) :: real_time
    
    real(8) :: w_l = 1.0d-2    ! RR threshold weight
    real(8) :: w_surv = 2.0d-2    ! RR survival weight
    real(8) :: total_weight
    integer :: npg  ! number of precursor groups
    
    ! material / geom change parameters
    logical :: mat_change = .false. 
    logical :: geom_change = .false. 
    integer :: n_interval
    
    type time_purterb
        real(8) :: start_time
        real(8) :: end_time
        integer :: idx1
        integer :: idx2
        character(100) :: fcn
    end type 
    type(time_purterb), allocatable :: purterb(:)
    character(20), allocatable :: move_surf(:) 
    
    
    ! tally 
    real(8), allocatable :: point_power(:,:) 
    
    
    ! source size control 
    integer :: prev_source_size, psize 
    real(8) :: src_ratio, src_ratio0, src_ratio1
    real(8) :: prec_ratio, prec_ratio0, prec_size0
    real(8) :: w_c
    
    integer :: n_microtimestep = 100
    real(8) :: dt_micro
    integer(8) :: n_col, n_col_avg = 0
    integer(8) :: n_cross, n_cross_avg = 0
    real(8) :: SSP  ! source sampling probability 
    
    contains 
    
!===============================================================================
! SET_TRANSIENT - Set variables for transient MC.
!===============================================================================
    subroutine set_transient()
        
        if (.not. do_transient) return 
        
        del_t = time_lag / real(n_timestep,8) 
        total_weight = real(ngen,8)
        real_time = 0 
        
        if (do_pcqs) dt_micro = del_t / real(n_microtimestep, 8) 
        
        if (icore == score) allocate(point_power(n_timestep,n_act))
    end subroutine
    
!===============================================================================
! SAMPLE_PRECURSOR - Sample precursor group & decay constant.
!===============================================================================
    function sample_precursor(iso, E) result (lambda)
        integer, intent(in) :: iso 
        real(8), intent(in) :: E 
        real(8) :: lambda
        integer :: i, pt1, pt2, pt3
        integer :: NE
        real(8) :: r, pdf, ipfac
        real(8) :: prec_group 
        
        
        r = rang(); pdf = 0 
        prec_group = ace(iso)%NXS(8)
        ! sample precursor group
        do i = 1, ace(iso)%NXS(8)-1
            NE = ace(iso) % prcr( i ) % NE
            
            pt1 = 1; pt2 =  NE
            BS: do 
                if(pt2-pt1 == 1) exit BS
                pt3 = (pt2+pt1)/2 
                if (E >= ace(iso)%prcr(i)%E(pt3)) then 
                    pt1 = pt3 
                else 
                    pt2 = pt3
                endif 
            enddo BS
            
            ipfac = max(0.d0, min(1.d0,(E-ace(iso)%prcr(i)%E(pt1))/(ace(iso)%prcr(i)%E(pt1+1)-ace(iso)%prcr(i)%E(pt1))))
            pdf = pdf + ace(iso)%prcr(i)%F(pt1) + ipfac*(ace(iso)%prcr(i)%F(pt1+1)-ace(iso)%prcr(i)%F(pt1))
            
            if (r < pdf) then 
                prec_group = i
                exit 
            endif
        enddo 
        lambda = ace(iso) % prcr( prec_group ) % decay_const
    end function

        
    
!===============================================================================
! COPY_SRC - Copy fission sources into prompt sources.
!===============================================================================
    subroutine copy_src (src_bank)
        type(bank),intent(in) :: src_bank(:)
        integer :: i
        integer :: size_bank
        
        if (allocated(prompt_bank)) deallocate(prompt_bank) 
        
        size_bank = size(src_bank)
        allocate(prompt_bank(1:size_bank))
        
        
        prompt_bank(:) = src_bank(:)
        prompt_bank(:)%delayed = .false. 
        
        !do i = 1, size_bank
        !    prompt_bank(i)%wgt     = src_bank(i)%wgt
        !    prompt_bank(i)%xyz(:)  = src_bank(i)%xyz(:)
        !    prompt_bank(i)%uvw(:)  = src_bank(i)%uvw(:)
        !    prompt_bank(i)%E       = src_bank(i)%E
        !    prompt_bank(i)%G       = src_bank(i)%G
        !    prompt_bank(i)%delayed = .false. 
        !    prompt_bank(i)%time    = src_bank(i)%time
        !enddo 
        
    end subroutine 
    
    
    
!===============================================================================
! FINALIZE_SRC - Deallocate prompt & delayed source banks for next simulation.
!===============================================================================
    subroutine finalize_src() 
        if(allocated(delayed_bank)) deallocate(delayed_bank) 
        if(allocated(prompt_bank))  deallocate(prompt_bank) 
        if(allocated(dynamic_bank)) deallocate(dynamic_bank) 
        if(allocated(source_bank))  deallocate(source_bank) 
        if(allocated(split_bank))  deallocate(split_bank) 
        if(allocated(split_bank_temp))  deallocate(split_bank_temp) 
        if(allocated(prec_bank))  deallocate(prec_bank) 
    end subroutine
    
    
    subroutine normalizeInitialSource() 
        use mpi 
        implicit none
        integer :: i, j, n
        integer :: isize, ierr
        real(8) :: w, norm
        
        
        ! Normalize prompt_bank & prec_bank in server core
        !if (icore == score) then 
            n = size(prec_bank) 
            w = sum(prec_bank(:)%wgt) 
            norm = ngen / w
            
            !if (icore==score) print *, 'normalizeInitialSource  ', size(prompt_bank) 
            
            !if (icore==score) print *, sum(prompt_bank(:)%wgt)/size(prompt_bank), size(prompt_bank)
            !if (icore==score) print *, sum(prec_bank(:)%wgt)/size(prec_bank), size(prec_bank)
            
            prompt_bank(:)%wgt = prompt_bank(:)%wgt * norm
            prec_bank(:)%wgt   = prec_bank(:)%wgt * norm
            
            !if (icore==score) print *, sum(prompt_bank(:)%wgt)/size(prompt_bank)
            !if (icore==score) print *, sum(prec_bank(:)%wgt)/size(prec_bank)
            !if (icore==score) print *, norm
            !stop 
            
        !endif 
        
    end subroutine 
    
    
    ! ==================================================================================
    !                     MG_XS adjust subroutines for C5G7-TD Problem 
    ! ==================================================================================
    
    subroutine adjust_MG_XS() 
        integer :: i 
        integer :: idx
        real(8) :: t 
        
        if (E_mode == 1) return 
        if (.not. do_transient) return 
        
        t = (curr_timestep-1) * del_t
        
        XS_MG(6) = XS_MG_temp(6) 
        MGD(6)   = MGD_temp(6) 
        
        ! ===================================================================================================
        !     C5G7-TD0-5 
        ! ===================================================================================================
        ! CR : 5 / GT : 6
        ! 0.0s ~ 0.5s steady 
        ! 0.5s ~ 1.5s step 1 
        ! 1.5s ~ 2.5s step 2 
        ! 2.5s ~ 3.0s step 3 
        
        if (t<0.5) then  
            XS_MG(6) = XS_MG_temp(6) 
            MGD(6)   = MGD_temp(6) 
        elseif (t>=0.5 .and. t<1.5) then  
            XS_MG(6)%sig_tr(:)     = XS_MG_temp(6)%sig_tr(:)     + 0.1 * (XS_MG_temp(5)%sig_tr(:)     - XS_MG_temp(6)%sig_tr(:)    )
            XS_MG(6)%sig_abs(:)    = XS_MG_temp(6)%sig_abs(:)    + 0.1 * (XS_MG_temp(5)%sig_abs(:)    - XS_MG_temp(6)%sig_abs(:)   )
            XS_MG(6)%sig_cap(:)    = XS_MG_temp(6)%sig_cap(:)    + 0.1 * (XS_MG_temp(5)%sig_cap(:)    - XS_MG_temp(6)%sig_cap(:)   )
            XS_MG(6)%sig_fis(:)    = XS_MG_temp(6)%sig_fis(:)    + 0.1 * (XS_MG_temp(5)%sig_fis(:)    - XS_MG_temp(6)%sig_fis(:)   )
            XS_MG(6)%nu(:)         = XS_MG_temp(6)%nu(:)         + 0.1 * (XS_MG_temp(5)%nu(:)         - XS_MG_temp(6)%nu(:)        )
            XS_MG(6)%chi(:)        = XS_MG_temp(6)%chi(:)        + 0.1 * (XS_MG_temp(5)%chi(:)        - XS_MG_temp(6)%chi(:)       )
            XS_MG(6)%sig_scat(:,:) = XS_MG_temp(6)%sig_scat(:,:) + 0.1 * (XS_MG_temp(5)%sig_scat(:,:) - XS_MG_temp(6)%sig_scat(:,:))
            
            MGD(6)%beta(:)      = MGD_temp(6)%beta(:)      + 0.1 * (MGD_temp(5)%beta(:)      - MGD_temp(6)%beta(:)     )
            MGD(6)%lambda(:)    = MGD_temp(6)%lambda(:)    + 0.1 * (MGD_temp(5)%lambda(:)    - MGD_temp(6)%lambda(:)   )
            MGD(6)%vel(:)       = MGD_temp(6)%vel(:)       + 0.1 * (MGD_temp(5)%vel(:)       - MGD_temp(6)%vel(:)      )
            MGD(6)%spectra(:,:) = MGD_temp(6)%spectra(:,:) + 0.1 * (MGD_temp(5)%spectra(:,:) - MGD_temp(6)%spectra(:,:))
        elseif (t>=1.5 .and. t<2.5) then 
            XS_MG(6)%sig_tr(:)     = XS_MG_temp(6)%sig_tr(:)     + 0.05 * (XS_MG_temp(5)%sig_tr(:)     - XS_MG_temp(6)%sig_tr(:)    )
            XS_MG(6)%sig_abs(:)    = XS_MG_temp(6)%sig_abs(:)    + 0.05 * (XS_MG_temp(5)%sig_abs(:)    - XS_MG_temp(6)%sig_abs(:)   )
            XS_MG(6)%sig_cap(:)    = XS_MG_temp(6)%sig_cap(:)    + 0.05 * (XS_MG_temp(5)%sig_cap(:)    - XS_MG_temp(6)%sig_cap(:)   )
            XS_MG(6)%sig_fis(:)    = XS_MG_temp(6)%sig_fis(:)    + 0.05 * (XS_MG_temp(5)%sig_fis(:)    - XS_MG_temp(6)%sig_fis(:)   )
            XS_MG(6)%nu(:)         = XS_MG_temp(6)%nu(:)         + 0.05 * (XS_MG_temp(5)%nu(:)         - XS_MG_temp(6)%nu(:)        )
            XS_MG(6)%chi(:)        = XS_MG_temp(6)%chi(:)        + 0.05 * (XS_MG_temp(5)%chi(:)        - XS_MG_temp(6)%chi(:)       )
            XS_MG(6)%sig_scat(:,:) = XS_MG_temp(6)%sig_scat(:,:) + 0.05 * (XS_MG_temp(5)%sig_scat(:,:) - XS_MG_temp(6)%sig_scat(:,:))
            
            MGD(6)%beta(:)      = MGD_temp(6)%beta(:)      + 0.05 * (MGD_temp(5)%beta(:)      - MGD_temp(6)%beta(:)     )
            MGD(6)%lambda(:)    = MGD_temp(6)%lambda(:)    + 0.05 * (MGD_temp(5)%lambda(:)    - MGD_temp(6)%lambda(:)   )
            MGD(6)%vel(:)       = MGD_temp(6)%vel(:)       + 0.05 * (MGD_temp(5)%vel(:)       - MGD_temp(6)%vel(:)      )
            MGD(6)%spectra(:,:) = MGD_temp(6)%spectra(:,:) + 0.05 * (MGD_temp(5)%spectra(:,:) - MGD_temp(6)%spectra(:,:))
        else
            XS_MG(6) = XS_MG_temp(6) 
            MGD(6)   = MGD_temp(6) 
        endif 
        
        ! ===================================================================================================
        !     C5G7-TD1-5 
        ! ===================================================================================================
        !if (t<0.5) then  
        !    XS_MG(6) = XS_MG_temp(6) 
        !    MGD(6)   = MGD_temp(6) 
        !elseif (t>=0.5 .and. t<1.5) then  
        !    XS_MG(6)%sig_tr(:)     = XS_MG_temp(6)%sig_tr(:)     + (t-0.5)*0.01 * (XS_MG_temp(5)%sig_tr(:)     - XS_MG_temp(6)%sig_tr(:)    )
        !    XS_MG(6)%sig_abs(:)    = XS_MG_temp(6)%sig_abs(:)    + (t-0.5)*0.01 * (XS_MG_temp(5)%sig_abs(:)    - XS_MG_temp(6)%sig_abs(:)   )
        !    XS_MG(6)%sig_cap(:)    = XS_MG_temp(6)%sig_cap(:)    + (t-0.5)*0.01 * (XS_MG_temp(5)%sig_cap(:)    - XS_MG_temp(6)%sig_cap(:)   )
        !    XS_MG(6)%sig_fis(:)    = XS_MG_temp(6)%sig_fis(:)    + (t-0.5)*0.01 * (XS_MG_temp(5)%sig_fis(:)    - XS_MG_temp(6)%sig_fis(:)   )
        !    XS_MG(6)%nu(:)         = XS_MG_temp(6)%nu(:)         + (t-0.5)*0.01 * (XS_MG_temp(5)%nu(:)         - XS_MG_temp(6)%nu(:)        )
        !    XS_MG(6)%chi(:)        = XS_MG_temp(6)%chi(:)        + (t-0.5)*0.01 * (XS_MG_temp(5)%chi(:)        - XS_MG_temp(6)%chi(:)       )
        !    XS_MG(6)%sig_scat(:,:) = XS_MG_temp(6)%sig_scat(:,:) + (t-0.5)*0.01 * (XS_MG_temp(5)%sig_scat(:,:) - XS_MG_temp(6)%sig_scat(:,:))
        !    
        !elseif (t>=1.5 .and. t<2.5) then 
        !    XS_MG(6)%sig_tr(:)     = XS_MG_temp(6)%sig_tr(:)     + (2.5-t)*0.01 * (XS_MG_temp(5)%sig_tr(:)     - XS_MG_temp(6)%sig_tr(:)    )
        !    XS_MG(6)%sig_abs(:)    = XS_MG_temp(6)%sig_abs(:)    + (2.5-t)*0.01 * (XS_MG_temp(5)%sig_abs(:)    - XS_MG_temp(6)%sig_abs(:)   )
        !    XS_MG(6)%sig_cap(:)    = XS_MG_temp(6)%sig_cap(:)    + (2.5-t)*0.01 * (XS_MG_temp(5)%sig_cap(:)    - XS_MG_temp(6)%sig_cap(:)   )
        !    XS_MG(6)%sig_fis(:)    = XS_MG_temp(6)%sig_fis(:)    + (2.5-t)*0.01 * (XS_MG_temp(5)%sig_fis(:)    - XS_MG_temp(6)%sig_fis(:)   )
        !    XS_MG(6)%nu(:)         = XS_MG_temp(6)%nu(:)         + (2.5-t)*0.01 * (XS_MG_temp(5)%nu(:)         - XS_MG_temp(6)%nu(:)        )
        !    XS_MG(6)%chi(:)        = XS_MG_temp(6)%chi(:)        + (2.5-t)*0.01 * (XS_MG_temp(5)%chi(:)        - XS_MG_temp(6)%chi(:)       )
        !    XS_MG(6)%sig_scat(:,:) = XS_MG_temp(6)%sig_scat(:,:) + (2.5-t)*0.01 * (XS_MG_temp(5)%sig_scat(:,:) - XS_MG_temp(6)%sig_scat(:,:))
        !    
        !else
        !    XS_MG(6) = XS_MG_temp(6) 
        !    MGD(6)   = MGD_temp(6) 
        !endif         
        
        
        
        ! ===================================================================================================
        !     C5G7-TD2-5 
        ! ===================================================================================================
        !if (t<0.5) then  
        !    XS_MG(6) = XS_MG_temp(6) 
        !    MGD(6)   = MGD_temp(6) 
        !elseif (t>=0.5 .and. t<1.5) then  
        !    XS_MG(6)%sig_tr(:)     = XS_MG_temp(6)%sig_tr(:)     + (t-0.5)*0.1 * (XS_MG_temp(5)%sig_tr(:)     - XS_MG_temp(6)%sig_tr(:)    )
        !    XS_MG(6)%sig_abs(:)    = XS_MG_temp(6)%sig_abs(:)    + (t-0.5)*0.1 * (XS_MG_temp(5)%sig_abs(:)    - XS_MG_temp(6)%sig_abs(:)   )
        !    XS_MG(6)%sig_cap(:)    = XS_MG_temp(6)%sig_cap(:)    + (t-0.5)*0.1 * (XS_MG_temp(5)%sig_cap(:)    - XS_MG_temp(6)%sig_cap(:)   )
        !    XS_MG(6)%sig_fis(:)    = XS_MG_temp(6)%sig_fis(:)    + (t-0.5)*0.1 * (XS_MG_temp(5)%sig_fis(:)    - XS_MG_temp(6)%sig_fis(:)   )
        !    XS_MG(6)%nu(:)         = XS_MG_temp(6)%nu(:)         + (t-0.5)*0.1 * (XS_MG_temp(5)%nu(:)         - XS_MG_temp(6)%nu(:)        )
        !    XS_MG(6)%chi(:)        = XS_MG_temp(6)%chi(:)        + (t-0.5)*0.1 * (XS_MG_temp(5)%chi(:)        - XS_MG_temp(6)%chi(:)       )
        !    XS_MG(6)%sig_scat(:,:) = XS_MG_temp(6)%sig_scat(:,:) + (t-0.5)*0.1 * (XS_MG_temp(5)%sig_scat(:,:) - XS_MG_temp(6)%sig_scat(:,:))
        !    
        !elseif (t>=1.5 .and. t<2.5) then 
        !    XS_MG(6)%sig_tr(:)     = XS_MG_temp(6)%sig_tr(:)     + (2.5-t)*0.1 * (XS_MG_temp(5)%sig_tr(:)     - XS_MG_temp(6)%sig_tr(:)    )
        !    XS_MG(6)%sig_abs(:)    = XS_MG_temp(6)%sig_abs(:)    + (2.5-t)*0.1 * (XS_MG_temp(5)%sig_abs(:)    - XS_MG_temp(6)%sig_abs(:)   )
        !    XS_MG(6)%sig_cap(:)    = XS_MG_temp(6)%sig_cap(:)    + (2.5-t)*0.1 * (XS_MG_temp(5)%sig_cap(:)    - XS_MG_temp(6)%sig_cap(:)   )
        !    XS_MG(6)%sig_fis(:)    = XS_MG_temp(6)%sig_fis(:)    + (2.5-t)*0.1 * (XS_MG_temp(5)%sig_fis(:)    - XS_MG_temp(6)%sig_fis(:)   )
        !    XS_MG(6)%nu(:)         = XS_MG_temp(6)%nu(:)         + (2.5-t)*0.1 * (XS_MG_temp(5)%nu(:)         - XS_MG_temp(6)%nu(:)        )
        !    XS_MG(6)%chi(:)        = XS_MG_temp(6)%chi(:)        + (2.5-t)*0.1 * (XS_MG_temp(5)%chi(:)        - XS_MG_temp(6)%chi(:)       )
        !    XS_MG(6)%sig_scat(:,:) = XS_MG_temp(6)%sig_scat(:,:) + (2.5-t)*0.1 * (XS_MG_temp(5)%sig_scat(:,:) - XS_MG_temp(6)%sig_scat(:,:))
        !    
        !else
        !    XS_MG(6) = XS_MG_temp(6) 
        !    MGD(6)   = MGD_temp(6) 
        !endif                 
        
    end subroutine
    
    subroutine save_MG_XS() 
        integer :: xs_size
        
        if (E_mode == 1) return 
        if (.not. do_transient) return 
        
        xs_size = size(XS_MG) 
        allocate(XS_MG_temp(1:xs_size))
        XS_MG_temp(:) = XS_MG(:)
        
        xs_size = size(MGD) 
        allocate(MGD_temp(1:xs_size))
        MGD_temp(:) = MGD(:)


    end subroutine
    
    
    subroutine restore_MG_XS() 
    
        if (E_mode == 1) return 
        if (.not. do_transient) return 
        
        deallocate(XS_MG) 
        call move_alloc(XS_MG_temp, XS_MG) 
        
        deallocate(MGD) 
        call move_alloc(MGD_temp, MGD) 
    end subroutine
    
    
    subroutine adjust_CE_MAT() 
        real(8) :: t 
        
        
        if (E_mode == 0) return 
        if (.not. do_transient) return 
        
        t = (curr_timestep-1) * del_t
        
        if (t > 5 .and. t < 15) then
            materials(1)%numden(1) = 0.045E+24
        else 
            materials(1)%numden(1) = 0.044744E+24
        endif
        
    end subroutine 
    
    subroutine scramble( array )
        implicit none 
        integer,intent(inout) :: array(:)
        integer :: i, j, k, n, m, itemp, number_of_values
        real(8) :: u 
        
        number_of_values = size(array)
        array=[(i,i=1,number_of_values)]
        n=1; m=number_of_values
        do k=1,2
            do i=1,m
                call random_number(u)
                j = n + FLOOR((m+1-n)*u)
                itemp=array(j); array(j)=array(i); array(i)=itemp
            enddo
        enddo
    end subroutine scramble
    
    
    
    subroutine condition_change()
        implicit none
        real(8) :: time
        integer :: i, idx_time, idx_surf, idx_parm, idx_mat1, idx_mat2
        real(8) :: newval 
        integer :: i_g, j_g
        character(20) :: surf
        
        time = curr_time
        
        if (geom_change) then 
            surfaces = surfaces_temp
        elseif (mat_change) then 
            if (E_mode == 0) then 
                XS_MG = XS_MG_temp
                MGD = MGD_temp
            else 
                materials = materials_temp
            endif 
        endif 
        
        idx_time = 0 
        do i = 1, n_interval
            if (purterb(i)%start_time <= time  .and. purterb(i)%end_time > time) then 
                idx_time = i
                
                if (geom_change) then 
                    !idx_surf = purterb(i)%idx1
                    surf = surfaces(purterb(i)%idx1)%surf_id
                    idx_surf = find_surf_idx(surfaces, adjustl(surf))
                    idx_parm = purterb(i)%idx2
                    call defparam('t',time)
                    call evalexpr(purterb(i)%fcn, newval)
                    surfaces(idx_surf)%parmtrs(idx_parm) = newval
                    !print *, adjustl(surf), newval 
                
                elseif (mat_change) then 
                        
                    idx_mat1 = purterb(i)%idx1
                    idx_mat2 = purterb(i)%idx2

                    if (E_mode == 0) then 
                        do i_g = 1, n_group
                            call defparam('a',XS_MG_temp(idx_mat1)%sig_tr(i_g))
                            call defparam('b',XS_MG_temp(idx_mat2)%sig_tr(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            XS_MG(idx_mat1)%sig_tr(i_g) = newval
                            
                            
                            call defparam('a',XS_MG_temp(idx_mat1)%sig_abs(i_g))
                            call defparam('b',XS_MG_temp(idx_mat2)%sig_abs(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            XS_MG(idx_mat1)%sig_abs(i_g) = newval
                            
                            call defparam('a',XS_MG_temp(idx_mat1)%sig_cap(i_g))
                            call defparam('b',XS_MG_temp(idx_mat2)%sig_cap(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            XS_MG(idx_mat1)%sig_cap(i_g) = newval
                            
                            call defparam('a',XS_MG_temp(idx_mat1)%sig_fis(i_g))
                            call defparam('b',XS_MG_temp(idx_mat2)%sig_fis(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            XS_MG(idx_mat1)%sig_fis(i_g) = newval
                            
                            call defparam('a',XS_MG_temp(idx_mat1)%chi(i_g))
                            call defparam('b',XS_MG_temp(idx_mat2)%chi(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            XS_MG(idx_mat1)%chi(i_g) = newval
                            
                            call defparam('a',XS_MG_temp(idx_mat1)%nu(i_g))
                            call defparam('b',XS_MG_temp(idx_mat2)%nu(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            XS_MG(idx_mat1)%nu(i_g) = newval
                            
                            do j_g = 1, n_group 
                                call defparam('a',XS_MG_temp(idx_mat1)%sig_scat(i_g, j_g))
                                call defparam('b',XS_MG_temp(idx_mat2)%sig_scat(i_g, j_g))
                                call defparam('t',time)
                                call evalexpr(purterb(i)%fcn, newval)
                                XS_MG(idx_mat1)%sig_scat(i_g, j_g) = newval
                            enddo 
                            
                            
                        enddo 
                        
                        do i_g = 1, n_group 
                            call defparam('a',MGD_temp(idx_mat1)%vel(i_g))
                            call defparam('b',MGD_temp(idx_mat2)%vel(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            MGD(idx_mat1)%vel(i_g) = newval
                        enddo 
                        
                        do i_g = 1, npg 
                            call defparam('a',MGD_temp(idx_mat1)%beta(i_g))
                            call defparam('b',MGD_temp(idx_mat2)%beta(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            MGD(idx_mat1)%beta(i_g) = newval

                            call defparam('a',MGD_temp(idx_mat1)%lambda(i_g))
                            call defparam('b',MGD_temp(idx_mat2)%lambda(i_g))
                            call defparam('t',time)
                            call evalexpr(purterb(i)%fcn, newval)
                            MGD(idx_mat1)%lambda(i_g) = newval

                            do j_g = 1, n_group 
                                call defparam('a',MGD_temp(idx_mat1)%spectra(i_g, j_g))
                                call defparam('b',MGD_temp(idx_mat2)%spectra(i_g, j_g))
                                call defparam('t',time)
                                call evalexpr(purterb(i)%fcn, newval)
                                MGD(idx_mat1)%spectra(i_g, j_g) = newval
                            enddo 
                        enddo 
                        
                    else 
                        call defparam('x0', materials_temp(idx_mat1)%numden(idx_mat2))
                        call defparam('t', time)
                        call evalexpr(purterb(i)%fcn, newval)
                        materials(idx_mat1)%numden(idx_mat2) = newval
                    endif 
                    
                endif
                
                
                
                
            endif
        enddo 
        
        
    end subroutine 
    
    
    subroutine save_condition() 
        implicit none 
        integer :: size_data
        
        if (geom_change) then 
            size_data = size(surfaces) 
            if (.not. allocated(Surfaces_temp)) allocate(Surfaces_temp(1:size_data)) 
            surfaces_temp(:) = surfaces(:)
            
        elseif (mat_change) then 
            if (E_mode == 0) then 
                size_data = size(XS_MG)
                if (.not. allocated(XS_MG_temp))  allocate(XS_MG_temp(1:size_data)) 
                XS_MG_temp(:) = XS_MG(:)
                
                size_data = size(MGD) 
                if (.not. allocated(MGD_temp)) allocate(MGD_temp(1:size_data))
                MGD_temp(:) = MGD(:)
                
            else 
                size_data = size(materials)
                if (.not. allocated(materials_temp)) allocate(materials_temp(1:size_data)) 
                materials_temp(:) = materials(:)
            endif 
        endif 
        
    end subroutine 
    
    
    subroutine restore_condition() 
        implicit none 
        
        if (geom_change) then 
            deallocate(surfaces)
            call move_alloc(surfaces_temp, surfaces)
        elseif (mat_change) then 
            if (E_mode == 0) then 
                deallocate(XS_MG)
                call move_alloc(XS_MG_temp, XS_MG)
                deallocate(MGD)
                call move_alloc(MGD_temp, MGD)
            else 
                deallocate(materials)
                call move_alloc(materials_temp, materials)
            endif 
        endif 
        
    end subroutine 
    
end module


