    Module calc_tb_multiview_AMSU

        use rtm_tables_AMSU, only : find_abs_q_AMSU, find_cld_abs_AMSU
        use rtm_tables_AMSU, only : find_emiss_AMSU, find_emiss_sea_ice_AMSU
        use rtm, only : atm_tran_p
        use AMSU_Constants
        use, intrinsic :: ieee_arithmetic, Only : ieee_is_nan


        integer(4),parameter        :: OCEAN = 1
        integer(4),parameter        :: LAND  = 2
        integer(4),parameter        :: SEA_ICE = 3
        
        real(8),parameter,private        ::     R_GAS    =  8.3145112d0                   !6.0221367d23 * 1.380658d-23   ;J/mol/K
        real(8),parameter,private        ::     M_W_AIR  =  2.8966d-2                     !kg/mol
        real(8),parameter,private        ::     M_W_H2O  =  1.8015324d-2                  !kg/mol
    
        real(8),parameter,private        ::     c_vap    =  M_W_H2O/M_W_AIR
        real(8),parameter,private        ::     c_vap2   =  1.0 - c_vap
    
        real(8),parameter,private        ::     c_air    =  R_GAS/M_W_AIR
        real(8),parameter,private        ::     c_h2o    =  R_GAS/M_W_H2O
   
        real(8),parameter,private        ::     g = 9.80665  !m/s^2

    contains

        subroutine calc_tb_multiview_table_AMSU_multi_profiles(  &
                                                num_profiles,    &  ! number of profiles
                                                num_levels,      &  ! number of levels
                                                num_views,       &  ! number of views
                                                t,               &  ! temperature
                                                p,               &    ! pressure
                                                q,               &    ! specific humidity
                                                cld,             &  ! cloud water
                                                cld_mix,         &  ! flag for how to intepret cld -- if gt 0 then mixing ratio
                                                t_surf,          &  ! surface temperature
                                                p_surf,          &  ! surface pressure
                                                q_surf,          &  ! surface specific humidity
                                                cld_surf,        &  ! surface cloud water
                                                amsu_channel,    &  ! amsu channel number (1-4)
                                                theta,           &  ! EIA angle (degrees)  array(1:num_views)
                                                wind,            &  ! windspeed
                                                land_emiss_in,   &  ! land_emissivity input....
                                                emissivity,      &    ! emissivity for each pol, surface, view
                                                surf_wt,         &  ! weight from surface emission for each surface,view
                                                space_wt,        &  ! weight from space
                                                tb,              &    ! brightness temperature at TOA for each surface,view
                                                tb_up,           &  ! upwelling brightness temperature for each view
                                                tb_dw,           &  ! downwelling brightness temperature for each view
                                                tau,             &  ! total transmittance for each view
                                                error)

            ! This routine interates through each profile and calls calc_tb_multiview_table_AMSU for each profile. 
            ! It also adds the surface as a level and adds a level above the top of the atmosphere, and checks for 
            ! NaN values and levels with pressure greater than surface pressure.
        
            implicit none
             
            integer(4),intent(in)                                           :: num_profiles
            integer(4),intent(in)                                           :: num_levels 
            integer(4),intent(in)                                           :: num_views
            real(4),dimension(num_levels,num_profiles),intent(in)           :: t        
            real(4),dimension(num_levels),intent(in)                        :: p        
            real(4),dimension(num_levels,num_profiles),intent(in)           :: q
            real(4),dimension(num_levels,num_profiles),intent(in)           :: cld
            integer(4),intent(in)                                           :: cld_mix
            real(4),dimension(num_profiles),intent(in)                      :: t_surf
            real(4),dimension(num_profiles),intent(in)                      :: p_surf
            real(4),dimension(num_profiles),intent(in)                      :: q_surf
            real(4),dimension(num_profiles),intent(in)                      :: cld_surf
            integer(4),intent(in)                                           :: amsu_channel
            real(4),dimension(num_views),intent(in)                         :: theta
            
            real(4),dimension(num_profiles),intent(in)                      :: wind
            real(4),dimension(num_views),intent(in)                         :: land_emiss_in

            real(4),dimension(3,num_views,num_profiles),intent(out)         :: emissivity 
            real(4),dimension(3,num_views,num_profiles),intent(out)         :: surf_wt        
            real(4),dimension(3,num_views,num_profiles),intent(out)         :: space_wt

            real(4),dimension(3,num_views,num_profiles),intent(out)         :: tb
            
            real(4),dimension(num_views,num_profiles),intent(out)           :: tb_up
            real(4),dimension(num_views,num_profiles),intent(out)           :: tb_dw
            real(4),dimension(num_views,num_profiles),intent(out)           :: tau            
                
            integer(4),intent(out)                                          :: error
            integer(4)                             :: i

            real(4),dimension(0:num_levels)      :: t_i
            real(4),dimension(0:num_levels)      :: p_i
            real(4),dimension(0:num_levels)      :: q_i
            real(4),dimension(0:num_levels)      :: cld_i
            integer(4)                           :: num_levels_ok
            integer(4)                           :: level

            error = 0
            
            do i = 1,num_profiles
                 
                num_levels_ok = 0

                do level = 1,num_levels
                    if (p(level) .ge. p_surf(i)) cycle
                    if (ieee_is_nan(t(level,i))) cycle
                    if (ieee_is_nan(q(level,i))) cycle
                    
                    t_i(num_levels_ok) = t(level,i)
                    p_i(num_levels_ok) = p(level)
                    q_i(num_levels_ok) = q(level,i)
                    cld_i(num_levels_ok) = cld(level,i)
                    if (ieee_is_nan(cld_i(num_levels_ok))) cld_i(num_levels_ok) = 0.0
                    num_levels_ok = num_levels_ok + 1
                end do

                t_i(num_levels_ok) = t_surf(i)
                p_i(num_levels_ok) = p_surf(i)
                q_i(num_levels_ok) = q_surf(i)
                cld_i(num_levels_ok) = cld_surf(i)
                
                call calc_tb_multiview_table_AMSU(t_i(0:num_levels_ok),p_i(0:num_levels_ok),q_i(0:num_levels_ok),cld_i(0:num_levels_ok), &
                     cld_mix,num_levels_ok,amsu_channel,theta,num_views,wind(i),land_emiss_in, &
                     emissivity(:,:,i),surf_wt(:,:,i),space_wt(:,:,i),tb(:,:,i),tb_up(:,i),tb_dw(:,i),tau(:,i),error)
            enddo
            end subroutine calc_tb_multiview_table_AMSU_multi_profiles
    
        subroutine calc_tb_multiview_table_AMSU(t,               &  ! temperature
                                                p,               &    ! pressure
                                                q,               &    ! specific humidity
                                                cld,             &  ! cloud water
                                                cld_mix,         &  ! flag for how to intepret cld -- if gt 0 then mixing ratio
                                                num_levels,      &  ! number of levels
                                                amsu_channel,    &  ! amsu channel number (1-4)
                                                theta,           &  ! EIA angle (degrees)  array(1:num_views)
                                                num_views,       &  ! number of views
                                                wind,            &  ! windspeed
                                                land_emiss_in,   &  ! land_emissivity input....
                                                emissivity,      &    ! emissivity for each pol, surface, view
                                                surf_wt,         &  ! weight from surface emission for each surface,view
                                                space_wt,        &  ! weight from space
                                                tb,              &    ! brightness temperature at TOA for each surface,view
                                                tb_up,           &  ! upwelling brightness temperature for each view
                                                tb_dw,           &  ! downwelling brightness temperature for each view
                                                tau,             &  ! total transmittance for each view
                                                error)

            implicit none
             
            integer(4)                           :: num_levels 
            real(4),dimension(0:num_levels)      :: t        
            real(4),dimension(0:num_levels)      :: p        
            real(4),dimension(0:num_levels)      :: q
            real(4),dimension(0:num_levels)      :: cld
            integer(4)                           :: cld_mix
            integer(4)                           :: amsu_channel
            integer(4)                           :: num_views
            integer(4)                           :: view_num
            real(4),dimension(num_views)         :: theta
            real(4),dimension(num_views)         :: theta_view
            real(4)                              :: wind
            real(4),dimension(3,num_views)       :: surf_wt        
            real(4),dimension(3,num_views)       :: space_wt

            real(4),dimension(3,num_views)       :: tb
            
            real(4),dimension(num_views)         :: tb_up
            real(4),dimension(num_views)         :: tb_dw
            real(4),dimension(num_views)         :: tau            
                
            real(4),dimension(0:num_levels)      :: transmittance_up   ! 
            real(4),dimension(0:num_levels)      :: transmittance_down ! from bottom of each level
            real(4)                              :: tbup
            real(4)                              :: tbdw
            integer(4)                           :: level
            integer(4)                           :: error

            real(4),dimension(3,num_views)       :: emissivity
            real(4),dimension(num_views)         :: land_emiss_in
            
            real(4)                              :: rho_dry
            real(4)                              :: rho_vap
            real(4)                              :: rho_cld
            real(4)                              :: pv
            real(4)                              :: p_dry
            real(4)                              :: cld_level
            real(4)                              :: cld_abs_np_m
            real(4)                              :: trans

            real(4),dimension(0:NUM_LEVELS)      :: cloud_dens
            real(4),dimension(0:NUM_LEVELS)      :: scl_interface_pressure
            real(4),dimension(0:NUM_LEVELS)      :: total_abs
            real(4),dimension(0:NUM_LEVELS)      :: dry_abs
            real(4),dimension(0:NUM_LEVELS)      :: cloud_abs
            real(4),dimension(0:NUM_LEVELS)      :: rho
            

            integer(4)                           :: surf_type
            real(4)                              :: sst


            ! find the absorbtion coefficients for each level
            ! these are independent of angle, so they only need to be calculated 
            ! once for each profile/frequency

            if (maxval(p) .gt. 2000) then
                scl_interface_pressure = p/100.0
            else
                scl_interface_pressure = p
            endif

            ! calculate absorption coefficients at each level
            
            ! absorption coefficients are returned in nepers per Km
            cloud_abs = 0.0
            total_abs = 0.0
            do level = 0,num_levels
                   dry_abs(level) = find_abs_q_AMSU(t(level),scl_interface_pressure(level),q(level),error)
     
                   ! find partial pressure of air and h2o
                   pv = p(level)*q(level)/(c_vap+c_vap2*q(level))
                   p_dry = p(level) - pv
      
                   ! find densities
                   rho_dry = 100.0*p_dry/(c_air*T(level))  ! the 100 converts pressures into Pa
                   rho_vap = 100.0*pv/(c_h2o*T(level))
                   rho(level) = rho_dry+rho_vap
        
                   if (cld(level) .gt. 0.0) then  ! cld(lev) is *liquid* cloud water in g/cm^3
                     cld_level = cld(level)
                     if (cld_mix .gt. 0) then 
                        cld_level = cld_level*rho(level)
                     endif
                     cloud_dens(level) = cld_level  ! this is in kg/m^3
                     cld_abs_np_m   = find_cld_abs_AMSU(t(level),cld_level,error)  ! nepers/m
                     cloud_abs(level) = 100.0*cld_Abs_np_m/(rho(level) * g) ! nepers/hPa
                   endif
                   total_abs(level) = dry_abs(level) + cloud_abs(level)

            enddo
            
            !do view_num = 1,num_views 
            do view_num = 1,num_views     
                call atm_tran_p(    num_levels,                &
                                    theta(view_num),        &
                                    t,                        &
                                    scl_interface_pressure, &
                                    total_Abs,                &
                                    transmittance_up,        &
                                    transmittance_down,        &
                                    tbdw,                    &
                                    tbup)

                ! calculate surface emissivities
                emissivity(OCEAN,view_num) =  find_emiss_AMSU(t(num_levels),wind,view_num,error)
                emissivity(LAND,view_num)  =  land_emiss_in(view_num)   
                emissivity(SEA_ICE,view_num)= find_emiss_sea_ice_AMSU(T(num_levels),view_num,error)
                trans = transmittance_down(0)
 
                tb(OCEAN,view_num) =   tbup + trans * ((tbdw*(1.0-emissivity(OCEAN,view_num)))+  (emissivity(OCEAN,view_num)*  t(num_levels))) 
                tb(LAND,view_num) =    tbup + trans * ((tbdw*(1.0-emissivity(LAND,view_num)))+   (emissivity(LAND,view_num)*   t(num_levels)))
                tb(SEA_ICE,view_num) = tbup + trans * ((tbdw*(1.0-emissivity(SEA_ICE,view_num)))+  (emissivity(SEA_ICE,view_num)*  t(num_levels)))
                 
                surf_wt(OCEAN,view_num) = trans * emissivity(OCEAN,view_num)
                surf_wt(LAND, view_num) = trans * emissivity(LAND ,view_num)
                surf_wt(SEA_ICE, view_num) = trans * emissivity(SEA_ICE ,view_num)

                space_wt(OCEAN,view_num) = (1.0-emissivity(OCEAN,view_num)) *trans*trans
                space_wt(LAND,view_num) =  (1.0-emissivity(LAND,view_num)) *trans*trans
                space_wt(SEA_ICE,view_num) =  (1.0-emissivity(SEA_ICE,view_num)) *trans*trans
                 
                tau(view_num) = trans
                tb_up(view_num) = tbup
                tb_dw(view_num) = tbdw
            enddo
            error = 0

        end subroutine calc_tb_multiview_table_AMSU
    end module calc_tb_multiview_AMSU