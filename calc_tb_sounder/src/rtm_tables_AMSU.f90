module rtm_tables_AMSU

    use netcdf

    implicit none

    integer(4),parameter,private            :: num_t = 200
    integer(4),parameter,private            :: num_p = 110
    integer(4),parameter,private            :: num_q = 150
    integer(4),parameter,private            :: num_w = 30
    integer(4),parameter,private            :: num_fov = 15

    real(4),dimension(0:num_t,0:num_p,0:num_q)    :: abs_table_q
    real(4),dimension(0:num_t)                    :: cld_abs_table
    real(4),dimension(0:num_t,0:num_w,num_fov)    :: ocean_emiss_table
    real(4),dimension(0:num_t,num_fov)          :: sea_ice_emiss_table
    real(4)                                        :: T0,delta_t
    real(4)                                     :: T0_sea_ice,delta_t_sea_ice
    real(4)                                        :: delta_p
    real(4)                                        :: delta_q
    real(4)                                        :: W0,delta_w

    integer(4)                                :: amsu_channel_loaded       = -1
    integer(4)                                :: amsu_channel_loaded_cloud = -1
    integer(4)                                :: amsu_channel_loaded_emiss = -1

contains

    subroutine read_abs_table_q_AMSU(amsu_channel,path_to_data,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err
        character(len = 200)                  :: file
        integer(4)                            :: numt,nump,numq
        integer(4)                            :: ivap,ioxy,channel
        

        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/abs_tables/amsu_',i2.2,'_abs_table_q_per_Pa.dat')
    
        open(unit= 15,file = file, access='stream', form='unformatted')
        read(15)numt,nump,numq,T0,Delta_t,Delta_p,Delta_q,ivap,ioxy,channel
        read(15)abs_table_q
        close(15)
        abs_table_q(:,0,:) = 0.0

        amsu_channel_loaded = amsu_channel
        err = 0
        return
    end subroutine read_abs_table_q_AMSU

    subroutine read_abs_table_q_AMSU_netcdf(amsu_channel,path_to_data,ivap,ioxy,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(in)                 :: ivap,ioxy
        integer(4),intent(out)                :: err
        character(len = 200)                  :: file
        integer(4)                            :: numt,nump,numq
        integer(4)                            :: channel
        integer(4)                            :: ncid,varid
        integer(4)                            :: status
        integer(4)                            :: t_index,p_index,q_index
        real(4),dimension(0:num_t)            :: t_vals
        real(4),dimension(0:num_p)            :: p_vals
        real(4),dimension(0:num_q)            :: q_vals
        real(4),dimension(0:num_q,0:num_p,0:num_t) :: abs_table_q_netcdf
        

        write(file,100) trim(path_to_data),amsu_channel,ivap,ioxy
100     format(a,'/abs_tables/amsu_',i2.2,'_abs_table_per_Pa_q.',i1.1,'.',i1.1,'.nc')
        print *, "Reading AMSU absorption table from netCDF file: "
        print *, trim(file)
    
        status = nf90_open(trim(file), nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            err = status
            return
        endif



        status = nf90_inq_varid(ncid, 'T', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif


        status = nf90_get_var(ncid, varid, t_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif


        status = nf90_inq_varid(ncid, 'p', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, p_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_inq_varid(ncid, 'q', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, q_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_inq_varid(ncid, 'absorption_per_Pa', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, abs_table_q_netcdf)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_close(ncid)

        numt = size(t_vals) - 1
        nump = size(p_vals) - 1
        numq = size(q_vals) - 1
        channel = amsu_channel

        T0 = t_vals(0)
        Delta_t = t_vals(1) - t_vals(0)
        Delta_p = p_vals(1) - p_vals(0)
        Delta_q = q_vals(1) - q_vals(0)

        ! NetCDF uses (q, p, t); map into (t, p, q)
        do t_index = 0, num_t
            do p_index = 0, num_p
                do q_index = 0, num_q
                    abs_table_q(t_index,p_index,q_index) = abs_table_q_netcdf(q_index,p_index,t_index)
                end do
            end do
        end do
        abs_table_q(:,0,:) = 0.0

        amsu_channel_loaded = amsu_channel
        err = 0
        return
    end subroutine read_abs_table_q_AMSU_netcdf

    
    subroutine read_cld_abs_table_AMSU(amsu_channel,path_to_data,err)


        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err

        character(len = 200)                :: file
        integer(4)                            :: numt,nump,numq
        integer(4)                            :: ivap,ioxy,channel
        

        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/abs_tables/amsu_',i2.2,'_cld_abs_table.dat')

        open(unit= 15,file = file, access='stream', form='unformatted')
        read(15)T0,Delta_t,numt
        read(15)cld_abs_table
        close(15)

        amsu_channel_loaded_cloud = amsu_channel
        err = 0
        return
    end subroutine read_cld_abs_table_AMSU

    subroutine read_cld_abs_table_AMSU_netcdf(amsu_channel,path_to_data,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err

        character(len = 200)                  :: file
        integer(4)                            :: ncid,varid
        integer(4)                            :: status
        real(4),dimension(0:num_t)            :: t_vals

        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/abs_tables/amsu_',i2.2,'_cld_abs_table.nc')

        print *, "Reading AMSU cloud absorption table from netCDF file: "
        print *, trim(file)

        status = nf90_open(trim(file), nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            err = status
            return
        endif

        status = nf90_inq_varid(ncid, 'temperature', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, t_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_inq_varid(ncid, 'absorptivity', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, cld_abs_table)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_close(ncid)

        T0 = t_vals(0)
        Delta_t = t_vals(1) - t_vals(0)

        amsu_channel_loaded_cloud = amsu_channel
        err = 0
        return

    end subroutine read_cld_abs_table_AMSU_netcdf

    subroutine read_ocean_emiss_table_AMSU(amsu_channel,path_to_data,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err
        character(len = 200)                :: file
        integer(4)                            :: numt,numw,channel
        
        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/emiss_tables/amsu_',i2.2,'_emiss_table_W.dat')

        open(unit= 15,file = file, access='stream', form='unformatted')
        read(15)numt,T0,Delta_t,numw,W0,Delta_W,channel
        read(15)ocean_emiss_table
        close(15)

        amsu_channel_loaded_emiss = amsu_channel
        err = 0
        return

    end subroutine read_ocean_emiss_table_AMSU

    subroutine read_ocean_emiss_table_AMSU_netcdf(amsu_channel,path_to_data,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err
        character(len = 200)                  :: file
        integer(4)                            :: channel
        integer(4)                            :: ncid,varid
        integer(4)                            :: status
        integer(4)                            :: t_index,w_index,fov_index
        real(4),dimension(0:num_t)            :: t_vals
        real(4),dimension(0:num_w)            :: w_vals
        real(4),dimension(0:num_fov-1)        :: fov_vals
        !real(4),dimension(0:num_t,0:num_w,0:num_fov-1) :: ocean_emiss_netcdf
        real(4),dimension(0:num_fov-1,0:num_w,0:num_t) :: ocean_emiss_netcdf

        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/ocean_emiss_tables/ocean_emissivity_table_AMSU_channel_',i2.2,'.nc')
        print *, "Reading AMSU ocean emissivity table from netCDF file: "
        print *, trim(file)

        status = nf90_open(trim(file), nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            err = status
            return
        endif

        status = nf90_inq_varid(ncid, 'temperature', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_get_var(ncid, varid, t_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_inq_varid(ncid, 'wind_speed', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_get_var(ncid, varid, w_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

         status = nf90_inq_varid(ncid, 'fov', varid)
        status = nf90_inq_varid(ncid, 'fov', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_get_var(ncid, varid, fov_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        
        status = nf90_inq_varid(ncid, 'emissivity', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_get_var(ncid, varid, ocean_emiss_netcdf)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_close(ncid)

        channel = amsu_channel
        T0 = t_vals(0)
        Delta_t = t_vals(1) - t_vals(0)
        W0 = w_vals(0)
        Delta_w = w_vals(1) - w_vals(0)

        ! NetCDF uses (temperature, wind_speed, fov); map fov index 0..14 -> 1..15
        do t_index = 0, num_t
            do w_index = 0, num_w
                do fov_index = 0, num_fov-1
                    ocean_emiss_table(t_index,w_index,fov_index+1) = ocean_emiss_netcdf(fov_index,w_index,t_index)
                end do
            end do
        end do

        amsu_channel_loaded_emiss = amsu_channel
        err = 0
        return

    end subroutine read_ocean_emiss_table_AMSU_netcdf
    
    subroutine read_sea_ice_emiss_table_AMSU(amsu_channel,path_to_data,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err
        character(len = 200)                  :: file
        integer(4)                            :: numt,channel
        
        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/emiss_tables/amsu_',i2.2,'_emiss_table_sea_ice.dat')

        open(unit= 15,file = file, access='stream', form='unformatted')
        read(15)numt,T0_sea_ice,Delta_t_sea_ice
        read(15)sea_ice_emiss_table
        close(15)

        amsu_channel_loaded_emiss = amsu_channel
        err = 0
        return

    end subroutine read_sea_ice_emiss_table_AMSU

    subroutine read_sea_ice_emiss_table_AMSU_netcdf(amsu_channel,path_to_data,err)

        implicit none

        integer(4),intent(in)                 :: amsu_channel
        character(len = *), intent(in)        :: path_to_data
        integer(4),intent(out)                :: err
        character(len = 200)                  :: file
        integer(4)                            :: ncid,varid
        integer(4)                            :: status
        integer(4)                            :: t_index,fov_index
        real(4),dimension(0:num_t)            :: t_vals
        real(4),dimension(0:num_fov-1)        :: fov_vals
        real(4),dimension(0:num_t,0:num_fov-1) :: sea_ice_emiss_netcdf

        write(file,100) trim(path_to_data),amsu_channel
100     format(a,'/emiss_tables/amsu_',i2.2,'_emiss_table_sea_ice.nc')

        print *, "Reading AMSU sea ice emissivity table from netCDF file: "
        print *, trim(file)

        status = nf90_open(trim(file), nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, nf90_strerror(status)
            err = status
            return
        endif

        status = nf90_inq_varid(ncid, 'temperature', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, t_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_inq_varid(ncid, 'fov', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif
        status = nf90_get_var(ncid, varid, fov_vals)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif


        status = nf90_inq_varid(ncid, 'emissivity', varid)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif

        status = nf90_get_var(ncid, varid, sea_ice_emiss_netcdf)
        if (status /= nf90_noerr) then
            err = status
            print *, nf90_strerror(status)
            status = nf90_close(ncid)
            return
        endif


        status = nf90_close(ncid)

        T0_sea_ice = t_vals(0)
        Delta_t_sea_ice = t_vals(1) - t_vals(0)

        ! NetCDF uses (fov, temperature); map fov index 0..14 -> 1..15
        do t_index = 0, num_t
            do fov_index = 0, num_fov-1
                sea_ice_emiss_table(t_index,fov_index+1) = sea_ice_emiss_netcdf(t_index,fov_index)
            end do
        end do

        amsu_channel_loaded_emiss = amsu_channel
        err = 0
        return

    end subroutine read_sea_ice_emiss_table_AMSU_netcdf
    
     
    real(4) function find_abs_q_AMSU(t,p,q,error)

        ! this routine finds the absorption coefficient (in nepers/Pa) as
        ! a function of temperature t (Kelvin), pressure p (hPa), and specific humidity q (kg/kg)
        ! by interpolating a table.  The subroutine read_abs_table_q must be called before this
        ! routine is used.  

        use, intrinsic :: ieee_arithmetic

        real(4)                :: t    ! temperature
        real(4)                :: p    ! pressure
        real(4)                :: q    ! spec. humidity (kg h20/kg (h20 + dry air))
        integer(4)            :: error

        integer(4)            :: t1,p1,q1
        real(4)                :: t_scaled,p_scaled,q_scaled
        real(4)                :: wt,wp,wq

        real(4)                :: abs_interp

        if (amsu_channel_loaded < 1) then
            error = -1
            find_abs_q_AMSU = ieee_value(0.0, ieee_quiet_nan)
            return
        endif

        error = 0

        t_scaled = (t-t0)/delta_t
        t1 = floor(t_scaled)
        wt = t_scaled - t1

        if (t1 < 0) then
            t1 = 0
            wt = 0.0
            error = ior(error,1)
        endif

        if (t1 > (num_t - 1)) then   ! extrapolate -- set error to +1
            t1 = num_t -1
            wt = t_scaled - t1
            error = ior(error,1)
        endif

        p_scaled = p/delta_p
        p1 = floor(p_scaled)
        wp = p_scaled - p1

        if (p1 < 0) then 
            p1 = 0
            wp = 0.0
            error = ior(error,2)
        endif

        if (p1 > (num_p-1)) then   ! extrapolate
            p1 =num_p-1
            wp = p_scaled - p1
            error = ior(error,2)
        endif

        q_scaled = q/delta_q
        q1 = floor(q_scaled)
        wq = q_scaled - q1

        if (q1 < 0) then 
            q1 = 0
            wq = 0.0
            error = ior(error,4)
        endif

        if (q1 > (num_q-1)) then   ! extrapolate
            q1 = num_q-1
            wq = q_Scaled - q1
            error = ior(error,4)
        endif

        abs_interp = (1.0-wt)*((1.0-wp)*((1.0-wq)*abs_table_q(t1,p1,    q1) + wq*abs_table_q(t1  ,p1,  q1+1)) +     &
                                     wp*((1.0-wq)*abs_table_q(t1,p1+1,  q1) + wq*abs_table_q(t1  ,p1+1,q1+1))) +    &
                           wt*((1.0-wp)*((1.0-wq)*abs_table_q(t1+1,p1,  q1) + wq*abs_table_q(t1+1,p1,  q1+1)) +     &
                                     wp*((1.0-wq)*abs_table_q(t1+1,p1+1,q1) + wq*abs_table_q(t1+1,p1+1,q1+1)))

        find_abs_q_AMSU = abs_interp

        return

    end function find_abs_q_AMSU

    real(4) function find_cld_abs_AMSU(T,rho_cld,error)
        
        real(4)                    :: T
        real(4)                    :: rho_cld
        integer(4)                :: error

        real(4)                    :: t_scaled
        integer(4)                :: t1
        real(4)                    :: wt

        real(4)                    :: abs_interp

        t_scaled = (t-t0)/delta_t
        t1 = floor(t_scaled)
        wt = t_scaled - t1

        if (t1 < 0) then
            t1 = 0 
            wt = 0.0
            error = 1
        endif

        if (t1 > num_t-1) then ! extrapolate
            t1 = num_t-1
            wt = t_scaled - t1
            error = 1
        endif 

        abs_interp = (1.0-wt)*cld_abs_table(t1) + wt*cld_abs_table(t1+1)

        ! now multiply by rho in g/cm^3 to get absorbtion in nepers/cm

        abs_interp = abs_interp*rho_cld*0.001

        find_cld_abs_AMSU = abs_interp*100.0

        return

    end function find_cld_abs_AMSU

    real(4) function find_emiss_AMSU(T,W,fov,error)

        real(4)                    :: T
        real(4)                    :: W
        integer(4)                :: fov
        integer(4)                :: error

        real(4)                    :: t_scaled
        integer(4)                :: t1
        real(4)                    :: wt
        real(4)                    :: w_scaled
        integer(4)                :: w1
        real(4)                    :: ww

        real(4)                    :: emiss_interp

        t_scaled = (t-t0)/delta_t
        t1 = floor(t_scaled)
        wt = t_scaled - t1

        if (t1 < 0) then
            t1 = 0 
            wt = 0.0
            error = 1
        endif

        if (t1 > num_t-1) then ! extrapolate
            t1 = num_t-1
            wt = t_scaled - t1
            error = 1
        endif

        w_scaled = (w-w0)/delta_w
        w1 = floor(w_scaled)
        ww = w_scaled - w1

        if (w1 < 0) then
            w1 = 0 
            ww = 0.0
            error = 1
        endif

        if (w1 > num_w-1) then ! extrapolate
            w1 = num_w-1
            ww = w_scaled - w1
            error = 1
        endif

        emiss_interp = (1.0-wt)*((1.0-ww)*ocean_emiss_table(t1,  w1,fov) + ww*ocean_emiss_table(t1  ,w1+1,fov)) +     &
                             wt*((1.0-ww)*ocean_emiss_table(t1+1,w1,fov) + ww*ocean_emiss_table(t1+1,w1+1,fov))

        find_emiss_AMSU = emiss_interp

    end function find_emiss_AMSU
    
    real(4) function find_emiss_sea_ice_AMSU(T,fov,error)

        real(4)                    :: T
        integer(4)                :: fov
        integer(4)                :: error

        real(4)                    :: t_scaled
        integer(4)                :: t1
        real(4)                    :: wt
        real(4)                    :: w_scaled
        integer(4)                :: w1
        real(4)                    :: ww

        real(4)                    :: emiss_interp

        t_scaled = (t-t0_sea_ice)/delta_t_sea_ice
        t1 = floor(t_scaled)
        wt = t_scaled - t1

        if (t1 < 0) then
            t1 = 0 
            wt = 0.0
            error = 1
        endif

        if (t1 > num_t-1) then ! extrapolate
            t1 = num_t-1
            wt = t_scaled - t1
            error = 1
        endif



        emiss_interp = (1.0-wt)*sea_ice_emiss_table(t1,fov) + wt*sea_ice_emiss_table(t1+1,fov) 

        find_emiss_sea_ice_AMSU = emiss_interp

    end function find_emiss_sea_ice_AMSU

end module rtm_tables_AMSU

        
