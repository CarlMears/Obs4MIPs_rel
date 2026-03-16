module rtm
    implicit none
    
    contains
        function cosd(x)
            real(8) :: cosd
            real(8), intent(in) :: x
            cosd = cos(x * 3.141592653589793d0 / 180.0d0)
        end function cosd

        subroutine atm_tran_2(nlev,tht,t_4,z_4,tabs_4,tran_up_4,tran_dwn_4,tbdw_4,tbup_4)
        !-----------------------------------------------------------------------------
        !    computer atmospheric downwelling and upwelling brightness temperatures
        !    and upward transmittance at each pressure level (altitude) 
        !-----------------------------------------------------------------------------
        !    input:
        !     nlev           number of atmosphere levels
        !     tht            earth incidence angle [in deg]
        !     tabs(0:nlev)   atmosphric absorptrion coefficients [nepers/m]
        !     t(0:nlev)      temperature profile[in k]
        !    z(0:nlev)      elevation (m) 
        !-----------------------------------------------------------------------------
        !     output:
        !     tran            total atmospheric transmission
        !     tbdw            downwelling brightness temperature t_bd [in k]
        !     tbup            upwelling   brightness temperature t_bu [in k]  
        !-----------------------------------------------------------------------------
            implicit none

            integer(4):: nlev,i
            real(4) :: t_4(0:nlev)
            real(4) :: z_4(0:nlev)
            real(4) :: tabs_4(0:nlev)
            real(4) :: tran_up_4(0:nlev)
            real(4) :: tran_dwn_4(0:nlev)
            real(8) :: t(0:nlev),z(0:nlev),tabs(0:nlev),tran_up(0:nlev),tran_dwn(0:nlev)
            real(8) :: tran,sectht,sumop,sumdw,sumup,tbavg,tbdw,tbup,tht_8
            real(4) :: tbdw_4,tbup_4,tht
            real(8) :: opacty(nlev),tavg(nlev),ems(nlev)

            t    = t_4
            z    = z_4
            tabs = tabs_4
            tht_8 = tht
            sectht=1.0d0/cosd(tht_8)

            do i=1,nlev
                opacty(i)= -sectht*0.5*(tabs(i-1)+tabs(i))*(z(i)-z(i-1))
                tavg(i)  = 0.5*(t(i-1)+t(i))
                ems(i)   = 1.-exp(opacty(i))  ! ems(i) is the total emissivity of level 1
                                            ! between temperature levels 0 and 1, etc
            enddo

            sumop=0.0 
            sumdw=0.0
            
            tran_dwn(0) = 1.00
            do i=1,nlev
                tran_dwn(i) = exp(sumop+opacty(i))
                sumdw=sumdw+(tavg(i)-t(1))*ems(i)*tran_dwn(i-1)
                sumop=sumop+opacty(i)
            enddo

            sumop=0 
            sumup=0.
            do i=nlev,1,-1
                tran_up(i)  = exp(sumop)
                sumup=sumup+(tavg(i)-t(1))*ems(i)*tran_up(i)
                sumop=sumop+opacty(i)
            enddo

            tran_up(0) = exp(sumop)

            tran=exp(sumop)
            tbavg=(1.-tran)*t(1)
            tbdw=tbavg+sumdw
            tbup=tbavg+sumup

        !   convert back to real(4)

            tran_up_4   = tran_up
            tran_dwn_4  = tran_dwn
            tbdw_4      = tbdw
            tbup_4        = tbup
            return
            end subroutine atm_tran_2


        subroutine atm_tran_p(nlev,tht,t_4,p_4,tabs_p_4,tran_up_4,tran_dwn_4,tbdw_4,tbup_4)
            !-----------------------------------------------------------------------------
            !    computer atmospheric downwelling and upwelling brightness temperatures
            !    and upward transmittance at each pressure level 
            !-----------------------------------------------------------------------------
            !    input:
            !     nlev           number of atmosphere levels
            !     tht            earth incidence angle [in deg]
            !     tabs_p(0:nlev)   atmospheric absorptrion coefficients [nepers/pa]
            !     t(0:nlev)      temperature profile[in k]
            !     p(0:nlev)      pressure (hpa)  This should be ascending in pressure, i.e., 
            !                                    p(0) is the lowest pressure and p(nlev) surface pressure.
            !-----------------------------------------------------------------------------
            !     output:
            !     tran            total atmospheric transmission
            !     tbdw            downwelling brightness temperature t_bd [in k]
            !     tbup            upwelling brightness temperature t_bu [in k]  
            !-----------------------------------------------------------------------------
            implicit none

            integer(4):: nlev,i
            real(4) :: t_4(0:nlev)
            real(4)    :: p_4(0:nlev)
            real(4) :: tabs_p_4(0:nlev)
            real(4) :: tran_up_4(0:nlev)
            real(4) :: tran_dwn_4(0:nlev)
            real(8) :: t(0:nlev),p(0:nlev),tabs_p(0:nlev),tran_up(0:nlev),tran_dwn(0:nlev)
            real(8) :: tran,sectht,sumop,sumdw,sumup,tbavg,tbdw,tbup,tht_8
            real(4) :: tbdw_4,tbup_4,tht
            real(8) :: opacty(nlev),tavg(nlev),ems(nlev)

            t    = t_4
            p    = p_4
            tabs_p = tabs_p_4
            tht_8 = tht
            sectht=1.0d0/cosd(tht_8)

            do i=1,nlev
                opacty(i)= -sectht*0.5*(tabs_p(i-1)+tabs_p(i))*(p(i)-p(i-1))
                tavg(i)  = 0.5*(t(i-1)+t(i))
                ems(i)   = 1.-exp(opacty(i))  ! ems(i) is the total emissivity of level 1
                                              ! between temperature levels 0 and 1, etc
            enddo

            sumop=0.0 
            sumup=0.0
            tran_up(0) = 1.00
            do i=1,nlev
                tran_up(i) = exp(sumop+opacty(i))
                sumup=sumup+(tavg(i)-t(1))*ems(i)*tran_up(i-1)
                sumop=sumop+opacty(i)
            enddo

            sumop=0 
            sumdw=0.
            do i=nlev,1,-1
                tran_dwn(i)  = exp(sumop)
                sumdw=sumdw+(tavg(i)-t(1))*ems(i)*tran_dwn(i)
                sumop=sumop+opacty(i)
            enddo

            tran_dwn(0) = exp(sumop)
            tran=exp(sumop)
            tbavg=(1.-tran)*t(1)
            tbdw=tbavg+sumdw
            tbup=tbavg+sumup

        !   convert back to real(4)
            tran_up_4   = tran_up
            tran_dwn_4  = tran_dwn
            tbdw_4      = tbdw
            tbup_4      = tbup
            return
    end subroutine atm_tran_p
end module rtm


