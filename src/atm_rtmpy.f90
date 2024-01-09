module atm_rtm

    use trig_degrees, only: cosd

    !use rtm_atmosphere, only: get_atm_components
    real(4), parameter :: rd=287.05,epsilon=0.622,rv=rd/epsilon
    real(4), parameter :: one_minus_epsilon=1.0-epsilon,epsilon_ratio=one_minus_epsilon/epsilon
    
contains
    

    subroutine abs_o2_rosen_2017(n,t,p,pv,f,o2abs_out_rosen)

        use o2abs_19, only: o2abs

        implicit none

        integer(4) :: n
        real(4), dimension(n), intent(in) :: t,p,pv,f
        real(4), dimension(n) :: o2abs_out_rosen

        integer(4) :: i
        real(4) :: vapor_density
        
        do i=1,n   
            ! convert from pv to vapor density in g/m**3
            ! which obabs expects
            vapor_density = (1.0e5/rv)*pv(i)/t(i)
            o2abs_out_rosen(i) = o2abs(t(i),p(i),vapor_density,f(i))
        enddo

    end subroutine abs_o2_rosen_2017

    subroutine abs_o2_rss_2022(n,t,p,pv,f,o2abs_rss_out)

        use rss_2022_absorption, only: absoxy_rss_2022_func

        implicit none

        integer(4), intent(in) :: n
        real(4), dimension(n), intent(in) :: t,p,pv,f
        real(4), dimension(n), intent(inout) :: o2abs_rss_out

        integer(4) :: i

        do i=1,n
            o2abs_rss_out(i) = absoxy_rss_2022_func(p(i),t(i),pv(i),f(i))
            !print *, 'i,t,p,pv,f,o2abs_rss_out',i, t(i),p(i),pv(i),f(i),o2abs_rss_out(i)
        enddo

    end subroutine abs_o2_rss_2022

    subroutine abs_h2o_rss_2022(n,t,p,pv,f,h2oabs_rss_out)

        use rss_2022_absorption, only: absh2o_rss_2022

        implicit none

        integer(4), intent(in) :: n
        real(4), dimension(n), intent(in) :: t,p,pv,f
        real(4), dimension(n), intent(inout) :: h2oabs_rss_out

        integer(4) :: i
        real(4) :: temp

        do i=1,n
            call absh2o_rss_2022(p(i),t(i),pv(i),f(i),  temp)
            h2oabs_rss_out(i) = temp
        enddo

    end subroutine abs_h2o_rss_2022

    subroutine goff_gratch_rss_2022(n,t,p,relhum,pv)

        use goff_gratch, only: goff_gratch_vap 

        implicit none

        integer(4), intent(in) :: n
        real(4), dimension(n), intent(in) :: t,p,relhum
        real(4), dimension(n), intent(inout) :: pv

        integer(4) :: i
        real(4) :: P_V,RHO_V

        do i=1,n
            call goff_gratch_vap(t(i),p(i),relhum(i),P_V,RHO_V)
            pv(i) = P_V
        enddo
    end subroutine goff_gratch_rss_2022

    subroutine get_atm_components(ioxy, nlev, ilevel0, freq, tht, z, &
        p, t, pv, rhol, tran, tbup, tbdw)
        ! calculates atmospheric rtm components tran, tbup, tbdw for given profiles
        ! no rain

        use o2abs_19, only: o2abs
        use rss_2022_absorption, only: absoxy_rss_2022_func,absh2o_rss_2022

        implicit none

        integer(4), intent(in)    :: ioxy        ! 1 = RSS, 2 = Rosenkranz 2017
        integer(4), intent(in)    :: nlev        ! number of levels. numbered: 0,1,2, ... nlev
        integer(4), intent(in)    :: ilevel0     ! index of 1st level to start integration     
        real(4),    intent(in)    :: freq        ! frequency [ghz]
        real(4),    intent(in)    :: tht         ! earth incidence angle [deg]

        real(4), dimension(0:nlev), intent(in)    ::    z           ! elevation above sea level         [m]    

        real(4), dimension(0:nlev), intent(in)    ::    p           ! air pressure profile                [mbar]    
        real(4), dimension(0:nlev), intent(in)    ::    t           ! air temperaturee profile            [K]    
        real(4), dimension(0:nlev), intent(in)    ::    pv          ! water vapor pressure profile        [mbar]
        real(4), dimension(0:nlev), intent(in)    ::    rhol        ! liquid cloud water density profile[g/m**3]


        real(4), intent(inout)                    ::    tran        ! atmospheric transmittance    
        real(4), intent(inout)                    ::    tbup        ! upwelling   atmospheric brightness temperature at TOA [K]    
        real(4), intent(inout)                    ::    tbdw        ! downwelling atmospheric brightness temperature at Surface [K]    

        real(4), dimension(0:nlev)                ::  alpha_v, alpha_o, alpha_l, alpha_t
        integer(4)                                ::  ipr
        real(4)                                   ::  vapor_density        



        !   profile of absorption coefficients
        alpha_o=0.0
        alpha_v=0.0
        alpha_l=0.0
        alpha_t=0.0

        do ipr  = ilevel0,nlev
            if (ioxy == 1) then
                alpha_o(ipr) =  absoxy_rss_2022_func(p(ipr),t(ipr),pv(ipr),freq)
            else if (ioxy == 2) then
                vapor_density = (1.0e5/rv)*pv(ipr)/t(ipr)
                alpha_o(ipr) = o2abs(t(ipr),p(ipr),vapor_density,freq)
            else
                print *, 'ioxy = ',ioxy
                stop 'ioxy must be 1 or 2'
            endif

            call absh2o_rss_2022(p(ipr),t(ipr),pv(ipr),freq, alpha_v(ipr))

            if (rhol(ipr) > 1.0e-7 .and. t(ipr)>230.) then
                call fdcldabs(freq,t(ipr),rhol(ipr),   alpha_l(ipr))
            else
                alpha_l(ipr) = 0.0
            endif

        enddo

        alpha_o = alpha_o/1000.0     ! [neper/m]
        alpha_v = alpha_v/1000.0     ! [neper/m]
        alpha_l = alpha_l/1000.0     ! [neper/m]

        alpha_t=alpha_o+alpha_v+alpha_l

        call atm_tran(nlev-ilevel0, tht, t(ilevel0:nlev), z(ilevel0:nlev), alpha_t(ilevel0:nlev),    tran, tbdw, tbup)
        
    
        return
    end subroutine get_atm_components

    subroutine get_atm_components_rss_2022(nlev, ilevel0, freq, tht, z, &
        p, t, pv, rhol, results)
    
        integer(4), intent(in) :: nlev
        integer(4), intent(in) :: ilevel0
        real(4), intent(in) :: freq
        real(4), intent(in) :: tht

        real(4), intent(in) :: z(0:nlev)
        real(4), intent(in) :: p(0:nlev)
        real(4), intent(in) :: t(0:nlev)
        real(4), intent(in) :: pv(0:nlev)
        real(4), intent(in) :: rhol(0:nlev)
    
        real(4), intent(inout) :: results(3)

        real(4) :: tran
        real(4) :: tbup
        real(4) :: tbdw
    
        integer(4), parameter :: ioxy = 1

        integer(4) :: ilev !debug only
    
        call get_atm_components(ioxy,nlev, ilevel0, freq, tht, z, &
                                p, t, pv, rhol, tran, tbup, tbdw)
        
        results(1) = tran
        results(2) = tbdw
        results(3) = tbup

    end subroutine get_atm_components_rss_2022

    subroutine get_atm_components_rss_2022_2d(nlev,nprofile, freq, tht, z, &
        p, t, pv, rhol,tran,tbup,tbdw)
    
        integer(4), intent(in) :: nlev
        integer(4), intent(in) :: nprofile
        real(4), intent(in) :: freq
        real(4), intent(in) :: tht

        real(4), intent(in) :: z(0:nlev,nprofile)
        real(4), intent(in) :: p(0:nlev,nprofile)
        real(4), intent(in) :: t(0:nlev,nprofile)
        real(4), intent(in) :: pv(0:nlev,nprofile)
        real(4), intent(in) :: rhol(0:nlev,nprofile)
    
        real(4), intent(inout) :: tran(nprofile)
        real(4), intent(inout) :: tbup(nprofile)
        real(4), intent(inout) :: tbdw(nprofile)

        real(4) :: z1(0:nlev)
        real(4) :: p1(0:nlev)
        real(4) :: t1(0:nlev)
        real(4) :: pv1(0:nlev)
        real(4) :: rhol1(0:nlev)
    
        integer(4), parameter :: ioxy = 1 ! This chooses the RSS 2022 absorption model

        integer(4) :: iprofile,ilev,ilevel0
        
        print *,'nlev = ',nlev
        print *,'nprofile = ',nprofile
        ilevel0 = 0
        do iprofile=1,nprofile
            if (modulo(iprofile,100000) == 0) then
                print *,'get_atm_components_rss_2022_2d: iprofile = ',iprofile
            endif
            ! make local copy of the profile to process
            z1 = z(:,iprofile)
            p1 = p(:,iprofile)
            t1 = t(:,iprofile)
            pv1 = pv(:,iprofile)
            rhol1 = rhol(:,iprofile)
            ilevel0 = 0

            ! insert the surface data at the right level and set
            ! ilevel0 to the index of the surface data and ignore levels
            ! below this level
            do ilev = nlev,1,-1
                if (z1(0)>z1(ilev)) then
                    !replace this level with the surface data
                    p1(ilev) = p1(0)
                    t1(ilev) = t1(0)
                    pv1(ilev) = pv1(0)
                    rhol1(ilev) = rhol1(0)
                    z1(ilev) = z1(0)
                    ilevel0 = ilev
                    !print *,'iprofile = ',iprofile,'replacing level ',ilev,' with surface data'
                    exit
                endif
            enddo  

            !do the atm components calculation for this profile
            call get_atm_components(ioxy,nlev, ilevel0, freq, tht, z1, &
                                    p1, t1, pv1, rhol1, &
                                    tran(iprofile), tbup(iprofile), tbdw(iprofile))
            

        enddo

    end subroutine get_atm_components_rss_2022_2d

    


    subroutine get_atm_components_rosen_2017_2d(nlev,nprofile, freq, tht, z, &
        p, t, pv, rhol,tran,tbup,tbdw)
    
        integer(4), intent(in) :: nlev
        integer(4), intent(in) :: nprofile
        real(4), intent(in) :: freq
        real(4), intent(in) :: tht

        real(4), intent(in) :: z(0:nlev,nprofile)
        real(4), intent(in) :: p(0:nlev,nprofile)
        real(4), intent(in) :: t(0:nlev,nprofile)
        real(4), intent(in) :: pv(0:nlev,nprofile)
        real(4), intent(in) :: rhol(0:nlev,nprofile)
    
        real(4), intent(inout) :: tran(nprofile)
        real(4), intent(inout) :: tbup(nprofile)
        real(4), intent(inout) :: tbdw(nprofile)

        real(4) :: z1(0:nlev)
        real(4) :: p1(0:nlev)
        real(4) :: t1(0:nlev)
        real(4) :: pv1(0:nlev)
        real(4) :: rhol1(0:nlev)
    
        integer(4), parameter :: ioxy = 2 ! this chooses the Rosenkranz 2017 model

        integer(4) :: iprofile,ilev,ilevel0
        
        print *,'nlev = ',nlev
        print *,'nprofile = ',nprofile
        ilevel0 = 0
        do iprofile=1,nprofile
            if (modulo(iprofile,100000) == 0) then
                print *,'get_atm_components_rosen_2017_2d: iprofile = ',iprofile
            endif
            ! make local copy of the profile to process
            z1 = z(:,iprofile)
            p1 = p(:,iprofile)
            t1 = t(:,iprofile)
            pv1 = pv(:,iprofile)
            rhol1 = rhol(:,iprofile)
            ilevel0 = 0

            ! insert the surface data at the right level and set
            ! ilevel0 to the index of the surface data and ignore levels
            ! below this level
            do ilev = nlev,1,-1
                if (z1(0)>z1(ilev)) then  !if the surface height (z1(0) is above the current level (z1(ilev):...
                    !replace this level with the surface data
                    p1(ilev) = p1(0)
                    t1(ilev) = t1(0)
                    pv1(ilev) = pv1(0)
                    rhol1(ilev) = rhol1(0)
                    z1(ilev) = z1(0)
                    ilevel0 = ilev
                    !print *,'iprofile = ',iprofile,'replacing level ',ilev,' with surface data'
                    
                    exit
                endif
            enddo  

            !do the atm components calculation for this profile
            call get_atm_components(ioxy,nlev, ilevel0, freq, tht, z1, &
                                    p1, t1, pv1, rhol1, &
                                    tran(iprofile), tbup(iprofile), tbdw(iprofile))
            

        enddo

    end subroutine get_atm_components_rosen_2017_2d
    
    subroutine get_atm_components_rosen_2017(nlev, ilevel0, freq, tht, z, &
        p, t, pv, rhol, results)
    
        integer(4), intent(in) :: nlev
        integer(4), intent(in) :: ilevel0
        real(4), intent(in) :: freq
        real(4), intent(in) :: tht
    
        real(4), intent(in) :: z(0:nlev)
        real(4), intent(in) :: p(0:nlev)
        real(4), intent(in) :: t(0:nlev)
        real(4), intent(in) :: pv(0:nlev)
        real(4), intent(in) :: rhol(0:nlev)
    
        real(4), intent(inout) :: results(3)
        
        real(4) :: tran
        real(4) :: tbup
        real(4) :: tbdw
    
        integer(4), parameter :: ioxy = 2
    
        call get_atm_components(ioxy, nlev, ilevel0, freq, tht, z, &
                                p, t, pv, rhol, tran, tbup, tbdw)

        results(1) = tran
        results(2) = tbdw
        results(3) = tbup
    
    end subroutine get_atm_components_rosen_2017

    subroutine atm_tran_multiple_profiles(nlev,nprofile,tht,t_arr,z_arr,tabs_arr,tran_arr,tbdw_arr,tbup_arr)
        implicit none

        integer(4), intent(in) :: nlev
        integer(4), intent(in) :: nprofile
        real(4), intent(in) :: tht
        real(4), intent(in) :: t_arr(0:nlev,nprofile)
        real(4), intent(in) :: z_arr(0:nlev,nprofile)
        real(4), intent(in) :: tabs_arr(0:nlev,nprofile)
        real(4), intent(inout) :: tran_arr(nprofile)
        real(4), intent(inout) :: tbdw_arr(nprofile)
        real(4), intent(inout) :: tbup_arr(nprofile)

        real(4) :: z1(0:nlev)
        real(4) :: t1(0:nlev)
        real(4) :: tabs1(0:nlev)

        integer(4) :: iprofile,ilevel0,ilev

        do iprofile=1,nprofile
            
            z1 = z_arr(:,iprofile)
            t1 = t_arr(:,iprofile)
            tabs1 = tabs_arr(:,iprofile)
            ilevel0 = 0
            
            ! insert the surface data at the right level and set
            ! ilevel0 to the index of the surface data and ignore levels
            ! below this level
            do ilev = nlev,1,-1
                if (z1(0)>z1(ilev)) then
                    !replace this level with the surface data
                    t1(ilev) = t1(0)
                    tabs1(ilev) = tabs1(0)
                    z1(ilev) = z1(0)
                    ilevel0 = ilev
                    exit
                endif
            enddo 
            
            call atm_tran(nlev-ilevel0, tht, &
                t1(ilevel0:nlev), z1(ilevel0:nlev), tabs1(ilevel0:nlev), &
                tran_arr(iprofile), tbdw_arr(iprofile), tbup_arr(iprofile))

            if (iprofile .eq. 2809) then
                print *,'iprofile = ',iprofile
                print *,'start_level = ',ilevel0
                print *,'z1 = ',z1
                print *,'t1 = ',t1
                print *,'tabs1 = ',tabs1
                print *
                print *,'tran = ',tran_arr(iprofile)
            endif
    

        enddo

    end subroutine atm_tran_multiple_profiles

    subroutine atm_tran(nlev,tht,t,z,tabs,tran,tbdw,tbup)
        implicit none
    
        !       17 january 2007
        !       in expression for dsdh i use the actual elevation rather than frank's delta=0.00035
    
        !        february 2006
        !        substitute sectht by dsdh in order to include grazing angles
    
    
        !         computer atmospheric downwelling and upwelling brightness temperatures
        !         and upward transmittance at each pressure level (altitude) 
    
        !
        !        input:
        !             nlev           number of atmosphere levels
        !             tht            earth incidence angle [in deg]
        !             tabs(0:nlev)   atmospheric absorptrion coefficients [nepers/m]
        !             t(0:nlev)      temperature profile[in k]
        !             z(0:nlev)      elevation (m) 
    
        !        output:
        !            tran            total atmospheric transmission
        !            tbdw            downwelling brightness temperature t_bd [in k]
        !            tbup            upwelling   brightness temperature t_bu [in k]  
    
    
        integer(4)            :: nlev,i
        real(4)                :: t(0:nlev),z(0:nlev),tabs(0:nlev)
        real(4)                :: opacty(nlev),tavg(nlev),ems(nlev)
        real(4), parameter    :: re=6378.135 ! earth radius [km]
        real(4)                :: sumop, sumdw, sumup, tran ,tbavg, tbdw, tbup, dsdh, tht, h, delta
    
        do i=1,nlev
            if (z(i) <= z(i-1)) then
                write(*,*) z
                stop 'atm_tran: error in profile order, pgm stopped.'    
            endif
            h = 0.5*(z(i)+z(i-1))/1000.  ! average elevation of profile layer in km
            delta=h/re
            if (delta <= 0) then
                delta=0.0
                ! print *,'h = ',h
                ! print *,'delta = ',delta
                ! stop ' error in computing delta in routine atm_tran. pgm stopped.'
            endif
            dsdh = (1.0+delta)/sqrt(cosd(tht)*cosd(tht) + delta*(2+delta))
            opacty(i)=-dsdh*0.5*(tabs(i-1)+tabs(i))*(z(i)-z(i-1))
            tavg(i)  =0.5*(t(i-1)+t(i))
            ems(i)   =1.-exp(opacty(i))
        enddo
    
        sumop=0 
        sumdw=0
        do i=1,nlev
            sumdw=sumdw+(tavg(i)-t(1))*ems(i)*exp(sumop)
            sumop=sumop+opacty(i)
        enddo
    
        sumop=0 
        sumup=0.
        do i=nlev,1,-1
            sumup=sumup+(tavg(i)-t(1))*ems(i)*exp(sumop)
            sumop=sumop+opacty(i)
        enddo
    
        tran=exp(sumop)
        tbavg=(1.-tran)*t(1)
        tbdw=tbavg+sumdw
        tbup=tbavg+sumup
    
        return
    end subroutine atm_tran


    subroutine fdcldabs(freq,t,rhol,al)
        !    liquid cloud water absorption
        !    rayleigh (no rain)
        !    freq:      frequency [ghz]
        !    t:         temperature [k]
        !    rhol:      liquid cloud water density [g/m**3]

        !    output:
        !    al:        cloud water absorption coefficient [neper/km]

        implicit none 

        real(4), intent(in)    :: freq,t,rhol
        real(4), intent(out)   :: al

        real(4)                :: rhol0, wavlen


        complex(4)             :: permit

        real(4), parameter     :: c=29.979, pi=3.14159
        ! bug in c fixed (see Richard's email 08/10/2023)

        rhol0 = 1.0e-6*rhol ![g/cm**3]
        call find_permittivity_meissner_wentz(freq,t,0.0,   permit) 

        wavlen = c/freq

        al = (6.0*pi*rhol0/wavlen)*aimag((1.0-permit)/(2.0+permit))  ! nepers/cm
        al = 1.0e5*al ! nepers/km

        return
    end subroutine fdcldabs

    subroutine fdcldabs_2d_ql(nlev,nprofile,freq,p,t,q,ql,al)
        !    liquid cloud water absorption
        !    rayleigh (no rain)
        !    freq:      frequency [ghz]
        !    p(0:nlev,nprofiles):         pressure [mbar]
        !    t(0:nlev,nprofiles):         temperature [k]
        !    q(0:nlev,nprofiles):         specific humidity [kg/kg]
        !    ql(0:nlev,nprofiles):        specific liquid cloud water content [kg/kg]


        !    output:
        !    al(0:nlev,nprofiles):        cloud water absorption coefficient [neper/km]

        implicit none 
        integer(4), intent(in) :: nlev
        integer(4), intent(in) :: nprofile
        real(4), intent(in) :: freq
        real(4), intent(in) :: p(0:nlev,nprofile)
        real(4), intent(in) :: t(0:nlev,nprofile)
        real(4), intent(in) :: q(0:nlev,nprofile)
        real(4), intent(in) :: ql(0:nlev,nprofile)
        real(4), intent(inout) :: al(0:nlev,nprofile)

        integer(4) :: iprofile,ilev
        real(4) :: rhol,rm

        do iprofile=1,nprofile
            do ilev=0,nlev
                ! convert from ql (in kg/kg) to vapor density (in g/m**3)
                ! which fdcldabs expects
                ! calculate R moist 
                rm = rd*(1+epsilon_ratio*q(ilev,iprofile))
                !calculate rho cloud liquid water
                rhol = ql(ilev,iprofile)*100.0*p(ilev,iprofile)/(rm*t(ilev,iprofile))
                call fdcldabs(freq,t(ilev,iprofile),rhol,al(ilev,iprofile))
            enddo
        enddo

    end subroutine fdcldabs_2d_ql


    subroutine column_integral(nlevel,z,rho,ip,  col,ierr)
        !     columnar integral
        !     input:  
        !     nlevel: number of profiles (nlevel = 0 -> sfc)
        !     z:      altitude levels [in m]
        !     rho:    density [in unit/m]      (unit arbitrary)
        !     ip:     =1 : linear varying profile -> arithmetic mean for integration
        !             =2 : approximately exponentially varying profile ->
        !                  use average of arithmetic and geometric mean for integration
        !             this has smaller error in integration than either 
        !             arithmetic (a) or geometric (g) mean
        !             (if profile is varying exactly exponentially with z, then the integration error
        !              is minimal for the combination: 2/3 g + 1/3 a)
        !
        !     output:
        !     col [in unit]       
        !
        implicit none
        integer(4)                        ::  ip, nlevel, i ,ierr
        real(4)                            ::  col, dz , avg
        real(4), dimension(0:nlevel)    ::  z,rho
        !
            ierr = 0
            col = 0.
            profiles: do i=1,nlevel
                if (z(i).le.z(i-1)) then
                    write(*,*) ' wrong order in column. stop.'
                    stop
                endif
                dz = z(i) - z(i-1)
                if (ip.eq.1) then
                    avg = 0.5* (rho(i) + rho(i-1) )
                else if (ip.eq.2) then
                    avg = 0.25* ( rho(i) + rho(i-1) + 2.*sqrt(rho(i-1)*rho(i)) )
                endif

                col = col + avg*dz
            enddo profiles

        return
    end    subroutine column_integral

    !   input:
    !   name   parameter  unit  range
    !   freq   frequency  [ghz] 1 to 400
    !   t      sst        [k]   248.16 k (-25 c) to 313.16 k (40 c) for pure water
    !                           271.16 k (-2  c) to 307.16 k (34 c) for saline water
    !   s      salinity   [ppt]  0 to 40
    !
    !   output:
    !   eps    complex dielectric constant 
    !          negative imaginary part to be consistent with wentz1 conventionc
    !
    !
    subroutine    find_permittivity_meissner_wentz(freq,t,s,   eps)
        implicit none
            
            real(4), parameter :: f0=17.97510
            real(4) freq,t,sst,s
            real(4) e0s,e1s,e2s,n1s,n2s,sig
            complex(4) :: j = (0.0,1.0), eps

            sst = t - 273.15 ! [celsius]
            call dielectric_meissner_wentz(sst,s,  e0s,e1s,e2s,n1s,n2s,sig)

        !     debye law (2 relaxation wavelengths)
            eps = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s +  j*sig*f0/freq
            eps = conjg(eps)

        return 
    end subroutine find_permittivity_meissner_wentz    


    subroutine dielectric_meissner_wentz(sst_in,s,   e0s,e1s,e2s,n1s,n2s,sig)
        !
        !     complex dielectric constant: eps
        !     [MW 2004, MW 2012].
        !     
        !     Changes from [MW 2012]:
        !     1. Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
        !     2. Changed SST behavior of coefficient b2 from:
        !     b2 = 1.0 + s*(z(10) + z(11)*sst) to
        !     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
        !
        !!
        !     input:
        !     name   parameter  unit  range
        !     sst      sst        [c]   -25 c to 40 c for pure water
        !                               -2  c to 34 c for saline water
        !     s      salinity   [ppt]  0 to 40
        !
        !     output:
        !     eps    complex dielectric constant
        !            negative imaginary part to be consistent with wentz1 convention
        !

    implicit none


        real(4), intent(in)  :: sst_in,s
        real(4), intent(out) :: e0s,e1s,e2s,n1s,n2s,sig
    
        real(4), dimension(11), parameter :: &
        x=(/ 5.7230e+00, 2.2379e-02, -7.1237e-04, 5.0478e+00, -7.0315e-02, 6.0059e-04, 3.6143e+00, &
            2.8841e-02, 1.3652e-01,  1.4825e-03, 2.4166e-04 /)
        
        real(4), dimension(13), parameter :: &
        z=(/ -3.56417e-03,  4.74868e-06,  1.15574e-05,  2.39357e-03, -3.13530e-05, &
                2.52477e-07, -6.28908e-03,  1.76032e-04, -9.22144e-05, -1.99723e-02, &
                1.81176e-04, -2.04265e-03,  1.57883e-04  /)  ! 2004

        real(4), dimension(3), parameter :: a0coef=(/ -0.33330E-02,  4.74868e-06,  0.0e+00/)
        real(4), dimension(5), parameter :: b1coef=(/0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08/)
    
        real(4) :: e0,e1,e2,n1,n2
        real(4) :: a0,a1,a2,b1,b2
        real(4) :: sig35,r15,rtr15,alpha0,alpha1

        real(4) :: sst,sst2,sst3,sst4,s2
        
        sst=sst_in
        if(sst.lt.-30.16) sst=-30.16  !protects against n1 and n2 going zero for very cold water
        
        sst2=sst*sst
        sst3=sst2*sst
        sst4=sst3*sst

        s2=s*s
    
        !     pure water
        e0    = (3.70886e4 - 8.2168e1*sst)/(4.21854e2 + sst) ! stogryn et al.
        e1    = x(1) + x(2)*sst + x(3)*sst2
        n1    = (45.00 + sst)/(x(4) + x(5)*sst + x(6)*sst2)
        e2    = x(7) + x(8)*sst
        n2    = (45.00 + sst)/(x(9) + x(10)*sst + x(11)*sst2)
        
        !     saline water
        !     conductivity [s/m] taken from stogryn et al.
        sig35 = 2.903602 + 8.60700e-2*sst + 4.738817e-4*sst2 - 2.9910e-6*sst3 + 4.3047e-9*sst4
        r15   = s*(37.5109+5.45216*s+1.4409e-2*s2)/(1004.75+182.283*s+s2)
        alpha0 = (6.9431+3.2841*s-9.9486e-2*s2)/(84.850+69.024*s+s2)
        alpha1 = 49.843 - 0.2276*s + 0.198e-2*s2
        rtr15 = 1.0 + (sst-15.0)*alpha0/(alpha1+sst)
        
        sig = sig35*r15*rtr15
        
        !    permittivity
        a0 = exp(a0coef(1)*s + a0coef(2)*s2 + a0coef(3)*s*sst)  
        e0s = a0*e0
        
        if(sst.le.30) then
            b1 = 1.0 + s*(b1coef(1) + b1coef(2)*sst + b1coef(3)*sst2 + b1coef(4)*sst3 + b1coef(5)*sst4)
        else
            b1 = 1.0 + s*(9.1873715e-04 + 1.5012396e-04*(sst-30))
        endif
        
        n1s = n1*b1
        
        a1  = exp(z(7)*s + z(8)*s2 + z(9)*s*sst)
        e1s = e1*a1

        b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30))
        n2s = n2*b2
        
        a2 = 1.0  + s*(z(12) + z(13)*sst)
        e2s = e2*a2
        
        return
    end subroutine  dielectric_meissner_wentz

    subroutine trilinear_interpolation(array, xvalues, yvalues, zvalues, interp_points, result)
        real(kind=4), intent(in) :: array(:,:,:)
        real(kind=4), intent(in) :: xvalues(:)
        real(kind=4), intent(in) :: yvalues(:)
        real(kind=4), intent(in) :: zvalues(:)
        real(kind=4), intent(in) :: interp_points(:,:)
        real(kind=4), intent(inout) :: result(:)
    
        integer :: i, j, k, p, num_pts,ipt
        real(kind=4) :: interp_point(3)
        real(kind=4) :: xd, yd, zd
        real(kind=4) :: c00, c01, c10, c11, c0, c1
        integer(kind=4) :: shp(2)

        shp = SHAPE(interp_points)
        print *, 'shp = ',shp
        num_pts = shp(2)
        do ipt = 1,num_pts
            interp_point(1) = interp_points(1,ipt)
            interp_point(2) = interp_points(2,ipt)
            interp_point(3) = interp_points(3,ipt)
            !print *, 'ipt = ',ipt,' interp_point = ',interp_point
            ! Find the indices of the nearest lower grid points
            i = 1
            do while (i < size(xvalues) - 1 .and. interp_point(1) > xvalues(i + 1))
                i = i + 1
            end do
    
            j = 1
            do while (j < size(yvalues) - 1 .and. interp_point(2) > yvalues(j + 1))
                j = j + 1
            end do
        
            k = 1
            do while (k < size(zvalues) - 1 .and. interp_point(3) > zvalues(k + 1))
                k = k + 1
            end do
    
            ! Compute fractional distances
            xd = (interp_point(1) - xvalues(i)) / (xvalues(i + 1) - xvalues(i))
            yd = (interp_point(2) - yvalues(j)) / (yvalues(j + 1) - yvalues(j))
            zd = (interp_point(3) - zvalues(k)) / (zvalues(k + 1) - zvalues(k))
            
            ! Trilinear interpolation with scaling and offset
            c00 = array(i,j,k) * (1.0d0 - xd) + array(i+1,j,k) * xd
            c01 = array(i,j+1,k) * (1.0d0 - xd) + array(i+1,j+1,k) * xd
            c10 = array(i,j,k+1) * (1.0d0 - xd) + array(i+1,j,k+1) * xd
            c11 = array(i,j+1,k+1) * (1.0d0 - xd) + array(i+1,j+1,k+1) * xd
        
            c0 = c00 * (1.0d0 - yd) + c01 * yd
            c1 = c10 * (1.0d0 - yd) + c11 * yd
        
            result(ipt) = c0 * (1.0d0 - zd) + c1 * zd
        enddo
    end subroutine trilinear_interpolation
end module atm_rtm



