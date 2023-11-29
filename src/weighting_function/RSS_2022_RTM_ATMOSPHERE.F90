module rtm_atmosphere

    use rss_2022_absorption, only : absh2o_rss_2022, absoxy_rss_2022
    use o2abs_19, only : O2ABS
    use trig_degrees, only : cosd
    private
    save

    public :: get_atm_components, get_opacity, get_opacity_c, get_teff

contains

    subroutine get_atm_components(ioxy, nlev, ilevel0, freq, tht, z, &
                                  p, t, pv, rhol, tran, tbup, tbdw)
    ! calculates atmospheric rtm components tran, tbup, tbdw for given profiles
    ! no rain
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

        
    real(4), intent(out)                    ::    tran        ! atmospheric transmittance    
    real(4), intent(out)                    ::    tbup        ! upwelling   atmospheric brightness temperature at TOA [K]    
    real(4), intent(out)                    ::    tbdw        ! downwelling atmospheric brightness temperature at Surface [K]    

    real(4), dimension(0:nlev)                ::  alpha_v, alpha_o, alpha_l, alpha_t
    integer(4)                                ::    ipr



    !   profile of absorption coefficients
    alpha_o=0.0
    alpha_v=0.0
    alpha_l=0.0
    alpha_t=0.0


    do ipr  = ilevel0,nlev
        call  fdabscoeff(ioxy, freq, p(ipr), t(ipr), pv(ipr), alpha_o(ipr), alpha_v(ipr))   ! [neper/km]
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



    subroutine get_opacity(ioxy, nlev, ilevel0, freq, z, p, t, pv, rhol,        ao, av, al, opacity)
    ! calculates atmospheric path length components (vertical integral) for given profiles
    ! no rain
    implicit none
        integer(4), intent(in)    :: ioxy  ! 1 = RSS, 2 = Rosenkranz 2016
        integer(4), intent(in)    :: nlev        ! number of levels. numered: 0,1,2, ... nlev
        integer(4), intent(in)    :: ilevel0     ! index of 1st level to start integration     
        real(4),    intent(in)    :: freq        ! frequency [ghz]
        

        real(4), dimension(0:nlev), intent(in)    ::    z           ! elevation above sea level         [m]    

        real(4), dimension(0:nlev), intent(in)    ::    p           ! air pressure profile                [mbar]    
        real(4), dimension(0:nlev), intent(in)    ::    t           ! air temperaturee profile            [k]    
        real(4), dimension(0:nlev), intent(in)  ::    pv          ! water vapor pressure profile        [mbar]
        real(4), dimension(0:nlev), intent(in)    ::    rhol        ! liquid cloud water density profile[g/m**3]

        
        real(4), intent(out)                    ::    ao            ! total vertical o2          absorption [neper]    
        real(4), intent(out)                    ::    av            ! total vertical water vapor absorption [neper]    
        real(4), intent(out)                    ::    al            ! total liquid cloud         absorption [neper]    
        real(4), intent(out)                    ::    opacity        ! total                      absorption [neper]    

        real(4), dimension(0:nlev)        ::  alpha_v, alpha_o, alpha_l
        integer(4)                        ::    ipr, ierr

    !   profile of absorption coefficients
        alpha_o=0.0
        alpha_v=0.0
        alpha_l=0.0

        do ipr  = ilevel0,nlev
            call  fdabscoeff(ioxy=ioxy, freq=freq, p=p(ipr), t=t(ipr), pv=pv(ipr), ao=alpha_o(ipr), av=alpha_v(ipr))   ! [ neper/km]
            if (rhol(ipr) > 1.0e-7 .and. t(ipr)>230.) then
                call fdcldabs(freq,t(ipr),rhol(ipr),   alpha_l(ipr))
            else
                alpha_l(ipr) = 0.0
            endif
        enddo
        
        alpha_o = alpha_o/1000.0     ! [neper/m]
        alpha_v = alpha_v/1000.0     ! [neper/m]
        alpha_l = alpha_l/1000.0     ! [neper/m]


        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_o(ilevel0:nlev),2,   ao, ierr)
        if (ierr /= 0) stop ' error in o2  abs profile integration'

        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_v(ilevel0:nlev),2,   av, ierr)
        if (ierr /= 0) stop ' error in vap abs profile integration'

        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_l(ilevel0:nlev),1,   al, ierr)
        if (ierr /= 0) stop ' error in cld abs profile integration'    

        opacity= ao+av+al
            
    return
    end subroutine get_opacity


    subroutine get_opacity_c(ioxy, nlev, ilevel0, freq, z, p, t, pv, rhol,        ao, av_comp, al, opacity)
    ! calculates atmospheric path length components (vertical integral) for given profiles
    ! no rain
    ! same as get_opacity but return separate H2O components (line, FB, SB)
    implicit none
        integer(4), intent(in)    :: ioxy  ! 1 = RSS, 2 = Rosenkranz 2016
        integer(4), intent(in)    :: nlev        ! number of levels. numered: 0,1,2, ... nlev
        integer(4), intent(in)    :: ilevel0     ! index of 1st level to start integration     
        real(4),    intent(in)    :: freq        ! frequency [ghz]
        

        real(4), dimension(0:nlev), intent(in)    ::    z           ! elevation above sea level         [m]    

        real(4), dimension(0:nlev), intent(in)    ::    p           ! air pressure profile                [mbar]    
        real(4), dimension(0:nlev), intent(in)    ::    t           ! air temperaturee profile            [k]    
        real(4), dimension(0:nlev), intent(in)  ::    pv          ! water vapor pressure profile        [mbar]
        real(4), dimension(0:nlev), intent(in)    ::    rhol        ! liquid cloud water density profile[g/m**3]

        
        real(4), intent(out)                    ::    ao            ! total vertical o2          absorption [neper]    
        real(4), dimension(3),intent(out)       ::    av_comp     ! total vertical water vapor absorption [neper] components, 1=line, 2=FB, 3=SB    
        real(4), intent(out)                    ::    al            ! total liquid cloud         absorption [neper]    
        real(4), intent(out)                    ::    opacity        ! total                      absorption [neper]    

        real(4), dimension(0:nlev)        ::  alpha_o, alpha_l, alpha_v
        real(4), dimension(3,0:nlev)    ::  alpha_v_c
        integer(4)                        ::    ipr, ierr, icomp

    !   profile of absorption coefficients
        alpha_o=0.0
        alpha_v=0.0
        alpha_l=0.0
        alpha_v_c=0.0


        do ipr  = ilevel0,nlev
            call  fdabscoeff(ioxy=ioxy, freq=freq, p=p(ipr), t=t(ipr), &
                             pv=pv(ipr), ao=alpha_o(ipr), av=alpha_v(ipr), &
                             av_c=alpha_v_c(1:3,ipr))   ! [ neper/km]
            if (rhol(ipr) > 1.0e-7 .and. t(ipr)>230.) then
                call fdcldabs(freq,t(ipr),rhol(ipr),   alpha_l(ipr))
            else
                alpha_l(ipr) = 0.0
            endif
        enddo
        
        alpha_o = alpha_o/1000.0         ! [neper/m]
        alpha_v = alpha_v/1000.0         ! [neper/m]
        alpha_l = alpha_l/1000.0         ! [neper/m]
        alpha_v_c = alpha_v_c/1000.0     ! [neper/m]


        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_o(ilevel0:nlev),2,   ao, ierr)
        if (ierr /= 0) stop ' error in o2  abs profile integration'    

        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_l(ilevel0:nlev),1,   al, ierr)
        if (ierr /= 0) stop ' error in cld abs profile integration'    
        
        do icomp=1,3
        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_v_c(icomp,ilevel0:nlev),2,   av_comp(icomp), ierr)
        if (ierr /= 0) stop ' error in vap abs profile integration'
        enddo

        opacity= ao+av_comp(1)+av_comp(2)+av_comp(3)+al
            
    return
    end subroutine get_opacity_c


    subroutine get_teff(ioxy, nlev, ilevel0, freq, z, p, t, pv, rhol, teff_o, teff_v, teff_l)
    ! calculates effective temperature of O2, vap, cld columns 
    ! vertical integrals weighted by absorption for given profiles
    implicit none
        integer(4), intent(in)    :: ioxy  ! 1 = RSS, 2 = Rosenkranz 2016
        integer(4), intent(in)    :: nlev        ! number of levels. numered: 0,1,2, ... nlev
        integer(4), intent(in)    :: ilevel0     ! index of 1st level to start integration     
        real(4),    intent(in)    :: freq        ! frequency [ghz]
        
        real(4), dimension(0:nlev), intent(in)    ::    z           ! elevation above sea level         [m]    

        real(4), dimension(0:nlev), intent(in)    ::    p           ! air pressure profile                [mbar]    
        real(4), dimension(0:nlev), intent(in)    ::    t           ! air temperaturee profile            [k]    
        real(4), dimension(0:nlev), intent(in)  ::    pv          ! water vapor pressure profile        [mbar]
        real(4), dimension(0:nlev), intent(in)    ::    rhol        ! liquid cloud water density profile[g/m**3]

        
        real(4), intent(out)                    ::    teff_o       ! effective O2 temperature           [Kelvin]         
        real(4), intent(out)                    ::    teff_v        ! effective water vapor temperature  [Kelvin]
        real(4), intent(out)                    ::    teff_l       ! effective cloud temperature        [Kelvin]

        real(4), dimension(0:nlev)        ::  alpha_v, alpha_o, alpha_l
        integer(4)                        ::    ipr, ierr
        
        real(4)                         ::  ao, av, al
        real(4)                         ::  yo, yv, yl

    !   profile of absorption coefficients
        alpha_o=0.0
        alpha_v=0.0
        alpha_l=0.0

        do ipr  = ilevel0,nlev
            call  fdabscoeff(ioxy=ioxy, freq=freq, p=p(ipr), t=t(ipr), pv=pv(ipr), ao=alpha_o(ipr), av=alpha_v(ipr))   ! [ neper/km]
            if (rhol(ipr) > 1.0e-7 .and. t(ipr)>230.) then
                call fdcldabs(freq,t(ipr),rhol(ipr),   alpha_l(ipr))
            else
                alpha_l(ipr) = 0.0
            endif
        enddo
        
        alpha_o = alpha_o/1000.0     ! [neper/m]
        alpha_v = alpha_v/1000.0     ! [neper/m]
        alpha_l = alpha_l/1000.0     ! [neper/m]


        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_o(ilevel0:nlev),2,   ao, ierr)
        if (ierr /= 0) stop ' error in o2  abs profile integration'

        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_v(ilevel0:nlev),2,   av, ierr)
        if (ierr /= 0) stop ' error in vap abs profile integration'

        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_l(ilevel0:nlev),1,   al, ierr)
        if (ierr /= 0) stop ' error in cld abs profile integration'    

        ! weighted temperature average
        ! absorption weighting     
        
        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_o(ilevel0:nlev)*t(ilevel0:nlev),2,   yo, ierr)
        if (ierr /= 0) stop ' error in o2  abs profile integration'

        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_v(ilevel0:nlev)*t(ilevel0:nlev),2,   yv, ierr)
        if (ierr /= 0) stop ' error in vap abs profile integration'
        
        call column_integral (nlev-ilevel0, z(ilevel0:nlev), alpha_l(ilevel0:nlev)*t(ilevel0:nlev),1,   yl, ierr)
        if (ierr /= 0) stop ' error in cld abs profile integration'        

        teff_o = yo/ao
        teff_v = yv/av
        teff_l = yl/al
        
            
    return
    end subroutine get_teff



    subroutine atm_tran(nlev,tht,t,z,tabs,  tran,tbdw,tbup)
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
    !             tabs(0:nlev)   atmosphric absorptrion coefficients [nepers/m]
    !             t(0:nlev)      temperature profile[in k]
    !              z(0:nlev)      elevation (m) 


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
            stop ' error in profile order in routine atm_tran, pgm stopped.'    
        endif
        h = 0.5*(z(i)+z(i-1))/1000.  ! average elevation of profile layer in km
        delta=h/re
        if (delta <= 0) stop ' error in computing delta in routine atm_tran. pgm stopped.'
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



    subroutine fdabscoeff(ioxy,freq,p,t,pv, ao,av,av_c)
    use RSS_2022_ABSORPTION

    !     freq  frequency [in ghz]
    !     p     pressure [in h pa]
    !     t     temperature [in k]
    !     pv    water vapor pressure  [in hpa]
    !
    !     output:
    !     ao          dry air absortption coefficient      [neper/km]    
    !     av          water vapor absorption coefficients [neper/km]
    !     av_c        water vapor absorption coefficients [neper/km] components 1=line, 2=FBC 2=SBC


    implicit none

        integer(4), intent(in)           :: ioxy  ! 1 = RSS, 2 = Rosenkranz 2016
        real(4), intent(in)              :: freq
        real(4), intent(in)              :: p,t,pv

        real(4), intent(out)             :: ao
        real(4), intent(out)             :: av
        real(4), dimension(3),optional, intent(out) :: av_c
        
        real(4), dimension(3)          :: xav_c
        
        if (ioxy == 1) then
            call absoxy_RSS_2022(p=p,t=t,pv=pv,freq=freq,  absoxy=ao) 
        else if (ioxy == 2) then
            ao = O2ABS(t,p,pv,freq)
        else
            stop 'ioxy must be 1 or 2'
        endif
        call absh2o_RSS_2022(p=p,t=t,pv=pv,freq=freq,  absh2o=av, absh2o_comp=xav_c)
        if (present(av_c)) av_c=xav_c
    
    return
    end subroutine fdabscoeff



    subroutine fdcldabs(freq,t,rhol,   al)
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
    end    subroutine find_permittivity_meissner_wentz    


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



end module rtm_atmosphere