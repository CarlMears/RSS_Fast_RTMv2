program make_absorption_table

	use nan_support
	use rtm
	use msu_constants
	use amsu_constants

	implicit none

	integer(4)			:: num_T 
	integer(4)			:: num_p 
	integer(4)			:: num_q 

	real(4)				:: T0 
	real(4)				:: Delta_T 
	real(4)				:: Delta_P 
	real(4)				:: Delta_q 

	real(4),dimension(0:200,0:110,0:150)  :: abs_table
	real(4),dimension(0:200,0:110,0:150)  :: abs_table_per_Pa

	integer(4)						:: T_index,P_index,q_index
	real(4)							:: T,P,PV,q

	integer(4)						:: ivap = 4
	integer(4)						:: ioxy = 1

	real(4)							:: ao,av,total_abs,freq

	integer(4)						:: channel
	integer(4),parameter			:: amsu_num_band = 2
	integer(4),parameter			:: amsu_num_freq = 13
	integer(4)						:: freq_index,band_index
	real(4),dimension(amsu_num_freq,amsu_num_band) :: amsu_freq_arr

	real(8),parameter				:: M_W_air = 2.8966D-2
	real(8),parameter				:: M_W_H2O  = 1.8015324d-2    ! kg/mol
	real(8),parameter				:: R_GAS    = 8.3145112d0
	real(8),parameter				:: g        = 9.80665

	real(8)							:: c_air
	real(8)							:: c_h2o
	real(8)							:: rho_dry
	real(8)							:: rho_vap

		
	character(len=120)		:: file

	c_air = R_GAS/M_W_AIR
	c_h2o = R_GAS/M_W_H2O

	channel = 10

	! construct an array of frequencies where the calc is performed
	do band_index = 0,amsu_num_band-1 
		do freq_index = 0,amsu_num_freq-1
			amsu_freq_arr(freq_index+1,band_index+1) =						&
					AMSU_A_Freq(channel) -									&
					2.0*(band_index-0.5)*AMSU_A_Freq_Split_1(channel) - 	&			
					AMSU_A_BANDWIDTH(channel)/2.0 +							&
					(1+2*freq_index)*AMSU_A_BANDWIDTH(channel)/(2.0*amsu_num_freq)
	    enddo
	enddo

    print *,amsu_freq_arr


	num_T = 200
	num_p = 110
	num_q = 150

	T0 = 140.0
	Delta_T = 1.0
	Delta_P = 10
	Delta_q = 0.001

	
    abs_table = 0.0
	abs_table_per_Pa = 0.0
	do band_index = 1,amsu_num_band
	do freq_index = 1,amsu_num_freq
		freq = amsu_freq_arr(freq_index,band_index)
		print*,freq	
		do T_index = 0,num_T
		   T = T0+Delta_T*T_index
		   do P_index = 0,num_p
				P = Delta_P*P_index
				do q_index = 0,num_q
					q = Delta_q*q_index
					total_abs = 0.0

					PV = p*q/(0.622+0.378*q)

					call FDABSCOEFF(ivap,ioxy,freq,P,T,PV, av,ao)
					total_abs = av + ao
					abs_table(T_index,P_index,q_index) = abs_table(T_index,P_index,q_index) + total_abs
					
					! find density

					rho_dry = (p-pv)/(c_air*T)
					rho_vap = pv/(c_h2o*T)
					abs_table_per_Pa(T_index,P_index,q_index) = abs_table_per_Pa(T_index,P_index,q_index) + &
																total_abs*0.001/((rho_dry + rho_vap) * g)
				
				enddo
			enddo
		enddo
	enddo
	enddo

	abs_table = abs_table/(amsu_num_freq*amsu_num_band)
	abs_table_per_Pa = abs_table_per_Pa/(amsu_num_freq*amsu_num_band)

	file = 'E:\AMSU_MSU_simulation\abs_tables\amsu_10_abs_table_q.dat'
	open(unit= 15,file = file,form = 'BINARY')
	write(15)num_t,num_p,num_q,T0,Delta_T,Delta_P,Delta_q,ivap,ioxy,channel
	write(15)abs_table
	close(15)

	file = 'E:\AMSU_MSU_simulation\abs_tables\amsu_10_abs_table_q_per_Pa.dat'
	open(unit= 15,file = file,form = 'BINARY')
	write(15)num_t,num_p,num_q,T0,Delta_T,Delta_P,Delta_q,ivap,ioxy,channel
	write(15)abs_table_per_Pa
	close(15)


end program