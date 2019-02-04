#!/usr/bin/env python


'''
                    UNLV

Synopsis:  

This routine carries out a series of thin shell python simulations.
The wind mdot is set to produce a hydrogen density of 1e7 - to change
this one had to change the wind_mdot parameter. The loop is carried out
over the 2-10kev luminosity of the central source. A luminosity of
1e21 give a mainly neutral plasma (IP=5e-10) whilst 1e37 is almost 
totally ionized (IP=5e6). The code strips out ion densities and constructs 
data files with the cloudy type ionization parameter as the first 
line.
If one wants a constant temperature file, then one needs to hack python
to stop the temperature changing
We currently write out H,He,C,N,O and Fe. We also output heating and cooling
mechanisms.


Command line usage (if any):

	usage: python_pl_loop *optional file suffix*

Description:  

Primary routines:

Notes:
									   
History:

081214 nsh Commented

'''


if __name__ == "__main__":		# allows one to run from command line without running automatically with write_docs.py

	import sys, subprocess
	import numpy as np

	#Use an optional suffix so file will be py_hydrogen_*suffix*.dat, 
	#If nothing is supplied, default is PL.
	
	alpha=-0.9
	npoints=101 			#this is the parameter which says how many runs to perform - we do 10 per order of magnitude in lum
	python_ver="py83"   	#This is where we define the version of python to use
	atomic="data/standard80"
	python_opts=" "
	ncycles=20
	nphot=100000
	t_e=100000
	nh=1e7
	lum_start=1e25
	nprocs=4
	wind_mdot=4.72694719429145E-20*(nh/1e7)


	if  len(sys.argv) > 1:
		run_name=sys.argv[1]
		inp=open(run_name+".param")
		for line in inp.readlines():
			data=line.split()
			if len(data)!=2:
				print ("Improperly structured param file - each line should be a param name and value seperated by a space")
			else:
				print data[0],data[1]
				if data[0]=='alpha': 
					alpha=data[1]
				if data[0]=='nprocs': 
					nprocs=int(data[1])
				if data[0]=='npoints':
					npoints=int(data[1])
				if data[0]=='t_e':
					t_e=float(data[1])
				if data[0]=='python_ver':
					python_ver=data[1]
				if data[0]=='atomic':
					atomic=data[1]
				if data[0]=='python_opts':
					python_opts=data[1]
				if data[0]=='ncycles':
					ncycles=int(data[1])
				if data[0]=='nphot':
					nphot=data[1]
				if data[0]=='nh':
					nh=data[1]
				if data[0]=='lum_start':
					lum_start=float(data[1])
	else:
		run_name='PL'


		
	

	#Open empty files to contain data

	Hout=open('py_hydrogen_'+run_name+'.dat','w')
	Hout.write('U xi T_e H1 H2\n')
	Heout=open('py_helium_'+run_name+'.dat','w')
	Heout.write('U xi T_e He1 He2 He3\n')
	Cout=open('py_carbon_'+run_name+'.dat','w')
	Cout.write('U xi T_e C1 C2 C3 C4 C5 C6 C7\n')
	Nout=open('py_nitrogen_'+run_name+'.dat','w')
	Nout.write('U xi T_e N1 N2 N3 N4 N5 N6 N7 N8\n')
	Oout=open('py_oxygen_'+run_name+'.dat','w')
	Oout.write('U xi T_e O1 O2 O3 O4 O5 O6 O7 O8 O9\n')
	Feout=open('py_iron_'+run_name+'.dat','w')
	Feout.write('IU xi T_e')
	for i in range(27):
			Feout.write(" Fe"+str(i+1))
	Feout.write('\n')
	Tout=open('py_temperature_'+run_name+'.dat','w')
	Tout.write("U xi T_e Lum N_e heat cool convergence\n")
	Heat=open('py_heat_'+run_name+'.dat','w')
	Heat.write('U xi T_e total photo ff compton ind_comp lines auger\n')
	Cool=open('py_cool_'+run_name+'.dat','w')
	Cool.write('U xi T_e total recomb ff compton DR DI lines Adiabatic\n')



	for i in range(npoints):
		lum=10**((float(i)+np.log10(lum_start)*10.0)/10.0)   #The 210 means the first luminosity is 21.0
		print 'Starting cycle '+str(i+1)+' of '+str(npoints)
		print 'Lum= '+str(lum)
		inp =open('input.pf','w')
		inp.write("System_type(star,binary,agn,previous)     agn\n")
		
		inp.write("### Parameters for the Central Object\n")
		inp.write("Central_object.mass(msol)                  0.8\n")
		inp.write("Central_object.radius(cm)                  1e10\n")
		
		inp.write("### Parameters for the Disk (if there is one)\n")
		inp.write("Disk.type(none,flat,vertically.extended)                 none\n")

		inp.write("### Parameters for BL or AGN\n")
		inp.write("BH.radiation(yes,no)                            yes\n")
		inp.write("BH.rad_type_to_make_wind(bb,models,power,cloudy,brems)               cloudy\n")
		inp.write("BH.lum(ergs/s)            "+str(lum)+"\n")
		inp.write("BH.power_law_index                   "+str(alpha)+"\n")
		inp.write("BH.geometry_for_pl_source(sphere,lamp_post)               sphere\n")
		inp.write("low_energy_break(ev)                       0.136\n")
		inp.write("high_energy_break(ev)                      20000\n")

		inp.write("### Parameters descibing the various winds or coronae in the system\n")
		inp.write("Wind.radiation(yes,no)                           no\n")
		inp.write("Wind.number_of_components                  1\n")
		inp.write("Wind.type(SV,star,hydro,corona,kwd,homologous,yso,shell,imported)                shell\n")
		inp.write("Wind.coord_system(spherical,cylindrical,polar,cyl_var)            spherical\n")
		inp.write("Wind.dim.in.x_or_r.direction               4\n")

		inp.write("### Parameters associated with photon number, cycles,ionization and radiative transfer options\n")
		inp.write("Photons_per_cycle                     "+str(nphot)+"\n")
		inp.write("Ionization_cycles                     "+str(ncycles)+"\n")
		inp.write("Spectrum_cycles                       0\n")
		inp.write("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)           matrix_pow\n")
		inp.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)     thermal_trapping\n")
		inp.write("Atomic_data                                "+atomic+"\n")
		inp.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)               absorb\n")
		inp.write("Thermal_balance_options(0=everything.on,1=no.adiabatic)   0\n")

		inp.write("### Parameters for Domain 0\n")
		inp.write("Shell.wind_mdot(msol/yr)                   "+str(wind_mdot)+"\n")
		inp.write("Shell.wind.radmin(cm)                      1e11\n")
		inp.write("Shell.wind.radmax(cm)                      1.00000000001e11\n")
		inp.write("Shell.wind_v_at_rmin(cm)                   1.00000\n")
		inp.write("Shell.wind.v_at_rmax(cm)                   1.000010\n")
		inp.write("Shell.wind.acceleration_exponent           1\n")
		inp.write("Wind.t.init                                "+str(t_e)+"\n")
		inp.write("Wind.filling_factor(1=smooth,<1=clumped)   1\n")

		inp.write("### Parameters for Reverberation Modeling (if needed)\n")
		inp.write("Reverb.type(none,photon,wind,matom)                    none\n")

		inp.write("### Other parameters\n")
		inp.write("Photon_sampling.approach(T_star,cv,yso,AGN,min_max_freq,user_bands,cloudy_test,wide,logarithmic)   cloudy_test\n")
		inp.write("Photon_sampling.low_energy_limit(eV)       0.0001\n")
		inp.write("Photon_sampling.high_energy_limit(eV)      100000000\n")
		inp.close()
		
		if nprocs==0:
			cmd="time "+python_ver+" "+python_opts+" input > output"
		else:
			cmd="time mpirun -n "+str(nprocs)+" "+python_ver+" "+python_opts+" input > output"
		print cmd
		subprocess.check_call(cmd,shell=True)	   #This is where we define the version of python to use
		subprocess.check_call("tail -n 80 output  | grep Summary > temp",shell=True)#Strip the last 60 lines from the output
		subprocess.check_call("tail -n 80 output  | grep OUTPUT >> temp",shell=True)#Strip the last 60 lines from the output
		inp=open('temp','r')
		for line in inp.readlines():
			data=line.split()
			if (data[0]=='Summary' and data[1]=='convergence'):
				if data[2]=='1':
					print "Simulation converged"
					conv='1'
				else:
					print "Simulation did not converge"
					conv='0'
			if (data[1]=='Lum_agn='):       #This marks the start of the data we want. Field 14 is the cloudy ionization parameter
				Hout.write(data[14]+' '+data[16]+' '+data[4])
				Heout.write(data[14]+' '+data[16]+' '+data[4])
				Cout.write(data[14]+' '+data[16]+' '+data[4])
				Nout.write(data[14]+' '+data[16]+' '+data[4])
				Oout.write(data[14]+' '+data[16]+' '+data[4])
				Feout.write(data[14]+' '+data[16]+' '+data[4])
				Heat.write(data[14]+' '+data[16]+' '+data[4])
				Cool.write(data[14]+' '+data[16]+' '+data[4])
				Tout.write(data[14]+' '+data[16]+' '+data[4]+' '+data[2]+' '+data[8]+' ')
			if (data[1]=='Absorbed_flux(ergs-1cm-3)'):
				Heat.write(' '+data[2]+' '+data[4]+' '+data[6]+' '+data[8]+' '+data[10]+' '+data[12]+' '+data[14]+'\n')
				Tout.write(data[2]+' ')
			if (data[1]=='Wind_cooling(ergs-1cm-3)'):
				Cool.write(' '+data[2]+' '+data[4]+' '+data[6]+' '+data[8]+' '+data[10]+' '+data[12]+' '+data[16]+' '+data[14]+'\n')
				Tout.write(data[2]+' '+conv+'\n')
			if (data[1]=='H'):
				for j in range(2):
					Hout.write(' '+data[j+2])
				Hout.write('\n')
			if (data[1]=='He'):
				for j in range(3):
					Heout.write(' '+data[j+2])
				Heout.write('\n')
			if (data[1]=='C'):
				for j in range(7):
					Cout.write(' '+data[j+2])
				Cout.write('\n')
			if (data[1]=='N'):
				for j in range(8):
					Nout.write(' '+data[j+2])
				Nout.write('\n')
			if (data[1]=='O'):
				for j in range(9):
					Oout.write(' '+data[j+2])
				Oout.write('\n')
			if (data[1]=='Fe'):
				for j in range(27):
					Feout.write(' '+data[j+2])
				Feout.write('\n')
	#Flush the output files so one can see progress and if the loop crashes all is not lost
		Hout.flush()
		Heout.flush()
		Cout.flush()
		Nout.flush()
		Oout.flush()
		Feout.flush()
		Tout.flush()
		Heat.flush()
		Cool.flush()	
		print 'Finished cycle '+str(i+1)+' of '+str(npoints)
	#Close the files.
	Hout.close()
	Heout.close()
	Cout.close()
	Nout.close()
	Oout.close()
	Feout.close()
	Tout.close()
	Heat.close()
	Cool.close()

