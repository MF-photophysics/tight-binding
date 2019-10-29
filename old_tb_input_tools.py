import numpy as np
import matplotlib.pyplot as plt 
import lammps_structure_tools as structure_tools
import lammps_analysis_tools as analysis_tools
import subprocess
from index_conversions import *
"""
def stoich_to_unit_index(nx,ny,nz,N = 8): 
	# takes unit cell indices nx,ny,nz for a Pb atom and orthorhombic box size N and gives the Pb orthorhombic handle
	cell = [i%N for i in [(nx-nz)//2, ny//2,(nx + nz)//2]]
	Pb_num = 1 + 2*abs(nx - nz)%4 + ny%2 
	return [str(cell[0]) + ":" + str(cell[1]) + ":" + str(cell[2]), "Pb" + str(Pb_num)]


def I_x_ortho(Pb_index, N = 8): 
	Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	if Pb_num == 1: 
		I_cell,I_num = Pb_cell, 2
	elif Pb_num == 2: 
		I_cell,I_num = Pb_cell, 3
	elif Pb_num == 3: 
		I_cell,I_num = [Pb_cell[0],(Pb_cell[1] - 1)%N, Pb_cell[2]], 10
	elif Pb_num == 4: 
		I_cell,I_num = Pb_cell, 11
	return [str(I_cell[0]) + ":" + str(I_cell[1]) + ":" + str(I_cell[2]), "I" + str(I_num)]

def I_y_ortho(Pb_index, N = 8):
	Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	if Pb_num == 1: 
		I_cell,I_num = [Pb_cell[0],Pb_cell[1],(Pb_cell[2]-1)% N], 1
	elif Pb_num == 2: 
		I_cell,I_num = [(Pb_cell[0]-1)%N, Pb_cell[1], Pb_cell[2]], 12
	elif Pb_num == 3: 
		I_cell,I_num = Pb_cell, 7
	elif Pb_num == 4: 
		I_cell,I_num = Pb_cell, 6
	return [str(I_cell[0]) + ":" + str(I_cell[1]) + ":" + str(I_cell[2]), "I" + str(I_num)]

def I_z_ortho(Pb_index, N = 8):
	Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	if Pb_num == 1: 
		I_cell,I_num = [(Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]], 8
	elif Pb_num == 2: 
		I_cell,I_num = [(Pb_cell[0]-1)%N, Pb_cell[1], Pb_cell[2]], 9
	elif Pb_num == 3: 
		I_cell,I_num = [Pb_cell[0],(Pb_cell[1]-1)%N,Pb_cell[2]], 4
	elif Pb_num == 4: 
		I_cell,I_num = Pb_cell, 5
	return [str(I_cell[0]) + ":" + str(I_cell[1]) + ":" + str(I_cell[2]), "I" + str(I_num)]


def periodic_fix(r, box_vect,tilting = None):#, verbose =False): 
	t = np.array([0.0,0.0,0.0]) # translation, so that function is not mutative
	if tilting: 
		xy,xz,yz = tilting
		v1 = np.array([box_vect[0], 0,0])
		v2 = np.array([xy,box_vect[1],0])
		v3 = np.array([xz,yz,box_vect[2]])
	else: 
		v1,v2,v3 = np.array([box_vect[0],0,0]),np.array([0,box_vect[1],0]),np.array([0,0,box_vect[2]])
	vects = [v1,v2,v3]
	for i in range(3): 
		if r[i] > box_vect[i]/2:
			t += -vects[i]
		elif r[i] < -box_vect[i]/2:
			t += vects[i]
	#for i in range(3): 
	#	if r[i] > box_vect[i]/2: 
	#		t[i] =  - box_vect[i]
	#	elif r[i] < -box_vect[i]/2: 
	#		t[i] = box_vect[i]
	#if verbose: 
	#	print(r,t)
	return r + t	
"""

def kcal_to_eV(E): 
	return 0.04343 * E

def write_hopping_data(dump, filename,timesteps = None, ortho_cell = False, n = None, n_vec = None, charges = "old"):
	if timesteps is None: 
		timesteps = dump.timesteps.keys()
	E_PB_S, E_PB_P, E_I_P = -2.5, 6.0, 2.5 # averages, in eV
	SP_SIG_a,SP_SIG_b, SP_SIG_c = 42.49, 0.964, -0.54
	PP_SIG_a, PP_SIG_b, PP_SIG_c = 21.61,3.070, -5.88
	PP_PI_a, PP_PI_b, PP_PI_c =  11.66, 1.030, -0.16
	if charges == "old":
		Q_PB, Q_I = 2.03, -1.13
	elif charges == "new":
		Q_PB, Q_I = 1.4012, -0.7006
	t_sp = lambda r: SP_SIG_a*np.exp(-np.linalg.norm(r)/SP_SIG_b) + SP_SIG_c
	t_pp_sig = lambda r: PP_SIG_a*np.exp(-np.linalg.norm(r)/PP_SIG_b) + PP_SIG_c
	t_pp_pi = lambda r: PP_PI_a*np.exp(-np.linalg.norm(r)/PP_PI_b) + PP_PI_c
	N_atoms = len(dump.structure.atoms)
	
	if n is None and n_vec is None: 
		n = round((N_atoms/12)**(1/3)) # dimnesions in stoichiometric cells
		# ^ don't think that will work
	if n_vec is None:
		nx_max,ny_max,nz_max = n,n,n
	else: 
		nx_max,ny_max,nz_max = n_vec

	if ortho_cell: 
		N = round((N_atoms/48)**(1/3)) # number of unit cells in MD run (4 stoich cells each)
		print(N)
		index_conversion = stoich_to_unit_index
		I_x,I_y,I_z = I_x_ortho, I_y_ortho,I_z_ortho
	elif ortho_cell is False: 
		print(n)
		N = round((N_atoms/12)**(1/3)) #n
		#index_conversion = lambda nx,ny,nz,N = None: [str(nx%N) + ":" + str(ny%N) + ":" + str(nz%N), "Pb"]
		index_conversion = lambda nx,ny,nz,N = None: [(nx%N, ny%N, nz%N), "Pb"]
		I_x,I_y,I_z = lambda Pb_index,N = None : [Pb_index[0], "I1"], lambda Pb_index, N = None : [Pb_index[0], "I2"], lambda Pb_index, N = None : [Pb_index[0], "I3"]
	Pb_coul_sum, I_coul_sum = 0,0 # to compute averages
	for step in timesteps: # this loop computes average coulomb potentials
		for nx in range(nx_max): 
			for ny in range(ny_max): 
				for nz in range(nz_max): 
					Pb_index = index_conversion(nx,ny,nz,N)
					Ix_plus_index, Iy_plus_index, Iz_plus_index = I_x(Pb_index,N), I_y(Pb_index, N), I_z(Pb_index, N)
					Pb_coul_sum += kcal_to_eV(2 * dump.atom_pot(Pb_index, step, id_mode = "handles") / Q_PB)
					I_coul_sum += kcal_to_eV(2 * dump.atom_pot(Ix_plus_index,step,id_mode = "handles")/ Q_I + 2 * dump.atom_pot(Iy_plus_index,step, id_mode = "handles" )/ Q_I +  2 * dump.atom_pot(Iz_plus_index, step, id_mode = "handles")/ Q_I)
	# numbers may need changing
	Pb_coul_ave = Pb_coul_sum / ((nx_max*ny_max*nz_max)* len(timesteps))
	print("Pb ave  " + str(Pb_coul_ave))
	I_coul_ave = I_coul_sum / ((3 * nx_max*ny_max*nz_max) * len(timesteps))
	print("I ave  " + str(I_coul_ave))
	lines = []
	for step in timesteps:
		box = dump.timesteps[step].box_vect
		if hasattr(dump.timesteps[step],'xy'): 
			tilt = [dump.timesteps[step].xy,dump.timesteps[step].xz,dump.timesteps[step].yz]
		else: 
			tilt = None
		for nx in range(nx_max): 
			for ny in range(ny_max): 
				for nz in range(nz_max): 
					Pb_index = index_conversion(nx,ny,nz,N)
					Ix_plus_index, Iy_plus_index, Iz_plus_index = I_x(Pb_index,N), I_y(Pb_index, N), I_z(Pb_index, N)
					# wrong, kept around for now.Ix_minus_index, Iy_minus_index, Iz_minus_index = I_x(index_conversion((nx-1)%n, ny, nz,N),N), I_y(index_conversion(nx, (ny-1)%n, nz,N),N), I_z(index_conversion(nx, ny, (nz-1)%n,N),N)
					Ix_minus_index, Iy_minus_index, Iz_minus_index = I_x(index_conversion((nx-1), ny, nz,N),N), I_y(index_conversion(nx, (ny-1), nz,N),N), I_z(index_conversion(nx, ny, (nz-1),N),N)
					# ^ probably fixed
					coul_Pb = kcal_to_eV(2 * dump.atom_pot(Pb_index, step, id_mode = "handles") / Q_PB)
					coul_Ix, coul_Iy, coul_Iz = kcal_to_eV(2 * dump.atom_pot(Ix_plus_index,step,id_mode = "handles")/ Q_I), kcal_to_eV(2 * dump.atom_pot(Iy_plus_index,step, id_mode = "handles" )/ Q_I), kcal_to_eV(2 * dump.atom_pot(Iz_plus_index, step, id_mode = "handles")/ Q_I)
					
					lines.append(str(E_PB_S - (coul_Pb - Pb_coul_ave)) + (" " + str(E_PB_P - (coul_Pb - Pb_coul_ave)))*3)
					lines.append(((" " + str(E_I_P - (coul_Ix - I_coul_ave)))*4)[1:])
					lines.append(((" " + str(E_I_P - (coul_Iy - I_coul_ave)))*4)[1:])
					lines.append(((" " + str(E_I_P - (coul_Iz - I_coul_ave)))*4)[1:])

					Pb_pos = dump.atom_pos(Pb_index, step, id_mode = "handles")
					Ix_plus_pos, Iy_plus_pos, Iz_plus_pos = dump.atom_pos(Ix_plus_index,step, id_mode = "handles"),dump.atom_pos(Iy_plus_index,step, id_mode = "handles"),dump.atom_pos(Iz_plus_index,step,id_mode = "handles")
					Ix_minus_pos, Iy_minus_pos, Iz_minus_pos = dump.atom_pos(Ix_minus_index,step, id_mode = "handles"),dump.atom_pos(Iy_minus_index,step, id_mode = "handles"),dump.atom_pos(Iz_minus_index,step,id_mode = "handles")
					t_sp_box = lambda r: t_sp(periodic_fix(r,box,tilt))
					t_pp_sig_box = lambda r: t_pp_sig(periodic_fix(r,box,tilt))
					t_pp_pi_box = lambda r: t_pp_pi(periodic_fix(r,box,tilt))
					lines.append(str(t_sp_box(Pb_pos - Ix_plus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Ix_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_plus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iy_plus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iy_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_plus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iz_plus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iz_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_plus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Ix_minus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Ix_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_minus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iy_minus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iy_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_minus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iz_minus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iz_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_minus_pos)))

					
	try:
		f = open(filename, "w")
		f.write("\n".join(lines))# apparently .join is the most efficient way to write large strings
	except MemoryError:
		f = open(filename, "a")
		for line in lines[:-1]:
			f.write(line + "\n")
		f.write(lines[-1])
		

def vdw_coul_potential(E_tot, E_coul, q, n, factor = 2): 
	# q is charge for coulomb, n is number of electrons for van der waals
	E_vdw = E_tot - E_coul
	V_coul = kcal_to_eV(2*E_coul) / q
	V_vdw = factor * kcal_to_eV(2*E_vdw) / n  # if vdw potential is linear in number of electrons, factor is 1
	# it seems like the factor should be 2 based on perturbation theory, because E = k*n^2 so dE/dn = 2E/n
	return V_vdw - V_coul # subtract V_coul because electron is negative


def write_hopping_data_vdw(couldump,totdump, filename,timesteps = None, ortho_cell = False, n = None,n_vec = None, charges = "old"):
	if timesteps is None: 
		timesteps = couldump.timesteps.keys()
	E_PB_S, E_PB_P, E_I_P = -2.5, 6.0, 2.5 # averages, in eV
	SP_SIG_a,SP_SIG_b, SP_SIG_c = 42.49, 0.964, -0.54
	PP_SIG_a, PP_SIG_b, PP_SIG_c = 21.61,3.070, -5.88
	PP_PI_a, PP_PI_b, PP_PI_c =  11.66, 1.030, -0.16
	if charges == "old":
		Q_PB, Q_I = 2.03, -1.13
	elif charges == "new":
		Q_PB, Q_I = 1.4012, -0.7006
	N_PB, N_I = 4 - Q_PB, 5 - Q_I # I am choosing to ignore s electrons on iodine
	t_sp = lambda r: SP_SIG_a*np.exp(-np.linalg.norm(r)/SP_SIG_b) + SP_SIG_c
	t_pp_sig = lambda r: PP_SIG_a*np.exp(-np.linalg.norm(r)/PP_SIG_b) + PP_SIG_c
	t_pp_pi = lambda r: PP_PI_a*np.exp(-np.linalg.norm(r)/PP_PI_b) + PP_PI_c
	N_atoms = len(couldump.structure.atoms)

	if n is None and n_vec is None: 
		n = round((N_atoms/12)**(1/3)) # dimnesions in stoichiometric cells
		# ^ don't think that will work
	if n_vec is None:
		nx_max,ny_max,nz_max = n,n,n
	else: 
		nx_max,ny_max,nz_max = n_vec

	if ortho_cell: 
		N = round((N_atoms/48)**(1/3)) # number of unit cells in MD run (4 stoich cells each)
		print(N)
		index_conversion = stoich_to_unit_index
		I_x,I_y,I_z = I_x_ortho, I_y_ortho,I_z_ortho
	elif ortho_cell is False: 
		print(n)
		N = round((N_atoms/12)**(1/3)) #n
		#index_conversion = lambda nx,ny,nz,N = None: [str(nx%N) + ":" + str(ny%N) + ":" + str(nz%N), "Pb"]
		index_conversion = lambda nx,ny,nz,N = None: [(nx%N, ny%N, nz%N), "Pb"]
		I_x,I_y,I_z = lambda Pb_index,N = None : [Pb_index[0], "I1"], lambda Pb_index, N = None : [Pb_index[0], "I2"], lambda Pb_index, N = None : [Pb_index[0], "I3"]
	Pb_pot_sum, I_pot_sum = 0,0 # to compute averages
	for step in timesteps: # this loop computes average coulomb potentials
		for nx in range(nx_max): 
			for ny in range(ny_max): 
				for nz in range(nz_max): 
					Pb_index = index_conversion(nx,ny,nz,N)
					Ix_plus_index, Iy_plus_index, Iz_plus_index = I_x(Pb_index,N), I_y(Pb_index, N), I_z(Pb_index, N)
					Pb_pot_sum += vdw_coul_potential(totdump.atom_pot(Pb_index, step, id_mode = "handles"), couldump.atom_pot(Pb_index, step, id_mode = "handles"), Q_PB, N_PB)
					I_pot_sum += vdw_coul_potential(totdump.atom_pot(Ix_plus_index,step,id_mode = "handles"), couldump.atom_pot(Ix_plus_index,step,id_mode = "handles"), Q_I, N_I)
					I_pot_sum +=  vdw_coul_potential(totdump.atom_pot(Iy_plus_index,step,id_mode = "handles"), couldump.atom_pot(Iy_plus_index,step,id_mode = "handles"), Q_I, N_I)
					I_pot_sum +=  vdw_coul_potential(totdump.atom_pot(Iz_plus_index,step,id_mode = "handles"), couldump.atom_pot(Iz_plus_index,step,id_mode = "handles"), Q_I, N_I)
	# numbers may need changing
	Pb_pot_ave = Pb_pot_sum / ((nx_max*ny_max*nz_max)* len(timesteps))
	print("Pb ave  " + str(Pb_pot_ave))
	I_pot_ave = I_pot_sum / ((3 * nx_max*ny_max*nz_max) * len(timesteps))
	print("I ave  " + str(I_pot_ave))
	lines = []
	for step in timesteps:
		box = couldump.timesteps[step].box_vect
		if hasattr(couldump.timesteps[step],'xy'): 
			tilt = [couldump.timesteps[step].xy,couldump.timesteps[step].xz,couldump.timesteps[step].yz]
		else: 
			tilt = None
		for nx in range(nx_max): 
			for ny in range(ny_max): 
				for nz in range(nz_max): 
					Pb_index = index_conversion(nx,ny,nz,N)
					Ix_plus_index, Iy_plus_index, Iz_plus_index = I_x(Pb_index,N), I_y(Pb_index, N), I_z(Pb_index, N)
					# wrong, kept around for now.Ix_minus_index, Iy_minus_index, Iz_minus_index = I_x(index_conversion((nx-1)%n, ny, nz,N),N), I_y(index_conversion(nx, (ny-1)%n, nz,N),N), I_z(index_conversion(nx, ny, (nz-1)%n,N),N)
					Ix_minus_index, Iy_minus_index, Iz_minus_index = I_x(index_conversion((nx-1), ny, nz,N),N), I_y(index_conversion(nx, (ny-1), nz,N),N), I_z(index_conversion(nx, ny, (nz-1),N),N)
					# ^ probably fixed

					V_Pb = vdw_coul_potential(totdump.atom_pot(Pb_index, step, id_mode = "handles"), couldump.atom_pot(Pb_index, step, id_mode = "handles"), Q_PB, N_PB)
					V_Ix = vdw_coul_potential(totdump.atom_pot(Ix_plus_index,step,id_mode = "handles"), couldump.atom_pot(Ix_plus_index,step,id_mode = "handles"), Q_I, N_I)
					V_Iy =  vdw_coul_potential(totdump.atom_pot(Iy_plus_index,step,id_mode = "handles"), couldump.atom_pot(Iy_plus_index,step,id_mode = "handles"), Q_I, N_I)
					V_Iz =  vdw_coul_potential(totdump.atom_pot(Iz_plus_index,step,id_mode = "handles"), couldump.atom_pot(Iz_plus_index,step,id_mode = "handles"), Q_I, N_I)
					
					
					lines.append(str(E_PB_S  + V_Pb - Pb_pot_ave) + (" " + str(E_PB_P + V_Pb - Pb_pot_ave))*3)
					lines.append(((" " + str(E_I_P + V_Ix - I_pot_ave))*4)[1:])
					lines.append(((" " + str(E_I_P + V_Iy - I_pot_ave))*4)[1:])
					lines.append(((" " + str(E_I_P + V_Iz - I_pot_ave))*4)[1:])

					Pb_pos = couldump.atom_pos(Pb_index, step, id_mode = "handles")
					Ix_plus_pos, Iy_plus_pos, Iz_plus_pos = couldump.atom_pos(Ix_plus_index,step, id_mode = "handles"),couldump.atom_pos(Iy_plus_index,step, id_mode = "handles"),couldump.atom_pos(Iz_plus_index,step,id_mode = "handles")
					Ix_minus_pos, Iy_minus_pos, Iz_minus_pos = couldump.atom_pos(Ix_minus_index,step, id_mode = "handles"),couldump.atom_pos(Iy_minus_index,step, id_mode = "handles"),couldump.atom_pos(Iz_minus_index,step,id_mode = "handles")
					t_sp_box = lambda r: t_sp(periodic_fix(r,box,tilt))
					t_pp_sig_box = lambda r: t_pp_sig(periodic_fix(r,box,tilt))
					t_pp_pi_box = lambda r: t_pp_pi(periodic_fix(r,box,tilt))
					lines.append(str(t_sp_box(Pb_pos - Ix_plus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Ix_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_plus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iy_plus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iy_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_plus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iz_plus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iz_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_plus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_plus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Ix_minus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Ix_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Ix_minus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iy_minus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iy_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iy_minus_pos)))
					lines.append(str(t_sp_box(Pb_pos - Iz_minus_pos)) + " " + str(t_pp_sig_box(Pb_pos - Iz_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_minus_pos)) + " " + str(t_pp_pi_box(Pb_pos - Iz_minus_pos)))

					
	f = open(filename, "w")
	f.write("\n".join(lines))# apparently .join is the most efficient way to write large strings


