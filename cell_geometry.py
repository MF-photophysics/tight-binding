import numpy as np
import matplotlib.pyplot as plt 
import lammps_structure_tools as structure_tools
import lammps_analysis_tools as analysis_tools
from index_conversions import *

def kcal_to_eV(E): 
	return 0.04343 * E

def angle(v1,v2):
	return np.arccos(np.dot(v1,v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def r_angle(r1,r2,r3,box_vect, tilting): 
	return angle(periodic_fix(r1 - r2, box_vect, tilting), periodic_fix(r3-r2, box_vect,tilting))

def r_distance(r1,r2,box_vect, tilting): 
	return np.linalg.norm(periodic_fix(r1-r2, box_vect, tilting))

def r_projection(r1,r2, v, box_vect, tilting):
	r = periodic_fix(r1-r2, box_vect,tilting)
	return np.dot(r,v) / np.linalg.norm(r)

"""def atom_distance(handle1,handle2,dump, timestep,tilt):
	r1,r2 = dump.atom_pos(handle1, timestep, id_mode = "handles"),dump.atom_pos(handle2, timestep, id_mode = "handles")
	return periodic_fix(r1-r2, dump.timesteps[timestep].box_vect, tilt)
"""

def write_param_file(dump, timesteps, n_vec, ortho_cell = False, N_ortho = None,param_sets = {},feature_sets = {}, coul = False, vdw = False, coul_dump = None,hessian = False):
	#param_sets = {"filename":[param1,param2, ...]}
	lines = {file:[] for file in param_sets}
	for step in timesteps:
		print(step)
		atom_positions = stoich_atom_positions(dump,step, n_vec, ortho_cell,N_ortho)
		box = dump.timesteps[step].box_vect
		if hasattr(dump.timesteps[step],'xy'): 
			tilt = [dump.timesteps[step].xy,dump.timesteps[step].xz,dump.timesteps[step].yz]
		else: 
			tilt = None
		if coul or vdw: 
			#print("potentials")
			atom_potentials = stoich_atom_positions(dump, step, n_vec, ortho_cell, N_ortho,mode = "potentials", coul_dump = coul_dump)
		else: 
			atom_potentials = None
		if hessian: 
			atom_hessians = stoich_atom_positions(dump, step, n_vec, ortho_cell, N_ortho,mode = "hessians",positions = atom_positions)
		else: 
			atom_hessians = None
		for file in param_sets:
			try:
				features = feature_sets[file]
			except KeyError:
				features = False
			parameters = param_sets[file]
			for ix in range(n_vec[0]):
				for iy in range(n_vec[1]):
					for iz in range(n_vec[2]):
						for param in parameters:
							#lines[file].append(" ".join([str(i) for i in param(np.array([ix,iy,iz]), n_vec, atom_positions,box,tilt, coul = coul, vdw = vdw, atom_potentials = atom_potentials)]))
							if features:
								try: 
									lines[file].append(" ".join([str(i) for i in param(np.array([ix,iy,iz]), n_vec, atom_positions,box,tilt, features = features, coul = coul, vdw = vdw, atom_potentials = atom_potentials,hessian = hessian, atom_hessians = atom_hessians)]))
								except TypeError or NameError:
									#"excepted"
									lines[file].append(" ".join([str(i) for i in param(np.array([ix,iy,iz]), n_vec, atom_positions,box,tilt,features = features)]))
							else: 
								try: 
									lines[file].append(" ".join([str(i) for i in param(np.array([ix,iy,iz]), n_vec, atom_positions,box,tilt,  coul = coul, vdw = vdw, atom_potentials = atom_potentials,hessian = hessian, atom_hessians = atom_hessians)]))
								except TypeError or NameError:
									#"excepted"
									lines[file].append(" ".join([str(i) for i in param(np.array([ix,iy,iz]), n_vec, atom_positions,box,tilt)]))
	for filename in param_sets:
		f = open(filename, "w")
		f.write("\n".join(lines[filename]))# apparently .join is the most efficient way to write large strings



def stoich_atom_positions(dump, timestep, n_vec, ortho_cell = False, N_ortho = None, include_hydrogen = False, mode = "positions",positions = None, coul_dump = None): # add in potentials functionality? 
	#probably should be in analysis_tools or index_conversions
	nx,ny,nz = n_vec
	if ortho_cell: 
		index_conversion = stoich_to_unit_index
		Pb = lambda Pb_index, N = None: Pb_index
		I_x,I_y,I_z = I_x_ortho, I_y_ortho,I_z_ortho
		C,N = C_ortho, N_ortho
		Hc1, Hc2, Hc3 = Hc1_ortho, Hc2_ortho, Hc3_ortho
		Hn1, Hn2, Hn3 = Hn1_ortho, Hn2_ortho, Hn3_ortho
		N = round((nx * ny * nz // 4)**(1/3)) if N_ortho is None else N_ortho
	else: 
		N = None
		index_conversion = lambda nx,ny,nz,N = None: [(nx, ny, nz), "Pb"]
		Pb = lambda Pb_index, N = None : Pb_index
		I_x,I_y,I_z = lambda Pb_index,N = None : [Pb_index[0], "I1"], lambda Pb_index, N = None : [Pb_index[0], "I2"], lambda Pb_index, N = None : [Pb_index[0], "I3"]
		C, N = lambda Pb_index,N = None : [Pb_index[0], "MA", "C"], lambda Pb_index,N = None : [Pb_index[0], "MA", "N"]
		Hc1, Hc2,Hc3 = 	lambda Pb_index,N = None : [Pb_index[0], "MA", "H_C1"],lambda Pb_index,N = None : [Pb_index[0], "MA", "H_C2"],lambda Pb_index,N = None : [Pb_index[0], "MA", "H_C3"]
		Hn1, Hn2,Hn3 = 	lambda Pb_index,N = None : [Pb_index[0], "MA", "H_N1"],lambda Pb_index,N = None : [Pb_index[0], "MA", "H_N2"],lambda Pb_index,N = None : [Pb_index[0], "MA", "H_N3"]
	if include_hydrogen: 
		atom_funcs = [Pb, I_x, I_y, I_z, C, N, Hc1, Hc2,Hc3, Hn1,Hn2, Hn3]
	else: 
		atom_funcs = [Pb, I_x, I_y, I_z, C, N]
	def cell_positions(ix,iy,iz, atom_funcs = atom_funcs): 
		return [dump.atom_pos(handles, timestep, id_mode = "handles") for handles in [atom_func(index_conversion(ix,iy,iz,N)) for atom_func in atom_funcs]]
	def cell_potentials(ix,iy,iz,atom_funcs = atom_funcs[:4], charges = [2.03, -1.13,-1.13,-1.13]):
		totpots = [dump.atom_pot(handles, timestep, id_mode = "handles") for handles in [atom_func(index_conversion(ix,iy,iz,N)) for atom_func in atom_funcs]]
		coulpots = [coul_dump.atom_pot(handles, timestep, id_mode = "handles") for handles in [atom_func(index_conversion(ix,iy,iz,N)) for atom_func in atom_funcs]]
		return [(kcal_to_eV(coulpots[i] / charges[i]) , kcal_to_eV(totpots[i] - coulpots[i])) for i in range(len(atom_funcs))]
	def cell_hessians(ix,iy,iz,atom_positions):
		step = dump.timesteps[timestep]
		box,tilt = step.box_vect, [step.xy,step.xz,step.yz] if hasattr(step,"xz") else None
		return [pbhess(np.array((ix,iy,iz)), n_vec, atom_positions, box, tilt), ihess(np.array((ix,iy,iz)), n_vec, 'x',atom_positions, box, tilt), ihess(np.array((ix,iy,iz)), n_vec, 'y',atom_positions, box, tilt), ihess(np.array((ix,iy,iz)), n_vec, 'z',atom_positions, box, tilt)]
	if mode == "positions": 
		return [[[cell_positions(ix,iy,iz) for iz in range(nz)] for iy in range(ny)] for ix in range(nx)]
	elif mode == "potentials": 
		return [[[cell_potentials(ix,iy,iz) for iz in range(nz)] for iy in range(ny)] for ix in range(nx)]
	elif mode == "hessians":
		return [[[cell_hessians(ix,iy,iz,positions) for iz in range(nz)] for iy in range(ny)] for ix in range(nx)]



def direction_info(direction, pdir = None): 
	if direction == "x" or "s": 
		bond_trans = np.array([1,0,0])
		side_trans1 = np.array([0,1,0])
		sdir1 = "y"
		side_trans2 = np.array([0,0,1])
		sdir2 = "z"
	if direction == "y": 
		bond_trans = np.array([0,1,0])
		side_trans1 = np.array([0,0,1])
		sdir1 = "z"
		side_trans2 = np.array([1,0,0])
		sdir2 = "x"
	if direction == "z": 
		bond_trans = np.array([0,0,1])
		side_trans1 = np.array([1,0,0])
		sdir1 = "x"
		side_trans2 = np.array([0,1,0])
		sdir2 = "y"	
	if direction == "x" and pdir == "z" or direction == "y" and pdir == "x" or direction == "z" and pdir == "y":
		side_trans1, side_trans2, sdir1, sdir2 = side_trans2, side_trans1, sdir2, sdir1
		# flip which are which for non-cyclic permutations
	return bond_trans,side_trans1,sdir1, side_trans2,sdir2


def Pb_I_bond_geometry(i_vec,n_vec,direction, atom_positions,box, tilt, pdir = None,negdir = False,fitting = False, features = {"bond_len" : True, "neighbor_len" : True, "bond_angle" : True, "neighbor_angle": True},split_len = 1): 
	bond_trans,side_trans1,sdir1, side_trans2,sdir2 = direction_info(direction,pdir)

	# relevant other cells
	dirvec = (i_vec + bond_trans) % n_vec
	mdirvec = (i_vec - bond_trans) % n_vec
	msidevec1 = (i_vec - side_trans1) % n_vec
	msidevec2 = (i_vec - side_trans2) % n_vec

	
	
	Idir = {"x":1,"y":2,"z":3}

	I0 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[direction]]
	Pb0 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][0]
	I2 = atom_positions[mdirvec[0]][mdirvec[1]][mdirvec[2]][Idir[direction]]
	
	if negdir:
		dirvec,mdirvec = mdirvec,dirvec
		I0,I2 = I2,I0
	
	Pb1 = atom_positions[dirvec[0]][dirvec[1]][dirvec[2]][0]
	I3 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[sdir1]]
	I4 = atom_positions[msidevec2[0]][msidevec2[1]][msidevec2[2]][Idir[sdir2]]
	I5 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[sdir2]]
	I6 = atom_positions[msidevec1[0]][msidevec1[1]][msidevec1[2]][Idir[sdir1]]

	
	if fitting:
		vals = [H_index1, H_index2, param_real,param_im]
	else:
		vals = [i_vec[0],i_vec[1],i_vec[2],Idir[direction]]
	if features["bond_len"]:
		r = periodic_fix(Pb0 - I0, box,tilt)
		if split_len == 1: 
			vals.append(np.linalg.norm(r))
			#vals.append(r_distance(Pb0,I0,box,tilt))
		elif split_len == 2: 
			r_par = abs(r[Idir[direction]-1])
			r_perp = (r[Idir[sdir1]-1]**2 + r[Idir[sdir2]-1]**2)**(1/2)
			vals.append(r_par)
			vals.append(r_perp)
		elif split_len == 3: 
			r_par = abs(r[Idir[direction]-1])
			r_side1, r_side2 = abs(r[Idir[sdir1]-1]), abs(r[Idir[sdir2]-1])
			vals.append(r_par)
			vals.append(r_side1)
			vals.append(r_side2)
	if features["neighbor_len"]: 
		vals.append(r_distance(I0,Pb1,box,tilt)) 
		vals.append(r_distance(Pb0,I2,box,tilt))
		# order of below sensitive to pi stuff
		vals.append(r_distance(Pb0,I3,box,tilt))
		vals.append(r_distance(Pb0,I4,box,tilt))
		vals.append(r_distance(Pb0,I6,box,tilt))
		vals.append(r_distance(Pb0,I5,box,tilt))
	if features["bond_angle"]:
		vals.append(r_angle(Pb0,I0,Pb1,box,tilt))
		vals.append(r_angle(I2,Pb0,I0,box,tilt))
	if features["neighbor_angle"]:
		# of of below sensitive to pi stuff
		vals.append(r_angle(I3,Pb0,I0,box,tilt))
		vals.append(r_angle(I4,Pb0,I0,box,tilt))
		vals.append(r_angle(I6,Pb0,I0,box,tilt))
		vals.append(r_angle(I5,Pb0,I0,box,tilt))
		
	return vals

MA_angle_default = True
MA_dist_default = True

def Pb_onsite_geometry(i_vec,n_vec,direction, atom_positions,box, tilt,fitting = False,soc = False, features = {"bond_len" : True, "bond_angle" : True, "MA_angles" : MA_angle_default,"MA_distances" : MA_dist_default},coul = False, vdw = False,atom_potentials = None,atom_hessians = None,hessian = False): 
	# here, direction refers to the p orbital. It's not ideal
	bond_trans,side_trans1,sdir1, side_trans2,sdir2 = direction_info(direction)

	# relevant other cells
	dirvec = (i_vec + bond_trans) % n_vec
	mdirvec = (i_vec - bond_trans) % n_vec
	msidevec1 = (i_vec - side_trans1) % n_vec
	msidevec2 = (i_vec - side_trans2) % n_vec
	
	Idir = {"s": 1,"x":1,"y":2,"z":3} # s uses the same order as x

	I0 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[direction]]
	Pb0 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][0]
	#Pb1 = atom_positions[dirvec[0]][dirvec[1]][dirvec[2]][0]
	I2 = atom_positions[mdirvec[0]][mdirvec[1]][mdirvec[2]][Idir[direction]]
	I3 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[sdir1]]
	I4 = atom_positions[msidevec2[0]][msidevec2[1]][msidevec2[2]][Idir[sdir2]]
	I5 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[sdir2]]
	I6 = atom_positions[msidevec1[0]][msidevec1[1]][msidevec1[2]][Idir[sdir1]]
	if features["MA_angles"] or features["MA_distances"]: 
		MA_cells = [i_vec, mdirvec, msidevec1, msidevec2] + [(i_vec + trans)% n_vec for trans in [-(mdirvec + msidevec1), -(mdirvec + msidevec2), -(msidevec1 + msidevec2), -(mdirvec + msidevec1+msidevec2)]]
		MA_pos = [(atom_positions[c[0]][c[1]][c[2]][4], atom_positions[c[0]][c[1]][c[2]][5]) for c in MA_cells]
	if fitting:
		vals = [H_index1, H_index2, param_real,param_im]
	else:
		vals = [i_vec[0],i_vec[1],i_vec[2],Idir[direction]]
	if soc: 
		spin1, spin2, pdir1, pdir2 = soc
		vals.append(spin1)
		vals.append(spin2)
		vals.append(pdir1)
		vals.append(pdir2)
	if features["bond_len"]:
		# order these two based on p direction
		vals.append(r_distance(Pb0,I0,box,tilt)) 
		vals.append(r_distance(Pb0,I2,box,tilt))

		vals.append(r_distance(Pb0,I3,box,tilt))
		vals.append(r_distance(Pb0,I4,box,tilt))
		vals.append(r_distance(Pb0,I6,box,tilt))
		vals.append(r_distance(Pb0,I5,box,tilt))
	if features["bond_angle"]:
		#vals.append(r_angle(Pb0,I0,Pb1,box,tilt))
		vals.append(r_angle(I0,Pb0,I2,box,tilt)) # both parallel

		vals.append(r_angle(I0,Pb0,I3,box,tilt)) # one parallel
		vals.append(r_angle(I0,Pb0,I4,box,tilt))
		vals.append(r_angle(I0,Pb0,I5,box,tilt))
		vals.append(r_angle(I0,Pb0,I6,box,tilt))
		vals.append(r_angle(I2,Pb0,I3,box,tilt))
		vals.append(r_angle(I2,Pb0,I4,box,tilt))
		vals.append(r_angle(I2,Pb0,I5,box,tilt))
		vals.append(r_angle(I2,Pb0,I6,box,tilt))

		vals.append(r_angle(I3,Pb0,I4,box,tilt))
		vals.append(r_angle(I3,Pb0,I5,box,tilt))
		vals.append(r_angle(I3,Pb0,I6,box,tilt))
		vals.append(r_angle(I4,Pb0,I5,box,tilt))
		vals.append(r_angle(I4,Pb0,I6,box,tilt))
		vals.append(r_angle(I5,Pb0,I6,box,tilt))
	if coul:
		vals.append(atom_potentials[i_vec[0]][i_vec[1]][i_vec[2]][0][0])
	if vdw:
		vals.append(atom_potentials[i_vec[0]][i_vec[1]][i_vec[2]][0][1])
	if features["MA_distances"]: 
		for ma in MA_pos: 
			vals.append(r_distance(periodic_fix((ma[0] + ma[1])/2, box,tilt), Pb0, box,tilt))
	if features["MA_angles"]: 
		for ma in MA_pos: 
			vals.append(r_projection(ma[1],ma[0], dirvec, box,tilt))
	if hessian:
		#for i in range(3): only one necessary
		#	index = (i + Idir[direction] - 1)%3
		#	vals.append(atom_hessians[i_vec[0]][i_vec[1]][i_vec[2]][0][index][index])
		index = Idir[direction] - 1 # 0 for x, 1 for y, 2 for z
		vals.append(atom_hessians[i_vec[0]][i_vec[1]][i_vec[2]][0][index][index])
		#vals.append(atom_hessians[i_vec[0]][i_vec[1]][i_vec[2]][0][(index+1)%3][(index+1)%3])
	return vals

def I_onsite_geometry(i_vec,n_vec,direction, atom_positions,box, tilt,pdir,fitting = False,features = {"bond_len" : True, "bond_angle" : True,"MA_angles" : MA_angle_default,"MA_distances" : MA_dist_default},coul = False, vdw = False,soc = False,atom_potentials = None,atom_hessians = None, hessian = False): 
	# direction is which direction it is relative to the Pb atom in its unit cell. pdir is which p orbital it is. Both take 'x','y',or 'z'
	bond_trans,side_trans1,sdir1, side_trans2,sdir2 = direction_info(direction)

	# relevant other cells
	dirvec = (i_vec + bond_trans) % n_vec
	mdirvec = (i_vec - bond_trans) % n_vec
	msidevec1 = (i_vec - side_trans1) % n_vec
	msidevec2 = (i_vec - side_trans2) % n_vec
	
	Idir = {"s":0, "x":1,"y":2,"z":3}

	I0 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][Idir[pdir]]
	Pb0 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][0]
	Pb1 = atom_positions[dirvec[0]][dirvec[1]][dirvec[2]][0]
	if features["MA_angles"] or features["MA_distances"]:
		msidevecboth = (i_vec - side_trans1 - side_trans2) % n_vec
		C1 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][4]
		N1 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][5]

		C2 = atom_positions[msidevec1[0]][msidevec1[1]][msidevec1[2]][4]
		N2 = atom_positions[msidevec1[0]][msidevec1[1]][msidevec1[2]][5]

		C3 = atom_positions[msidevec2[0]][msidevec2[1]][msidevec2[2]][4]
		N3 = atom_positions[msidevec2[0]][msidevec2[1]][msidevec2[2]][5]

		C4 = atom_positions[msidevecboth[0]][msidevecboth[1]][msidevecboth[2]][4]
		N4 = atom_positions[msidevecboth[0]][msidevecboth[1]][msidevecboth[2]][5]
	if fitting:
		vals = [H_index1, H_index2, param_real,param_im]
	else:
		vals = [i_vec[0],i_vec[1],i_vec[2],Idir[pdir]]
	if soc: 
		spin1, spin2, pdir1, pdir2 = soc
		vals.append(spin1)
		vals.append(spin2)
		vals.append(pdir1)
		vals.append(pdir2)
	vals.append(int(pdir == direction))
	if features["bond_len"]:
		# order these two based on p direction
		vals.append(r_distance(Pb0,I0,box,tilt)) 
		vals.append(r_distance(Pb1,I0,box,tilt))
		
	if features["bond_angle"]:
		#vals.append(r_angle(Pb0,I0,Pb1,box,tilt))
		vals.append(r_angle(Pb0,I0,Pb1,box,tilt)) # both parallel
	if coul:
		vals.append(atom_potentials[i_vec[0]][i_vec[1]][i_vec[2]][Idir[direction]][0])
	if vdw:
		vals.append(atom_potentials[i_vec[0]][i_vec[1]][i_vec[2]][Idir[direction]][1])
	if features["MA_distances"]: 
		vals.append(r_distance(MA_midpoint(C1,N1, box,tilt), I0, box,tilt))
		vals.append(r_distance(MA_midpoint(C2,N2, box,tilt), I0, box,tilt))
		vals.append(r_distance(MA_midpoint(C3,N3, box,tilt), I0, box,tilt))
		vals.append(r_distance(MA_midpoint(C4,N4, box,tilt), I0, box,tilt))
	if features["MA_angles"]: 
		v = {"x":np.array([1.0,0,0]), "y" : np.array([0,1.0,0]), "z" : np.array([0,0,1.0])}[pdir]
		vals.append(r_projection(N1,C1, v, box,tilt))
		vals.append(r_projection(N2,C2, v, box,tilt))
		vals.append(r_projection(N3,C3, v, box,tilt))
		vals.append(r_projection(N4,C4, v, box,tilt))
	if hessian: 
		#for i in range(3):
		#	index = (i + Idir[pdir] -1)%3
		#	vals.append(atom_hessians[i_vec[0]][i_vec[1]][i_vec[2]][Idir[direction]][index][index])
		index = Idir[pdir] - 1
		vals.append(atom_hessians[i_vec[0]][i_vec[1]][i_vec[2]][Idir[direction]][index][index])

	return vals

# not sure if these are really necessary 
Pb_I_bonds = {}
"""
bonds_vars = []
counter = 0

for Idir in ["x","y","z"]:
	bonds_vars.append(("sig",Idir))
	Pb_I_bonds[("sig",Idir)] =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, bond_len = True, neighbor_len = True, bond_angle = True, neighbor_angle = True : Pb_I_bond_geometry(i_vec,n_vec,Idir, atom_positions,box,tilt, fitting = fitting, bond_len = bond_len, neighbor_len = neighbor_len , bond_angle = bond_angle, neighbor_angle = neighbor_angle)
	for pdir in ["x","y","z"]: 
		if pdir != Idir: 
			Pb_I_bonds[("pi",Idir,pdir)] = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, bond_len = True, neighbor_len = True, bond_angle = True, neighbor_angle = True : Pb_I_bond_geometry(i_vec,n_vec,Idir, atom_positions,box,tilt, pdir = pdir,fitting = fitting, bond_len = bond_len, neighbor_len = neighbor_len , bond_angle = bond_angle , neighbor_angle = neighbor_angle)

Pb_onsites = {}
for pdir in ["s","x","y","z"]:
	Pb_onsites[(pdir)] = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, bond_len = True, bond_angle = True, coul = False, vdw = False : Pb_onsite_geometry(i_vec,n_vec,pdir, atom_positions,box,tilt,fitting = fitting, bond_len = bond_len , bond_angle = bond_angle, coul = coul, vdw = vdw)

I_onsites = {}
for Idir in ["x","y","z"]: 
	for pdir in ["x","y","z"]:
		I_onsites[(Idir,pdir)] = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, bond_len = True, bond_angle = True, coul = False, vdw = False : I_onsite_geometry(i_vec,n_vec,Idir, atom_positions,box,tilt,pdir,fitting = fitting, bond_len = bond_len , bond_angle = bond_angle, coul = coul, vdw = vdw)
"""
bond_features = {"bond_len" : True, "neighbor_len" : True, "bond_angle" : True, "neighbor_angle" :True }
bond_min_features = {"bond_len" : True, "neighbor_len" : False, "bond_angle" : False, "neighbor_angle" :False }
bond_features = bond_min_features
#sigma bonds, doesn't matter if they are sp or pp
sigsplitn = 2
Pb_I_bond_x_sig = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features: Pb_I_bond_geometry(i_vec,n_vec,"x", atom_positions,box,tilt, fitting = fitting, features = features,split_len = sigsplitn)
Pb_I_bond_y_sig = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"y", atom_positions,box,tilt, fitting = fitting, features = features,split_len = sigsplitn)
Pb_I_bond_z_sig = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"z", atom_positions,box,tilt, fitting = fitting,  features = features,split_len = sigsplitn)
Pb_I_bond_mx_sig = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"x", atom_positions,box,tilt, negdir = True,fitting = fitting, features = features,split_len = sigsplitn)
Pb_I_bond_my_sig = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"y", atom_positions,box,tilt, negdir = True,fitting = fitting, features = features, split_len = sigsplitn)
Pb_I_bond_mz_sig = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"z", atom_positions,box,tilt, negdir = True,fitting = fitting, features = features, split_len = sigsplitn)
# NEED NEGATIVE DIRECTION

# pi bonds
pisplitn = 3
Pb_I_bond_x_pi_y = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False,  features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"x", atom_positions,box,tilt, pdir = "y",fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_x_pi_z = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"x", atom_positions,box,tilt, pdir = "z",fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_y_pi_x = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"y", atom_positions,box,tilt, pdir = "x",fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_y_pi_z = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"y", atom_positions,box,tilt, pdir = "z",fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_z_pi_x = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"z", atom_positions,box,tilt, pdir = "x",fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_z_pi_y = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"z", atom_positions,box,tilt, pdir = "y",fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_mx_pi_y = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"x", atom_positions,box,tilt, pdir = "y",negdir = True, fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_mx_pi_z = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"x", atom_positions,box,tilt, pdir = "z",negdir = True,fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_my_pi_x = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"y", atom_positions,box,tilt, pdir = "x",negdir = True,fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_my_pi_z = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"y", atom_positions,box,tilt, pdir = "z",negdir = True,fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_mz_pi_x = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"z", atom_positions,box,tilt, pdir = "x",negdir = True,fitting = fitting, features = features,split_len = pisplitn)
Pb_I_bond_mz_pi_y = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = bond_features : Pb_I_bond_geometry(i_vec,n_vec,"z", atom_positions,box,tilt, pdir = "y",negdir = True,fitting = fitting, features = features,split_len = pisplitn)


onsite_features = {"bond_len" : True, "bond_angle" : True,"MA_angles" : MA_angle_default,"MA_distances" : MA_dist_default}
onsite_no_features = {"bond_len" : False, "bond_angle" : False,"MA_angles" : False,"MA_distances" : False}
# Pb onsites
Pb_onsite_s =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None: Pb_onsite_geometry(i_vec,n_vec,"s", atom_positions,box,tilt,fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
Pb_onsite_px =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None: Pb_onsite_geometry(i_vec,n_vec,"x", atom_positions,box,tilt,fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
Pb_onsite_py =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None: Pb_onsite_geometry(i_vec,n_vec,"y", atom_positions,box,tilt,fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
Pb_onsite_pz =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None: Pb_onsite_geometry(i_vec,n_vec,"z", atom_positions,box,tilt,fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)

# I onsites
I_x_onsite_px =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"x", atom_positions,box,tilt,"x",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_x_onsite_py =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"x", atom_positions,box,tilt,"y",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_x_onsite_pz =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"x", atom_positions,box,tilt,"z",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_y_onsite_px =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"y", atom_positions,box,tilt,"x",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_y_onsite_py =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"y", atom_positions,box,tilt,"y",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_y_onsite_pz =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"y", atom_positions,box,tilt,"z",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_z_onsite_px =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"z", atom_positions,box,tilt,"x",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_z_onsite_py =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"z", atom_positions,box,tilt,"y",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)
I_z_onsite_pz =  lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, features = onsite_features, coul = False, vdw = False, atom_potentials = None, hessian = False, atom_hessians = None : I_onsite_geometry(i_vec,n_vec,"z", atom_positions,box,tilt,"z",fitting = fitting, features = features, coul = coul, vdw = vdw,atom_potentials = atom_potentials, hessian = hessian, atom_hessians = atom_hessians)

"""
# spin orbit
spin_orbitals = [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3)] #first is 1-up, 2-down, second is 1-x,2-y,3-z
Pb_soc = {}
I_soc = {}
p_orb_key = {1:"x", 2:"y", 3:"z"}
for i in range(len(spin_orbitals)): 
	for j in range(i, len(spin_orbitals)): 
		if i != j: 
			orb1, orb2 = spin_orbitals[i],spin_orbitals[j]
			Pb_soc[(orb1,orb2)] = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, bond_len = True, bond_angle = True, coul = False, vdw = False : Pb_onsite_geometry(i_vec,n_vec,"s", atom_positions,box,tilt,fitting = fitting, bond_len = bond_len , bond_angle = bond_angle, coul = coul, vdw = vdw, soc = [orb1[0],orb2[0],orb1[1],orb2[1]])
			for direction in ["x","y","z"]: 
				I_soc[(direction, orb1,orb2)] = lambda i_vec, n_vec, atom_positions,box,tilt, fitting = False, bond_len = True, bond_angle = True, coul = False, vdw = False : Pb_onsite_geometry(i_vec,n_vec,"s", atom_positions,box,tilt,fitting = fitting, bond_len = bond_len , bond_angle = bond_angle, coul = coul, vdw = vdw, soc = [orb1[0],orb2[0],orb1[1],orb2[1]])


spsig = [Pb_I_bonds[("sig","x")], Pb_I_bonds[("sig","y")], Pb_I_bonds[("sig","z")]]
ppsig = spsig
pppi = [Pb_I_bonds[("pi","x","y")],Pb_I_bonds[("pi","x","z")],Pb_I_bonds[("pi","y","x")],Pb_I_bonds[("pi","y","z")],Pb_I_bonds[("pi","z","x")],Pb_I_bonds[("pi","z","y")]]
pbs_onsites = [Pb_onsites[("s")]]
pbp_onsites = [Pb_onsites[("x")],Pb_onsites[("y")], Pb_onsites[("z")]]
ip_onsites = [I_onsites[("x","x")],I_onsites[("x","y")],I_onsites[("x","z")],I_onsites[("y","x")],I_onsites[("y","y")],I_onsites[("y","z")],I_onsites[("z","x")],I_onsites[("z","y")],I_onsites[("z","z")]]
"""
spsig = [Pb_I_bond_x_sig, Pb_I_bond_y_sig, Pb_I_bond_z_sig,Pb_I_bond_mx_sig, Pb_I_bond_my_sig, Pb_I_bond_mz_sig]
ppsig = spsig
pppi = [Pb_I_bond_x_pi_y, Pb_I_bond_x_pi_z, Pb_I_bond_y_pi_x, Pb_I_bond_y_pi_z, Pb_I_bond_z_pi_x, Pb_I_bond_z_pi_y,Pb_I_bond_mx_pi_y, Pb_I_bond_mx_pi_z, Pb_I_bond_my_pi_x, Pb_I_bond_my_pi_z, Pb_I_bond_mz_pi_x, Pb_I_bond_mz_pi_y]
pbs_onsites = [Pb_onsite_s]
pbp_onsites = [Pb_onsite_px,Pb_onsite_py, Pb_onsite_pz]
ip_onsites = [I_x_onsite_px, I_x_onsite_py,I_x_onsite_pz, I_y_onsite_px, I_y_onsite_py, I_y_onsite_pz, I_z_onsite_px, I_z_onsite_py, I_z_onsite_pz]


def hess(dipoles, monopoles, diagonly = True): 
	H = np.zeros([3,3])
	def helement(i,j): 
		tot = 0.0
		for dip in dipoles:
			tot += diphess(dip[0],dip[1], i,j)
		for mon in monopoles:
			tot += monhess(mon[0],mon[1], i,j)
		return tot
	for i in range(3):
		H[i][i] = helement(i,i)
	if not diagonly:
		for i in range(3):
			for j in range(i+1,3):
				H[i][j] = helement(i,j)
				H[j][i] = H[i][j]
	return H


def monhess(q, rvec, i, j):
	# hessian element d^2V / dxidxj for a monopole charge a vector r away from the atom
	# i,j = 0,1,2 for x,y,z
	k = 14.39964812 # units of eV*Angsrom / electron charge^2
	if i == j:
		x,y,z = rvec[i],rvec[(i+1)%3],rvec[(i+2)%3]
		r = np.linalg.norm(rvec)
		return k*q*(3*x**2 - r**2) / r**5
	else:
		x,y,z = rvec[i],rvec[j],rvec[3 - i - j]
		r = np.linalg.norm(rvec)
		return k*3*q*x*y/r**5

def diphess(p, rvec, i, j): 
	# hessian element d^2V / dxidxj for a dipole with moment p a vector r away from the atom
	# i, j = 0,1,2 for x,y,z
	k = 14.39964812 # units of eV*Angsrom / electron charge^2
	if i == j:
		x,y,z = rvec[i],rvec[(i+1)%3],rvec[(i+2)%3]
		px,py,pz = p[i],p[(i+1)%3],p[(i+2)%3]
		r = np.linalg.norm(rvec)
		return k*3*(px*x*(5*x**2 - 3*r**2) + (5*x**2 - r**2)*(py*y + pz*z)) / r**7
	else:
		x,y,z = rvec[i],rvec[j],rvec[3 - i - j]
		px,py,pz = p[i], p[j], p[3 - i - j]
		r = np.linalg.norm(rvec)
		return k*3*(px*y*(5*x**2 - r**2) + py*x*(5*y*2 - r**2) + 5*pz*x*y*z) / r**7

pbnnsetting, innsetting,madipsetting,inchsetting = True,True, True, False
def pbhess(i_vec, n_vec,atom_positions,box,tilt, charges = "default",diagonly = True,pb_next_nearest = pbnnsetting, i_next_nearest = innsetting,madip = madipsetting,inch = inchsetting):
	ix,iy,iz = i_vec
	nx,ny,nz = n_vec
	pb = atom_positions[ix][iy][iz][0]
	ixp = atom_positions[ix][iy][iz][1]
	ixm = atom_positions[(ix-1)%nx][iy][iz][1]
	iyp = atom_positions[ix][iy][iz][2]
	iym = atom_positions[ix][(iy-1)%ny][iz][2]
	izp = atom_positions[ix][iy][iz][3]
	izm = atom_positions[ix][iy][(iz-1)%nz][3]
	if pb_next_nearest: 
		nnpb = [atom_positions[(ix+tx)%nx][(iy+ty)%ny][(iz+tz)%nz][0] for tx,ty,tz in [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]]
	if i_next_nearest:
		nnI = [] # 24 of these
		for pbdr in range(3): 
			for idr in range(3): 
				if pbdr != idr:
					for pbtrans in [-1,1]:
						for itrans in [0,-1]:
							pbvec, ivec = np.zeros(3,np.int32), np.zeros(3,np.int32)
							pbvec[pbdr] = pbtrans
							ivec[idr] = itrans
							cell = (i_vec + pbvec + ivec) % n_vec
							nnI.append(atom_positions[cell[0]][cell[1]][cell[2]][idr+1])						
	MA_pos = []
	for tx in range(-1,1):
		for ty in range(-1,1):
			for tz in range(-1,1):
				MA_pos.append(atom_positions[(ix+tx)%nx][(iy+ty)%ny][(iz+tz)%nz][4:])
				#MA_pos.append((atom_positions[(ix+tx)%nx][(iy+ty)%ny][(iz+tz)%nz][4],atom_positions[(ix+tx)%nx][(iy+ty)%ny][(iz+tz)%nz][5]))
	if charges == "default": # add in charges in the same time as positions if these become variable. 
		qi = -1.13
		monopoles = [(qi,periodic_fix(pb - i, box,tilt)) for i in [ixp,ixm,iyp,iym,izp,izm]] + ([(MA_charge,periodic_fix(pb - MA_midpoint(c,n,box,tilt), box,tilt)) for c,n in MA_pos] if madip else [])
		if pb_next_nearest: 
			qpb = 2.03 
			monopoles += [(qpb, periodic_fix(pb - i, box, tilt)) for i in nnpb]
		if i_next_nearest: 
			monopoles += [(qi, periodic_fix(pb - i, box, tilt)) for i in nnI]
	# method for varying charges not currently implemented
	def ma_dipole(c, n,pb):
		d = periodic_fix(c-n, box,tilt)
		phat = d / np.linalg.norm(d)
		p = MA_DIPOLE * phat
		r = periodic_fix(pb - (d/2 + c), box, tilt)
		return (p,r)
	if madip:
		dipoles = [ma_dipole(ma[0],ma[1],pb) for ma in MA_pos]
	else: 
		charges = [0.771, -1.1, 0.023,0.023,0.023,0.54,0.54,0.54] if len(MA_pos[0]) > 2 else [0.84, 0.52]
		for ma in MA_pos:
			for i in range(len(ma)):
				monopoles.append((charges[i], periodic_fix(pb - ma[i], box,tilt)))
		dipoles = []
	#print(len(monopoles))
	return hess(dipoles, monopoles, diagonly)

def ihess(i_vec, n_vec,direction,atom_positions,box,tilt, charges = "default",diagonly = True, inc_I_nn = True):
	ix,iy,iz = i_vec
	nx,ny,nz = n_vec
	bond_trans,side_trans1,sdir1, side_trans2,sdir2 = direction_info(direction)

	# relevant other cells
	dirvec = (i_vec + bond_trans) % n_vec
	mdirvec = (i_vec - bond_trans) % n_vec
	msidevec1 = (i_vec - side_trans1) % n_vec
	msidevec2 = (i_vec - side_trans2) % n_vec
	Idir = {"x":1,"y":2,"z":3}

	I = atom_positions[ix][iy][iz][Idir[direction]]
	Pb0 = atom_positions[ix][iy][iz][0]
	Pb1 = atom_positions[dirvec[0]][dirvec[1]][dirvec[2]][0]
	nnI = [atom_positions[(ix+t[0])%nx][(iy+t[1])%ny][(iz+t[2])%2][Idir[idr]] for t,idr in [(side_trans1,sdir1),(-side_trans1,sdir1),(side_trans2,sdir2),(-side_trans2,sdir2),(dirvec + side_trans1,sdir1),(dirvec-side_trans1,sdir1),(dirvec + side_trans2,sdir2),(dirvec-side_trans2,sdir2)]]
	# nearest neigbor iodines, separated by about sqrt(2)/2 unit cells

	msidevecboth = (i_vec - side_trans1 - side_trans2) % n_vec
	C1 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][4]
	N1 = atom_positions[i_vec[0]][i_vec[1]][i_vec[2]][5]

	C2 = atom_positions[msidevec1[0]][msidevec1[1]][msidevec1[2]][4]
	N2 = atom_positions[msidevec1[0]][msidevec1[1]][msidevec1[2]][5]

	C3 = atom_positions[msidevec2[0]][msidevec2[1]][msidevec2[2]][4]
	N3 = atom_positions[msidevec2[0]][msidevec2[1]][msidevec2[2]][5]

	C4 = atom_positions[msidevecboth[0]][msidevecboth[1]][msidevecboth[2]][4]
	N4 = atom_positions[msidevecboth[0]][msidevecboth[1]][msidevecboth[2]][5]

	if charges == "default":
		qpb = 2.03
		qi = -1.13
		monopoles = [(qpb,periodic_fix(I - i, box,tilt)) for i in [Pb0,Pb1]] + [(MA_charge,periodic_fix(I - MA_midpoint(c,n,box,tilt),box,tilt)) for c,n in [(C1,N1),(C2,N2),(C3,N3),(C4,N4)]]
		if inc_I_nn:
			monopoles += [(qi,periodic_fix(I - i, box,tilt)) for i in nnI]
	# method for varying charges not currently implemented
	def ma_dipole(c, n,I):
		d = periodic_fix(c-n, box,tilt)
		phat = d / np.linalg.norm(d)
		p = MA_DIPOLE * phat
		r = periodic_fix(I - (d/2 + c), box, tilt)
		return (p,r)
	dipoles = [ma_dipole(C1,N1,I),ma_dipole(C2,N2,I),ma_dipole(C3,N3,I),ma_dipole(C4,N4,I)]
	return hess(dipoles, monopoles, diagonly)


def calc_MA_dipole():
	C, N, Hc1,Hc2,Hc3, Hn1,Hn2,Hn3 = [structure_tools.tilted_n2.handle_lookup([(0,0,0),"MA",i]).position for i in ["C","N","H_C1","H_C2","H_C3","H_N1","H_N2","H_N3"]]
	qC, qN, qHc, qHn = 0.771, -1.1, 0.023, 0.54
	center = (qC*C + qN*N + qHc*(Hc1+Hc2+Hc3) + qHn*(Hn1+Hn2+Hn3))/(qC + qN + 3*qHc + 3*qHn)
	return np.linalg.norm(center - (C+N)/2)

def MA_midpoint(c,n,box,tilt): 
	return c + periodic_fix(n-c,box,tilt)/2

#MA_DIPOLE = calc_MA_dipole()
MA_charge = 0.771 + -1.1 + 3*(0.023 + 0.54)

MA_DIPOLE = -0.479 #2.3 debye, converted to e*A. From Frost et al, Atomistic Origins of High-Performance in Hybrid Halide Perovskite Solar Cells
# N is positive, C is negative so that gives it a negative sign. 


def mond4(q,r):
	# 4th derivative element d^4V/dxi^2dxj^2 for a distance r from a monopole. 
	# i,j = 0,1,2 for x,y,z
	k = 14.39964812 # units of eV*Angsrom / electron charge^2
	if i == j:
		x,y,z = rvec[i],rvec[(i+1)%3],rvec[(i+2)%3]
		r = np.linalg.norm(rvec)
		return (24*x**4 - 72*x**2*(r**2 - x**2) + 9*(r**2 - x**2)**2)/(r**9)
	else:
		x,y,z = rvec[i],rvec[j],rvec[3 - i - j]
		r = np.linalg.norm(rvec)
		return 3*(-4*x**2 + -4*y**2 + 27*x**2*y**2 - 3*(x**2 + y**2)*z**2 + z**4)/(r**9)
