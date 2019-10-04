import numpy as np


def str_to_tup(cellstr): 
	return tuple([int(i) for i in cellstr.split(":")])

def tup_to_str(celltup):
	return ":".join(str(i) for i in y)

def stoich_to_unit_index(nx,ny,nz,N = 8): 
	# takes unit cell indices nx,ny,nz for a Pb atom and orthorhombic box size N and gives the Pb orthorhombic handle
	cell = (((nx-nz)//2)%N, (ny//2)%N,((nx + nz)//2)%N)
	Pb_num = 1 + 2*abs(nx - nz)%4 + ny%2 
	#return [str(cell[0]) + ":" + str(cell[1]) + ":" + str(cell[2]), "Pb" + str(Pb_num)]
	return [cell, "Pb" + str(Pb_num)]

def I_x_ortho(Pb_index, N = 8): 
	#Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	Pb_cell, Pb_num = Pb_index[0], int(Pb_index[1][2])
	if Pb_num == 1: 
		I_cell,I_num = Pb_cell, 2
	elif Pb_num == 2: 
		I_cell,I_num = Pb_cell, 3
	elif Pb_num == 3: 
		#I_cell,I_num = [Pb_cell[0],(Pb_cell[1] - 1)%N, Pb_cell[2]], 10
		I_cell,I_num = (Pb_cell[0],(Pb_cell[1] - 1)%N, Pb_cell[2]), 10
	elif Pb_num == 4: 
		I_cell,I_num = Pb_cell, 11
	#return [str(I_cell[0]) + ":" + str(I_cell[1]) + ":" + str(I_cell[2]), "I" + str(I_num)]
	return [I_cell, "I" + str(I_num) ]

def I_y_ortho(Pb_index, N = 8):
	#Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	Pb_cell, Pb_num = Pb_index[0], int(Pb_index[1][2])
	if Pb_num == 1: 
		#I_cell,I_num = [Pb_cell[0],Pb_cell[1],(Pb_cell[2]-1)% N], 1
		I_cell,I_num = (Pb_cell[0],Pb_cell[1],(Pb_cell[2]-1)% N), 1
	elif Pb_num == 2: 
		#I_cell,I_num = [(Pb_cell[0]-1)%N, Pb_cell[1], Pb_cell[2]], 12
		I_cell, I_num = ((Pb_cell[0]-1)%N, Pb_cell[1], Pb_cell[2]), 12
	elif Pb_num == 3: 
		I_cell,I_num = Pb_cell, 7
	elif Pb_num == 4: 
		I_cell,I_num = Pb_cell, 6
	#return [str(I_cell[0]) + ":" + str(I_cell[1]) + ":" + str(I_cell[2]), "I" + str(I_num)]
	return [I_cell, "I" + str(I_num)]

def I_z_ortho(Pb_index, N = 8):
	#Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	Pb_cell, Pb_num = Pb_index[0], int(Pb_index[1][2])
	if Pb_num == 1: 
		#I_cell,I_num = [(Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]], 8
		I_cell,I_num = ((Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]), 8
	elif Pb_num == 2: 
		#I_cell,I_num = [(Pb_cell[0]-1)%N, Pb_cell[1], Pb_cell[2]], 9
		I_cell,I_num = ((Pb_cell[0]-1)%N, Pb_cell[1], Pb_cell[2]), 9
	elif Pb_num == 3: 
		#I_cell,I_num = [Pb_cell[0],(Pb_cell[1]-1)%N,Pb_cell[2]], 4
		I_cell,I_num = (Pb_cell[0],(Pb_cell[1]-1)%N,Pb_cell[2]), 4
	elif Pb_num == 4: 
		I_cell,I_num = Pb_cell, 5
	#return [str(I_cell[0]) + ":" + str(I_cell[1]) + ":" + str(I_cell[2]), "I" + str(I_num)]
	return [I_cell, "I" + str(I_num)]

def MA_ortho(Pb_index,N = 8): 
	#Pb_cell, Pb_num = [int(i) for i in Pb_index[0].split(":")], int(Pb_index[1][2])
	Pb_cell, Pb_num = Pb_index[0], int(Pb_index[1][2])
	if Pb_num == 1:
		#MA_cell,MA_num = [(Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]],1
		MA_cell,MA_num = ((Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]),1
	elif Pb_num == 2: 
		#MA_cell, MA_num = [(Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]], 4
		MA_cell, MA_num = ((Pb_cell[0]-1)%N,Pb_cell[1],Pb_cell[2]), 4
	elif Pb_num == 3: 
		MA_cell,MA_num = Pb_cell, 3
	elif Pb_num == 4: 
		MA_cell, MA_num = Pb_cell, 2
	#return [str(MA_cell[0]) + ":" + str(MA_cell[1]) + ":" + str(MA_cell[2]), "MA" + str(MA_num)]
	return [MA_cell, "MA" + str(MA_num)]

N_ortho,C_ortho = lambda Pb_index,N : MA_ortho(Pb_index, N) + ["N"],lambda Pb_index,N : MA_ortho(Pb_index, N) + ["C"]
Hc1_ortho, Hc2_ortho,Hc3_ortho = lambda Pb_index,N : MA_ortho(Pb_index, N) + ["H_C1"], lambda Pb_index,N :  MA_ortho(Pb_index, N) + ["H_C2"],lambda Pb_index,N : MA_ortho(Pb_index, N) + ["H_C3"]
Hn1_ortho, Hn2_ortho,Hn3_ortho = lambda Pb_index,N : MA_ortho(Pb_index, N) + ["H_N1"], lambda Pb_index,N :  MA_ortho(Pb_index, N) + ["H_N2"],lambda Pb_index,N : MA_ortho(Pb_index, N) + ["H_N3"]
			

"""def periodic_fix(r, box_vect):#, verbose =False): 
	t = np.array([0.0,0.0,0.0]) # translation, so that function is not mutative
	for i in range(3): 
		if r[i] > box_vect[i]/2: 
			t[i] =  - box_vect[i]
		elif r[i] < -box_vect[i]/2: 
			t[i] = box_vect[i]
	#if verbose: 
	#	print(r,t)
	return r + t"""

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
	"""for i in range(3): 
		if r[i] > box_vect[i]/2: 
			t[i] =  - box_vect[i]
		elif r[i] < -box_vect[i]/2: 
			t[i] = box_vect[i]"""
	#if verbose: 
	#	print(r,t)
	return r + t

def stoich_to_unit_conversion(stoich_handle, N = 6, xshift = 0, yshift = 0, zshift = 0): 
	#print(stoich_handle)
	#nx,ny,nz = [int(i) for i in stoich_handle[0].split(":")]
	nx,ny,nz = stoich_handle[0]
	nx += xshift
	ny += yshift
	nz += zshift
	Pb_index = stoich_to_unit_index(nx,ny,nz, N)
	if stoich_handle[1] == "Pb": 
		return Pb_index
	elif stoich_handle[1] == "I1": 
		return I_x_ortho(Pb_index, N)
	elif stoich_handle[1] == "I2": 
		return I_y_ortho(Pb_index,N)
	elif stoich_handle[1] == "I3":
		return I_z_ortho(Pb_index,N)
	elif stoich_handle[1] == "MA":
		return MA_ortho(Pb_index,N) + stoich_handle[2:]

def ortho_to_stoich_transform(ortho_vect,cos): 
	sin = np.sqrt(1 - cos**2)
	return np.array([[cos, 0, sin],[0, 1, 0],[-sin, 0, cos]]) @ ortho_vect

