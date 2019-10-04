import numpy as np
import lammps_analysis_tools as analysis_tools
import lammps_structure_tools as structure_tools
from index_conversions import * 

spins = ["up","down"]
wan_types = {"Pb:p":structure_tools.Pb,"I:p":structure_tools.I,"Pb:s":structure_tools.Pb,"C:p":structure_tools.C,"N:p" : structure_tools.N}


def make_wan_map(structure, spins = spins, wan_types = wan_types,handle_map = True): 
	wan_map = {}
	counter = 1
	
	for typ in wan_types:
		for atom in structure.atoms:
			if wan_types[typ] is atom.atom_type:
				if typ[-1] == 'p':
					for i in ['z','x','y']: # why does wannier do it this way!?!
						for spin in spins:
							wan_map[counter] = [spin,typ,atom,i]
							counter += 1
				elif typ[-1] == 's':
					for spin in spins:
						wan_map[counter] = [spin,typ,atom,'s']
						counter += 1
	if handle_map:
		hmap = {}
		for wan in wan_map: 
			info = wan_map[wan]
			hmap[(info[0],info[1],tuple(info[2].hpath),info[3])] = wan
		return wan_map,hmap
	else: 
		return wan_map


def read_hr(hr_file): # need dump or structure ? 
	lines = open(hr_file,'r').read()
	lines = lines.split("\n")
	num_wan,num_r = int(lines[1]),int(lines[2])
	mults = sum([[int(i) for i in line.split("    ") if i != ''] for line in lines[3:4 + (num_r//15)]], [])

	nr = int(round(num_r**(1.0/3)))
	H = [[[np.zeros((num_wan + 1,num_wan + 1), dtype = np.complex64) for nz in range(nr)] for ny in range(nr)] for nx in range(nr)]

	for line in lines[4 + (num_r //15): -1]: #last line is ""
		nx,ny,nz,i,j,real,im = [float(i) for i in line.split(" ") if i != ""]
		nx,ny,nz,i,j = int(nx),int(ny),int(nz),int(i),int(j)
		H[nx][ny][nz][i][j] = real + 1j*im
	return H

def write_hops_dat(hr_file, structure,outfile,ortho_cell = False,n = None, n_vec = None): 
	H = read_hr(hr_file)
	wan_map,hmap = make_wan_map(structure)
	lines = []
	data = []
	N_atoms = len(structure.atoms)
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
	for nx in range(nx_max): 
			for ny in range(ny_max): 
				for nz in range(nz_max): 
					for spin in spins:
						Pb_index = tuple(index_conversion(nx,ny,nz,N))
						Ix_plus_index, Iy_plus_index, Iz_plus_index = tuple(I_x(Pb_index,N)), tuple(I_y(Pb_index, N)), tuple(I_z(Pb_index, N))
						# wrong, kept around for now.Ix_minus_index, Iy_minus_index, Iz_minus_index = I_x(index_conversion((nx-1)%n, ny, nz,N),N), I_y(index_conversion(nx, (ny-1)%n, nz,N),N), I_z(index_conversion(nx, ny, (nz-1)%n,N),N)
						Ix_minus_index, Iy_minus_index, Iz_minus_index = tuple(I_x(index_conversion((nx-1), ny, nz,N),N)), tuple(I_y(index_conversion(nx, (ny-1), nz,N),N)), tuple(I_z(index_conversion(nx, ny, (nz-1),N),N))
						# ^ probably fixed
						#coul_Pb = kcal_to_eV(2 * dump.atom_pot(Pb_index, step, id_mode = "handles") / Q_PB)
						#coul_Ix, coul_Iy, coul_Iz = kcal_to_eV(2 * dump.atom_pot(Ix_plus_index,step,id_mode = "handles")/ Q_I), kcal_to_eV(2 * dump.atom_pot(Iy_plus_index,step, id_mode = "handles" )/ Q_I), kcal_to_eV(2 * dump.atom_pot(Iz_plus_index, step, id_mode = "handles")/ Q_I)

						Pb_s = hmap[(spin,'Pb:s',Pb_index,'s')]
						Pb_px,Pb_py,Pb_pz = hmap[(spin,'Pb:p',Pb_index,'x')], hmap[(spin,'Pb:p',Pb_index,'y')], hmap[(spin,'Pb:p',Pb_index,'z')]
						Ix_plus_px,Ix_plus_py,Ix_plus_pz = hmap[(spin,"I:p",Ix_plus_index,'x')],hmap[(spin,"I:p",Ix_plus_index,'y')],hmap[(spin,"I:p",Ix_plus_index,'z')]
						Iy_plus_px,Iy_plus_py,Iy_plus_pz = hmap[(spin,"I:p",Iy_plus_index,'x')],hmap[(spin,"I:p",Iy_plus_index,'y')],hmap[(spin,"I:p",Iy_plus_index,'z')]
						Iz_plus_px,Iz_plus_py,Iz_plus_pz = hmap[(spin,"I:p",Iz_plus_index,'x')],hmap[(spin,"I:p",Iz_plus_index,'y')],hmap[(spin,"I:p",Iz_plus_index,'z')]
						Ix_minus_px,Ix_minus_py,Ix_minus_pz = hmap[(spin,"I:p",Ix_minus_index,'x')],hmap[(spin,"I:p",Ix_minus_index,'y')],hmap[(spin,"I:p",Ix_minus_index,'z')]
						Iy_minus_px,Iy_minus_py,Iy_minus_pz = hmap[(spin,"I:p",Iy_minus_index,'x')],hmap[(spin,"I:p",Iy_minus_index,'y')],hmap[(spin,"I:p",Iy_minus_index,'z')]
						Iz_minus_px,Iz_minus_py,Iz_minus_pz = hmap[(spin,"I:p",Iz_minus_index,'x')],hmap[(spin,"I:p",Iz_minus_index,'y')],hmap[(spin,"I:p",Iz_minus_index,'z')]
						

						#lines.append(str(E_PB_S - (coul_Pb - Pb_coul_ave)) + (" " + str(E_PB_P - (coul_Pb - Pb_coul_ave)))*3)
						def onsite(n):
							return np.real(H[0][0][0][n][n])

						lines.append(str(onsite(Pb_s)) + " " + str(onsite(Pb_px)) + " " + str(onsite(Pb_py)) + " " + str(onsite(Pb_pz)))
						lines.append(str(onsite(Ix_plus_px)) + " " + str(onsite(Ix_plus_px)) + " " + str(onsite(Ix_plus_py)) + " " + str(onsite(Ix_plus_pz)))
						lines.append(str(onsite(Iy_plus_px)) + " " + str(onsite(Iy_plus_px)) + " " + str(onsite(Iy_plus_py)) + " " + str(onsite(Iy_plus_pz)))
						lines.append(str(onsite(Iz_plus_px)) + " " + str(onsite(Iz_plus_px)) + " " + str(onsite(Iz_plus_py)) + " " + str(onsite(Iz_plus_pz)))
						
						def hopping(n1,n2):
							hops = []
							for nx in range(-1,2):
								for ny in range(-1,2):
									for nz in range(-1,2):
										hops.append(H[nx][ny][nz][n1][n2].real)
							return(max(hops))


						lines.append(str(hopping(Pb_s,Ix_plus_px)) + " " + str(hopping(Pb_px,Ix_plus_px)) + " " + str(hopping(Pb_py,Ix_plus_py)) + " " + str(hopping(Pb_pz,Ix_plus_pz)))
						lines.append(str(hopping(Pb_s,Iy_plus_py)) + " " + str(hopping(Pb_py,Iy_plus_py)) + " " + str(hopping(Pb_px,Iy_plus_px)) + " " + str(hopping(Pb_pz,Iy_plus_pz)))
						lines.append(str(hopping(Pb_s,Iz_plus_pz)) + " " + str(hopping(Pb_pz,Iz_plus_pz)) + " " + str(hopping(Pb_px,Iz_plus_px)) + " " + str(hopping(Pb_py,Iz_plus_py)))
						lines.append(str(hopping(Pb_s,Ix_minus_px)) + " " + str(hopping(Pb_px,Ix_minus_px)) + " " + str(hopping(Pb_py,Ix_minus_py)) + " " + str(hopping(Pb_pz,Ix_minus_pz)))
						lines.append(str(hopping(Pb_s,Iy_minus_py)) + " " + str(hopping(Pb_py,Iy_minus_py)) + " " + str(hopping(Pb_px,Iy_minus_px)) + " " + str(hopping(Pb_pz,Iy_minus_pz)))
						lines.append(str(hopping(Pb_s,Iz_minus_pz)) + " " + str(hopping(Pb_pz,Iz_minus_pz)) + " " + str(hopping(Pb_px,Iz_minus_px)) + " " + str(hopping(Pb_py,Iz_minus_py)))
	f = open(outfile, "w")
	f.write("\n".join(lines))# apparently .join is the most efficient way to write large strings

def hopping(H,n1,n2,returnn = True):
	hops = []
	for nx in range(-1,2):
		for ny in range(-1,2):
			for nz in range(-1,2):
				hops.append(H[nx][ny][nz][n1][n2])
	return (max(hops,key = lambda i : i.real), n1, n2) if returnn else max(hops, key = lambda i: i.real)

def onsite(H,n,returnn = True):
	return (np.real(H[0][0][0][n][n]), n,n) if returnn else np.real(H[0][0][0][n][n])

index_conversion = lambda nx,ny,nz,N = None: [(nx%N, ny%N, nz%N), "Pb"]
I_x,I_y,I_z = lambda Pb_index,N = None : [Pb_index[0], "I1"], lambda Pb_index, N = None : [Pb_index[0], "I2"], lambda Pb_index, N = None : [Pb_index[0], "I3"]

#sigma bonds, sp
Pb_I_bond_x_sigsp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'x')])
Pb_I_bond_y_sigsp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'y')])
Pb_I_bond_z_sigsp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'z')])
Pb_I_bond_mx_sigsp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx-1,ny,nz,2))),'x')])
Pb_I_bond_my_sigsp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny-1,nz,2))),'y')])
Pb_I_bond_mz_sigsp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz-1,2))),'z')])


#sigma bonds, pp
Pb_I_bond_x_sigpp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'x')])
Pb_I_bond_y_sigpp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'y')])
Pb_I_bond_z_sigpp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'z')])
Pb_I_bond_mx_sigpp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx-1,ny,nz,2))),'x')])
Pb_I_bond_my_sigpp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny-1,nz,2))),'y')])
Pb_I_bond_mz_sigpp = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz-1,2))),'z')])
# need bonds in negative direction

# pi bonds
Pb_I_bond_x_pi_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'y')])
Pb_I_bond_x_pi_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'z')])
Pb_I_bond_y_pi_x = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'x')])
Pb_I_bond_y_pi_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'z')])
Pb_I_bond_z_pi_x = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'x')])
Pb_I_bond_z_pi_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'y')])
Pb_I_bond_mx_pi_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx-1,ny,nz,2))),'y')])
Pb_I_bond_mx_pi_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'I:p',tuple(I_x(index_conversion(nx-1,ny,nz,2))),'z')])
Pb_I_bond_my_pi_x = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny-1,nz,2))),'x')])
Pb_I_bond_my_pi_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny-1,nz,2))),'z')])
Pb_I_bond_mz_pi_x = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz-1,2))),'x')])
Pb_I_bond_mz_pi_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz-1,2))),'y')])

# Next nearest neighbor hopping
# I I pp sigma
I_I_nn_sig_x = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'x')], hmap[(spin,'I:p',tuple(I_x(index_conversion((nx+1)%2,ny,nz,2))),'x')])
I_I_nn_sig_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,(ny+1)%2,nz,2))),'y')])
I_I_nn_sig_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,(nz+1)%2,2))),'z')])

# Pb Pb pp sigma 
Pb_Pb_nn_sig_x = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')], hmap[(spin,'Pb:p',tuple(index_conversion((nx+1)%2,ny,nz,2)),'x')])
Pb_Pb_nn_sig_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')], hmap[(spin,'Pb:p',tuple(index_conversion(nx,(ny+1)%2,nz,2)),'y')])
Pb_Pb_nn_sig_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')], hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,(nz+1)%2,2)),'z')])

# I I p p ortho
I_I_nn_x_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'x')])
I_I_nn_px_y = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion((nx+1)%2,ny,nz,2))),'x')])
I_I_nn_x_my = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,(ny-1)%2,nz,2))),'x')])
I_I_nn_px_my = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'y')], hmap[(spin,'I:p',tuple(I_y(index_conversion((nx+1)%2,(ny-1)%2,nz,2))),'x')])

I_I_nn_x_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'x')])
I_I_nn_px_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion((nx+1)%2,ny,nz,2))),'x')])
I_I_nn_x_mz = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,(nz-1)%2,2))),'x')])
I_I_nn_px_mz = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion((nx+1),ny,(nz-1)%2,2))),'x')])

I_I_nn_y_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'y')])
I_I_nn_py_z = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,(ny+1)%2,nz,2))),'y')])
I_I_nn_y_mz = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,(nz-1)%2,2))),'y')])
I_I_nn_py_mz = lambda H,hmap,nx,ny,nz,spin : hopping(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'z')], hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,(ny+1)%2,(nz-1)%2,2))),'y')])


# Pb onsites
Pb_onsite_s =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'Pb:s',tuple(index_conversion(nx,ny,nz,2)),'s')])
Pb_onsite_px =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'x')])
Pb_onsite_py =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'y')])
Pb_onsite_pz =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'Pb:p',tuple(index_conversion(nx,ny,nz,2)),'z')])

# I onsites
I_x_onsite_px =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'x')])
I_x_onsite_py =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'y')])
I_x_onsite_pz =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_x(index_conversion(nx,ny,nz,2))),'z')])
I_y_onsite_px =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'x')])
I_y_onsite_py =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'y')])
I_y_onsite_pz =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_y(index_conversion(nx,ny,nz,2))),'z')])
I_z_onsite_px =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'x')])
I_z_onsite_py =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'y')])
I_z_onsite_pz =  lambda H,hmap,nx,ny,nz,spin : onsite(H, hmap[(spin,'I:p',tuple(I_z(index_conversion(nx,ny,nz,2))),'z')])

spsig = [Pb_I_bond_x_sigsp, Pb_I_bond_y_sigsp, Pb_I_bond_z_sigsp,Pb_I_bond_mx_sigsp, Pb_I_bond_my_sigsp, Pb_I_bond_mz_sigsp]
ppsig = [Pb_I_bond_x_sigpp, Pb_I_bond_y_sigpp, Pb_I_bond_z_sigpp,Pb_I_bond_mx_sigpp, Pb_I_bond_my_sigpp, Pb_I_bond_mz_sigpp]
pppi = [Pb_I_bond_x_pi_y, Pb_I_bond_x_pi_z, Pb_I_bond_y_pi_x, Pb_I_bond_y_pi_z, Pb_I_bond_z_pi_x, Pb_I_bond_z_pi_y,Pb_I_bond_mx_pi_y, Pb_I_bond_mx_pi_z, Pb_I_bond_my_pi_x, Pb_I_bond_my_pi_z, Pb_I_bond_mz_pi_x, Pb_I_bond_mz_pi_y]
pbs_onsites = [Pb_onsite_s]
pbp_onsites = [Pb_onsite_px,Pb_onsite_py, Pb_onsite_pz]
ip_onsites = [I_x_onsite_px, I_x_onsite_py,I_x_onsite_pz, I_y_onsite_px, I_y_onsite_py, I_y_onsite_pz, I_z_onsite_px, I_z_onsite_py, I_z_onsite_pz]
nn_ii_sig = [I_I_nn_sig_x, I_I_nn_sig_y, I_I_nn_sig_z]
nn_pbpb_sig = [Pb_Pb_nn_sig_x, Pb_Pb_nn_sig_y, Pb_Pb_nn_sig_z]
nn_ii_cage = [I_I_nn_x_y, I_I_nn_px_y, I_I_nn_x_my, I_I_nn_px_my,I_I_nn_x_z, I_I_nn_px_z, I_I_nn_x_mz, I_I_nn_px_mz,I_I_nn_y_z, I_I_nn_py_z, I_I_nn_y_mz, I_I_nn_py_mz]

def write_param_file(H,hmap, n_vec = np.array([2,2,2]), ortho_cell = False, N_ortho = None,param_sets = {}):
	#param_sets = {"filename":[param1,param2, ...]}
	lines = {file:[] for file in param_sets}
	
	for file in param_sets:
		parameters = param_sets[file]
		for spin in ["up","down"]:
			for ix in range(n_vec[0]):
				for iy in range(n_vec[1]):
					for iz in range(n_vec[2]):
						counter = 1
						for param in parameters:
							#lines[file].append(" ".join([str(i) for i in param(np.array([ix,iy,iz]), n_vec, atom_positions,box,tilt, coul = coul, vdw = vdw, atom_potentials = atom_potentials)]))
							p,n1,n2 = param(H,hmap,ix,iy,iz,spin)
							cell = str(ix) + str(iy) + str(iz)

							#lines[file].append(" ".join([str(i) for i in [cell,counter,p.real,p.imag]]))
							lines[file].append(" ".join([str(i) for i in [n1,n2,p.real,p.imag]]))
	for filename in param_sets:
		f = open(filename, "w")
		f.write("\n".join(lines[filename]))# apparently .join is the most efficient way to write large strings

#def pbpbhopping_