import numpy as np
import lammps_analysis_tools as analysis_tools
from index_conversions import * 
#import random

def param_list_to_str(params): 
	string = ""
	for i in params: 
		string += str(i) + " "
	return string

#def rand_velocity_for_temp(atom, T): 

class atom(): 
	def __init__(self,atom_type,molecule_tag = None, charge = "default",position = None,hpath = []):
		self.atom_type = atom_type
		#self.type_N = self.atom_type.number
		if charge == "default": 
			self.charge = atom_type.default_charge
		else: 
			self.charge = charge
		self.molecule_tag = molecule_tag
		self.position = position
		self.hpath = hpath

class atom_type(): 
	def __init__(self,mass,default_charge = 0): 
		self.mass = mass
		self.default_charge = default_charge

class pair_type(): 
	def __init__(self,style, *style_args): 
		self.style = style
		self.style_args = style_args
		self.pairs = []

	def param_str(self): 
		return param_list_to_str(self.style_args)

	def add_pair(self,pair): 
		self.pairs.append(pair)

class pair(): 
	def __init__(self,atom_type1,atom_type2, *args): 
		self.atom_type1 = atom_type1
		self.atom_type2 = atom_type2
		self.args = args

	def param_str(self): 
		return param_list_to_str(self.args)

class bond_type(): 
	def __init__(self,style,*params): 
		self.style = style
		self.params = params

	def param_str(self): 
		return param_list_to_str(self.params)

class bond():
	def __init__(self,bond_type,atom1,atom2):
		self.bond_type = bond_type 
		self.atom1 = atom1
		self.atom2 = atom2

class angle_type(): 
	def __init__(self,style,*params): 
		self.style = style
		self.params = params

	def param_str(self): 
		return param_list_to_str(self.params)

class angle():
	def __init__(self,angle_type,atom1,atom2,atom3):
		self.angle_type = angle_type 
		self.atom1 = atom1
		self.atom2 = atom2
		self.atom3 = atom3

class dihedral_type(): 
	def __init__(self,style,*params): 
		self.style = style
		self.params = params

	def param_str(self): 
		return param_list_to_str(self.params)

class dihedral():
	def __init__(self,dihedral_type,atom1,atom2,atom3,atom4):
		self.dihedral_type = dihedral_type 
		self.atom1 = atom1
		self.atom2 = atom2
		self.atom3 = atom3
		self.atom4 = atom4

class molecule(): 
	
	def __init__(self, atoms, bonds, angles, dihedrals, com = np.array([0,0,0]), molecule_tag = None, handles = {},hpath =[]):
		self.atoms = atoms
		if molecule_tag: 
			for a in self.atoms: 
				a.molecule_tag = molecule_tag
		init_com = 0
		total_mass = 0
		for a in self.atoms: 
			init_com += a.atom_type.mass * a.position
			total_mass += a.atom_type.mass
		init_com = init_com / total_mass
		for a in self.atoms: 
			a.position = a.position - init_com + com
		self.bonds = bonds
		self.angles = angles
		self.dihedrals = dihedrals
		self.handles = handles
		self.hpath = hpath
		for a in self.atoms:
			a.hpath = self.hpath + a.hpath

	def translate(self, r): 
		for a in self.atoms: 
			a.position += r


class unit_cell(): 

	def __init__(self, atoms = [], molecules = [],handles = {},hpath = []): 
		self.atoms = atoms
		self.molecules = molecules
		self.handles = handles
		self.hpath = hpath
		for a in self.atoms:
			a.hpath = self.hpath + a.hpath
		for m in self.molecules: 
			m.hpath = self.hpath + m.hpath


	def translate(self,r): 
		for a in self.atoms: 
			a.position += r
		for m in self.molecules: 
			m.translate(r)


"""
class lattice(): 
	def __init__(self, atoms, v1,v2,v3): 
"""

class structure(): 
	def __init__(self, atom_types, bond_types, angle_types, dihedral_types,pair_repulsions, xlo,xhi,ylo,yhi,zlo,zhi, dielectric = 1, xz_tilt = 0):
		self.atom_types = {atom_types[n]:n + 1 for n in range(len(atom_types))} # to match lammps indexing
		self.atoms = {}
		self.bond_types = {bond_types[n]:n + 1 for n in range(len(bond_types))}
		self.bonds = {}
		self.angle_types = {angle_types[n]:n + 1 for n in range(len(angle_types))}
		self.angles = {}
		self.dihedral_types = {dihedral_types[n]:n + 1 for n in range(len(dihedral_types))}
		self.dihedrals = {}
		self.pair_repulsions = pair_repulsions
		self.dielectric = dielectric
		self.xlo, self.xhi, self.ylo, self.yhi, self.zlo,self.zhi = xlo,xhi,ylo,yhi,zlo,zhi
		self.xz = xz_tilt
		self.handles = {}

	def zero_atom(self,atom_handle = ["0:0:0", "Pb"]):
		# translates all atoms so that the designated atom is at 0
		translation = self.handle_lookup(atom_handle).position
		for a in self.atoms:  
			a.position = a.position - translation

	def translate_to_first_cell(self): 
		a,b,c = np.array([self.xhi - self.xlo,0,0]), np.array([0,self.yhi -self.ylo,0]), np.array([self.xz, 0, self.zhi - self.zlo])
		for atom in self.atoms:
			atom.position = atom.position % (a + b + c)


	def handle_lookup(self,handle_list,output = "object"): 
		# return the atom object or number of an atom, specified by handles
		location = self
		for handle in handle_list: 
			# ie identifier = ["1:2:0","MA2","C"]
			location = location.handles[handle]
		if output == "object": 
			return location
		elif output == "number": 
			return self.atoms[location]
		
	def add_atom(self,atom):
		#type_num = self.atom_types[atom.atom_type]
		self.atoms[atom] = len(self.atoms) + 1

	def add_bond(self,bond): 
		#type_num = self.bond_types[bond.bond_type]
		self.bonds[bond] = len(self.bonds) + 1

	def add_angle(self,angle): 
		#type_num = self.angle_types[angle.angle_type]
		self.angles[angle] = len(self.angles) + 1

	def add_dihedral(self,dihedral): 
		#type_num = self.dihedral_types[dihedral.type]
		self.dihedrals[dihedral] = len(self.dihedrals) + 1

	def add_molecule(self, molecule,handle = None): 
		for a in molecule.atoms: 
			self.add_atom(a)
		for b in molecule.bonds: 
			self.add_bond(b)
		for a in molecule.angles: 
			self.add_angle(a)
		for d in molecule.dihedrals: 
			self.add_dihedral(d)
		if handle: 
			self.handles[handle] = molecule

	def add_cell(self,cell,handle = None): 
		for a in cell.atoms: 
			self.add_atom(a)
		for m in cell.molecules: 
			self.add_molecule(m)
		if handle: 
			self.handles[handle] = cell

	def write_pair_repulsions(self,return_str = False, file_name = "pairs.dat", include_style = False): 
		# need to change this to write potential for all interactions and incldue "hybrids"
		string = ""
		if include_style:
			if len(self.pair_repulsions) == 0: 
				string += "pair_style none" + "\n"
			elif len(self.pair_repulsions) == 1: 
				string += "pair_style " + self.pair_repulsions[0].style + " " + self.pair_repulsions[0].param_str() + "\n"
			else: 
				string += "pair_style hybrid "
				for p_style in self.pair_repulsions: 
					string += p_style.style + " " + p_style.param_str()
				string += "\n"
		string += "Pair Coeffs" + "\n"  + "\n"
		for p_style in self.pair_repulsions: 
				#string += "pair_style " + p_style.style + " " + p_style.param_str() + "\n"
				for p in p_style.pairs: 
					atom_num1 = self.atom_types[p.atom_type1]
					atom_num2 = self.atom_types[p.atom_type2]
					string += str(min(atom_num1, atom_num2)) + " " + str(max(atom_num1,atom_num2)) + " " + p_style.style + " " +  p.param_str() + "\n"
		if return_str: 
			return string
		f = open(file_name,'w')
		f.write(string)

	def positions_from_dump(self,dump,timestep, conversion = stoich_to_unit_conversion,coord_transform = "auto" , handle_list = [], scale_vect = np.array([1,1,1])): #6*np.sqrt(2) / 2, 6*2/2, 6*np.sqrt(2)/ 2])):  
		if coord_transform == "auto": 
			a,c = np.array([self.xhi - self.xlo, 0,0]), np.array([self.xz, 0, self.zhi - self.zlo])
			cos = 1/2 * np.linalg.norm(a + c) / np.linalg.norm(a)
			coord_transform = lambda r: ortho_to_stoich_transform(r,cos)
		update_positions_from_dump(self,dump, timestep, conversion, coord_transform, handle_list, scale_vect )

	def write_structure_file(self, file_name = "structure.dat", description = "", molecule_tags = True, include_pairs = False):
		# writes a structure file in the LAMMPS format, optionally including pair potential coefficients. 
		# header
		string = "LAMMPS Structure File" + description + '\n' + '\n' # first two lines ignored 
		string += str(len(self.atoms)) + " atoms" + '\n'
		string += str(len(self.bonds)) + " bonds" + '\n'
		string += str(len(self.angles)) + " angles" + '\n'
		string += str(len(self.dihedrals)) + " dihedrals" + '\n'
		string += "0" + " impropers" + '\n' # can add in impropers if a system has them 
		string += "\n"
		string += str(len(self.atom_types)) + " atom types" + "\n"
		string += str(len(self.bond_types)) + " bond types" + "\n"
		string += str(len(self.angle_types)) + " angle types" + "\n"
		string += str(len(self.dihedral_types)) + " dihedral types" + "\n"
		string += "0" + " improper types" + "\n"
		string += str(self.xlo) + " " + str(self.xhi) + " xlo xhi" + "\n" 
		string += str(self.ylo) + " " + str(self.yhi) + " ylo yhi" + "\n"
		string += str(self.zlo) + " " + str(self.zhi) + " zlo zhi" + "\n"
		if self.xz: 
			string += "0.0 " + str(self.xz) + " 0.0 xy xz yz" + "\n"
		string += "\n"

		# masses
		string += "Masses" + "\n" + "\n"
		for tp in self.atom_types: 
			string += str(self.atom_types[tp]) + " " + str(tp.mass) + '\n'
		string += "\n" 


		# Bond Coeffs 
		string += "Bond Coeffs" + "\n" + "\n"
		for tp in self.bond_types: 
			string += str(self.bond_types[tp]) + " " + tp.param_str() + "\n"
		string += "\n"

		# Angle Coeffs 
		string += "Angle Coeffs" + "\n" + "\n"
		for tp in self.angle_types: 
			string += str(self.angle_types[tp]) + " " + tp.param_str() + "\n"
		string += "\n"

		# Dihedral Coeffs 
		string += "Dihedral Coeffs" + "\n" + "\n"
		for tp in self.dihedral_types: 
			string += str(self.dihedral_types[tp]) + " " + tp.param_str() + "\n"
		string += "\n"

		# I skip some sections that don't show up in my system but they can be added later

		# Atoms
		string += "Atoms" + "\n" + "\n"
		for a in self.atoms: 
			atom_num = self.atoms[a]
			type_num = self.atom_types[a.atom_type]
			q = a.charge
			x, y, z = a.position
			if a.molecule_tag and molecule_tags: 
				mol_tag = a.molecule_tag
				string += str(atom_num) + " " + str(mol_tag) + " " + str(type_num) + " " + str(q) + " " + str(x) + " "  + str(y) + " " + str(z) + "\n" 
				# could add nx,ny, nz if needed
			else: 
				string += str(atom_num) + " " + str(type_num) + " " + str(q) + " " + str(x) + " "  + str(y) + " " + str(z) + "\n" 
		string += "\n"
		# can add velocities if needed

		# Bonds
		string += "Bonds" + "\n" + "\n"
		for b in self.bonds: 
			bond_num = self.bonds[b]
			type_num = self.bond_types[b.bond_type]
			atom_1_num, atom_2_num = self.atoms[b.atom1],self.atoms[b.atom2]
			atom_1_num,atom_2_num = min(atom_1_num,atom_2_num), max(atom_1_num,atom_2_num) # not sure if I need this
			string += str(bond_num) + " " + str(type_num) + " " + str(atom_1_num) + " " + str(atom_2_num) + "\n" 
		string += "\n"

		# angles 
		string += "Angles" + "\n" + "\n" 
		for a in self.angles: 
			angle_num = self.angles[a]
			type_num = self.angle_types[a.angle_type]
			atom_1_num, atom_2_num, atom_3_num = self.atoms[a.atom1],self.atoms[a.atom2], self.atoms[a.atom3]
			string += str(angle_num) + " " + str(type_num) + " " + str(atom_1_num) + " " + str(atom_2_num) + " " + str(atom_3_num) + "\n" 
		string += "\n"

		# dihedrals
		string += "Dihedrals" + "\n" + "\n" 
		for d in self.dihedrals: 
			dihedral_num = self.dihedrals[d]
			type_num = self.dihedral_types[d.dihedral_type]
			atom_1_num, atom_2_num, atom_3_num,atom_4_num = self.atoms[d.atom1],self.atoms[d.atom2], self.atoms[d.atom3], self.atoms[d.atom4]
			string += str(dihedral_num) + " " + str(type_num) + " " + str(atom_1_num) + " " + str(atom_2_num) + " " + str(atom_3_num) + " " + str(atom_4_num) +  "\n" 
		string += "\n"

		# pair repulsions
		if include_pairs: 
			string += self.write_pair_repulsions(return_str = True)
		f = open(file_name,'w')
		f.write(string)



def update_positions_from_dump(set, dump,timestep, conversion = stoich_to_unit_conversion,coord_transform = lambda r: ortho_to_stoich_transform(r, 1/np.sqrt(2)), handle_list = [], scale_vect = np.array([6*np.sqrt(2) / 2, 6*2/2, 6*np.sqrt(2)/ 2])): 
	# stoich to unit means the new structure has a "cubic" stoich cell basis and the old has an orthorhombic unit cell basis
	# set can be structure, unit cell, etc. 

	# use with the structure to get the right coordintate transformation. 
	for handle in set.handles: 
		#print(handle)
		#print(handle_list)
		item = set.handles[handle]
		if isinstance(item,atom):
			item.position = coord_transform(dump.atom_pos(conversion(handle_list + [handle]), timestep, id_mode = "handles") % (dump.timesteps[timestep].box_vect / scale_vect))
		else: 
			update_positions_from_dump(item, dump, timestep, conversion, coord_transform, handle_list + [handle], scale_vect)

# the rest of the program is for MAPbI3 perovskites, but it is also an example for how to use the tools

# atom types 
params = "old"
if params == "old": 
	H_N = atom_type(1.0080, 0.54)
	H_C = atom_type(1.0080, 0.023)
	N = atom_type(14.0100, -1.1)
	C = atom_type(12.0100, 0.771)
	Pb = atom_type(207.2,2.03)
	I = atom_type(126.9, -1.13)
elif params == "new": 
	H_N = atom_type(1.0080, 0.316975)
	H_C = atom_type(1.0080, -0.003565)
	N = atom_type(14.0100, -0.699825)
	C = atom_type(12.0100, 0.460195)
	Pb = atom_type(207.2,1.4012)
	I = atom_type(126.9, -0.7006)

# bond types 
C_H_bond = bond_type("harmonic",677.4/2,1.091) #note, some potentials may be defined differently in DL_POLY, check potential
N_H_bond = bond_type("harmonic",738/2,1.033)
C_N_bond = bond_type("harmonic",587.2/2,1.499)

# angle types 
H_C_H_angle = angle_type("harmonic",78/2,110.74)
H_N_H_angle = angle_type("harmonic",81/2,108.11)
N_C_H_angle = angle_type("harmonic",98/2,107.91)
H_N_C_angle = angle_type("harmonic",92.4/2,110.11)

# dihedral types 
# H_N_C_H_dihedral = dihedral_type("charmm",0.1556,3,0,0.833,0.5) # charmm is lammps' equivalent of cos in DL_Poly
H_N_C_H_dihedral = dihedral_type("charmm",0.1556,3,0,0.0)

# DL_POLY lists the electrostatic scale factor and van der waals scale factor separately, but in lammps they must be the same
# I think electrostatic is probably more important (back of the envelope says it's 81 - 204 / dielectric const. times larger)

# pair styles
buck_cutoff = 10 # figure out later
coul_cutoff = 20 # also later
buck = pair_type("buck/coul/long", buck_cutoff, coul_cutoff)
buck_Pb_Pb = pair(Pb,Pb,70359906.629702,0.131258,0.0)
buck_Pb_I = pair(Pb,I,103496.13301,0.321737, 0.0)
buck_I_I = pair(I,I,22793.338582,0.482217,696.949542)
buck_Pb_N = pair(Pb,N,32690390.937995,0.150947,0.0)
buck_Pb_C = pair(Pb,C,32690390.937995, 0.150947, 0.000000)
buck_I_N = pair(I,N, 112936.714213, 0.342426, 0.000000)
buck_I_C = pair(I,C, 112936.714213, 0.342426, 0.000000)

for p in [buck_Pb_Pb,buck_Pb_I,buck_I_I,buck_Pb_N,buck_Pb_C,buck_I_N,buck_I_C]:
	buck.add_pair(p)

lj_cutoff = 10 # also later
lj = pair_type("lj/cut/coul/long", lj_cutoff,coul_cutoff)
lj_Pb_H_N = pair(Pb,H_N,0.0140,2.26454)
lj_Pb_H_C = pair(Pb,H_C,0.0140, 2.70999)
lj_I_H_N = pair(I,H_N,0.0574, 2.75000)
lj_I_H_C = pair(I,H_C,0.0574, 3.10000)
lj_N_N = pair(N,N,0.1700, 3.25000)
lj_N_H_N = pair(N,H_N,0.0517, 2.15950)
lj_N_C = pair(N,C,0.1364, 3.32480)
lj_N_H_C = pair(N,H_C,0.0517, 2.60500)
lj_H_N_H_N = pair(H_N,H_N,0.0157, 1.06910)
lj_H_N_C = pair(H_N,C,0.0414, 2.23440)

lj_C_C = pair(C,C,0.1094, 3.39970)
lj_C_H_C = pair(C,H_C,0.0414, 2.67980)
lj_H_C_H_C = pair(H_C,H_C,0.0157, 1.96000)

for p in [lj_Pb_H_N, lj_Pb_H_C, lj_I_H_N, lj_I_H_C, lj_N_N, lj_N_H_N, lj_N_C, lj_N_H_C, lj_H_N_H_N, lj_H_N_C, lj_C_C, lj_C_H_C, lj_H_C_H_C]:
	lj.add_pair(p)

# LAMMPS requires a separate potential type for 1-4 atoms interacting via a dihedral
lj_charmm = pair_type("lj/charmm/coul/long", lj_cutoff, lj_cutoff, coul_cutoff,coul_cutoff)
lj_H_N_H_C = pair(H_N,H_C,0.0157, 1.51450)

lj_charmm.add_pair(lj_H_N_H_C)



# methyl-ammonium 
def create_MA_mol(molecule_tag = "ch3nh3", position = np.array([0,0,0]),hpath = [], position_mode = "auto"):
	atoms = [atom(H_N, 0.54,hpath = ["H_N1"]),atom(H_N,0.54,hpath = ["H_N2"]),atom(H_N,0.54,hpath = ["H_N3"]),atom(N,-1.1,hpath = ["N"]),atom(C,0.771,hpath = ["C"]),atom(H_C,0.023,hpath = ["H_C1"]),atom(H_C,0.023,hpath =["H_C2"]),atom(H_C,0.023,hpath = ["H_C3"])]
	H_C1, H_C2, H_C3 = atoms[5], atoms[6],atoms[7]
	H_N1, H_N2, H_N3 = atoms[0],atoms[1],atoms[2]
	C1 = atoms[4]
	N1 = atoms[3]
	bonds = [bond(C_H_bond,C1,H_C1), bond(C_H_bond,C1,H_C2),bond(C_H_bond,C1,H_C3),bond(C_N_bond,N1,C1),bond(N_H_bond,N1,H_N1),bond(N_H_bond,N1,H_N2),bond(N_H_bond,N1,H_N3)]
	angles = [angle(H_C_H_angle,H_C1,C1,H_C2),angle(H_C_H_angle,H_C1,C1,H_C3),angle(H_C_H_angle,H_C2,C1,H_C3),angle(H_N_H_angle,H_N1,N1,H_N2),angle(H_N_H_angle,H_N1,N1,H_N3),angle(H_N_H_angle,H_N2,N1,H_N3),angle(N_C_H_angle,N1,C1,H_C1),angle(N_C_H_angle,N1,C1,H_C2),angle(N_C_H_angle,N1,C1,H_C3),angle(H_N_C_angle,H_N1,N1,C1),angle(H_N_C_angle,H_N2,N1,C1),angle(H_N_C_angle,H_N3,N1,C1)]
	dihedrals = [dihedral(H_N_C_H_dihedral, H_N1,N1,C1,H_C1),dihedral(H_N_C_H_dihedral, H_N2,N1,C1,H_C1),dihedral(H_N_C_H_dihedral, H_N3,N1,C1,H_C1),dihedral(H_N_C_H_dihedral, H_N1,N1,C1,H_C2),dihedral(H_N_C_H_dihedral, H_N2,N1,C1,H_C2),dihedral(H_N_C_H_dihedral, H_N3,N1,C1,H_C2),dihedral(H_N_C_H_dihedral, H_N1,N1,C1,H_C3),dihedral(H_N_C_H_dihedral, H_N2,N1,C1,H_C3),dihedral(H_N_C_H_dihedral, H_N3,N1,C1,H_C3)]
	# these positions are kind of rough, can be corrected when put into a crystal 
	H_N1.position = np.array([0.93322,0,1.81047])
	H_N2.position = np.array([-0.46661,0.80819,1.8105])
	H_N3.position = np.array([-0.46661,-0.80819,1.8105])
	N1.position = np.array([0,0,1.48])
	C1.position = np.array([0,0,0])
	H_C1.position = np.array([0.51374,0.88982,-0.36385])
	H_C2.position = np.array([-1.0275,0,-0.36385])
	H_C3.position = np.array([0.51374,-0.88982,-0.36385])
	rot = np.array([[1/2**0.5, -1/2**0.5, 0],[1/2**0.5, 1/2**0.5, 0],[0,0,1]])
	rot2 = np.array([[0,1,0],[-1,0,0],[0,0,1]])
	defposfunc = lambda r :rot @  (lambda r: np.array([r[1], r[2], r[0]]))(r)
	if position_mode == "rand2":
		if np.random.random() < 0.5:
			posfunc = lambda r : -defposfunc(r)
		else: 
			posfunc = lambda r : defposfunc(r)
		H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position = [posfunc(r) for r in [H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position]]
	elif position_mode == "ortho1":
		posfunc = lambda r : -defposfunc(r)
		H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position = [posfunc(r) for r in [H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position]]
	elif position_mode == "ortho2": 
		posfunc = lambda r: defposfunc(r)
		H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position = [posfunc(r) for r in [H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position]]
	elif position_mode == "rand4": 
		rand = np.random.random()
		if rand < 0.25: 
			posfunc = lambda r : - rot2 @ defposfunc(r)
		elif rand < 0.5: 
			posfunc = lambda r : rot2 @ defposfunc(r)
		elif rand < 0.75: 
			posfunc = lambda r : - defposfunc(r)
		else: 
			posfunc = defposfunc
		H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position = [posfunc(r) for r in [H_N1.position, H_N2.position, H_N3.position, N1.position, C1.position, H_C1.position, H_C2.position, H_C3.position]]
	atom_handles = {"H_C1": H_C1, "H_C2": H_C2, "H_C3": H_C3, "H_N1":H_N1, "H_N2": H_N2, "H_N3": H_N3, "C": C1, "N": N1}
	ch3nh3 = molecule(atoms,bonds,angles, dihedrals, position, molecule_tag, atom_handles,hpath = hpath)
	return ch3nh3


#lx,ly,lz = 8.557,9.250,12.964 # unit cell dimensions, from materials project for now determine later

lx,ly,lz = 8.557, 12.964, 9.250
V = lx*ly*lz
structure_style = "orthorhombic"
if structure_style == "tetragonal": 
	ratio = 1.06
	lx = (V / ratio)**(1/3)
	lz,ly = lx,lx * ratio
elif structure_style == "square":
	lx = V**(1/3)
	ly,lz = lx,lx
l = (V/4)**(1/3)

def create_MAPbI3_cubic_cell(nx,ny,nz, lx = l, ly = l, lz = l, xz_cell = 0.0,hpath = [], randMA = False): 
	# I am fairly sure this structure is wrong
	# I fixed it, roughly. It's not very close to equilibrium but should be okay at high temperatures. 
	if not xz_cell:
		translation = np.array([nx*lx, ny * ly, nz * lz]) 
		cell_vect = np.array([lx,ly,lz])
	else: 
		xz_cell = 0.0 if xz_cell is True else xz_cell
		lxz = np.sqrt(lx**2 + lz**2) /2 
		costheta = (lx**2 - lz**2) / (lx**2 + lz**2)#(lx*lz)
		sintheta = (1 - costheta**2)**(1/2)
		a,b,c = np.array([lxz,0,0]), np.array([0,ly/2,0]), np.array([xz_cell, 0,lxz * sintheta ])
		cell_vect = a + b + c 
		translation = nx * a + ny * b + nz * c
	cell_str = str(nx) + "," + str(ny) + "," + str(nz)
	Pb_atom = atom(Pb,molecule_tag = "Pb:" + cell_str, position = np.array([0.0,0.0,0.0]) * cell_vect,hpath = ["Pb"])
	I_atom1 = atom(I, molecule_tag = "I:" + cell_str, position = np.array([0.5,0.0,0.0]) * cell_vect,hpath = ["I1"])
	I_atom2 = atom(I, molecule_tag = "I:" + cell_str, position = np.array([0.0,0.5 ,0.0]) * cell_vect,hpath = ["I2"])
	I_atom3 = atom(I,molecule_tag = "I:" + cell_str, position = np.array([0.0,0.0,0.5]) * cell_vect,hpath = ["I3"])
	if randMA:
		MA_molecule = create_MA_mol(molecule_tag = "MA:" + cell_str, position = np.array([0.5,0.5,0.5])*cell_vect,hpath = ["MA"],position_mode = "rand2" if xz_cell else "rand4")
	else: 
		MA_molecule = create_MA_mol(molecule_tag = "MA:" + cell_str, position = np.array([0.5,0.5,0.5])*cell_vect,hpath = ["MA"],position_mode = "ortho1" if (nx + ny + nz)%2 == 0 else "ortho2")
	handles = {"Pb": Pb_atom, "I1": I_atom1, "I2": I_atom2,"I3":I_atom3,"MA":MA_molecule}
	cell = unit_cell([Pb_atom,I_atom1,I_atom2,I_atom3],[MA_molecule],handles = handles,hpath = hpath)
	cell.translate(translation)
	return cell


def create_MAPbI3_supercell(nx,ny,nz,lx = lx, ly = ly, lz = lz,hpath = []): 
	# structure from materials project https://materialsproject.org/materials/mp-995214/
	translation = np.array([nx*lx, ny * ly, nz * lz])
	cell_str = str(nx) + "," + str(ny) + "," + str(nz)
	scale_vect = np.array([lx,ly,lz])
	Pb_atom1 = atom(Pb,molecule_tag = "Pb1:" + cell_str, position = np.array([0.0,0.0,0.0]) * scale_vect,hpath = ["Pb1"])
	Pb_atom2 = atom(Pb, molecule_tag = "Pb2:" + cell_str, position = np.array([0.0,0.5,0.0]) * scale_vect,hpath = ["Pb2"])
	Pb_atom3 = atom(Pb, molecule_tag = "Pb3:" + cell_str, position = np.array([0.5,0.0,0.5 ])*scale_vect,hpath = ["Pb3"])
	Pb_atom4 = atom(Pb, molecule_tag = "Pb4:" + cell_str, position = np.array([0.5,0.5,0.5]) * scale_vect,hpath = ["Pb4"])
	I_atom1 = atom(I, molecule_tag = "I1:" + cell_str, position = np.array([0.0377,0.25,0.9874]) * scale_vect,hpath = ["I1"])
	I_atom2 = atom(I, molecule_tag = "I2:" + cell_str, position = np.array([0.1839,0.0221,0.3072]) * scale_vect,hpath = ["I2"])
	I_atom3 = atom(I,molecule_tag = "I3:" + cell_str, position = np.array([0.1839,0.4779,0.3072]) * scale_vect,hpath = ["I3"])
	I_atom4 = atom(I, molecule_tag = "I4:" + cell_str, position = np.array([0.3161,0.9779,0.8072]) * scale_vect,hpath = ["I4"])
	I_atom5 = atom(I, molecule_tag = "I5:" + cell_str, position = np.array([0.3161,0.5221,0.8072]) * scale_vect,hpath = ["I5"])
	I_atom6 = atom(I,molecule_tag = "I6:" + cell_str, position = np.array([0.4623,0.75,0.4874]) * scale_vect,hpath = ["I6"])
	I_atom7 = atom(I, molecule_tag = "I7:" + cell_str, position = np.array([0.5377,0.25,0.5126]) * scale_vect,hpath = ["I7"])
	I_atom8 = atom(I, molecule_tag = "I8:" + cell_str, position = np.array([0.6839,0.0221,0.1928]) * scale_vect,hpath = ["I8"])
	I_atom9 = atom(I,molecule_tag = "I9:" + cell_str, position = np.array([0.6839,0.4779,0.1928 ]) * scale_vect,hpath = ["I9"])
	I_atom10 = atom(I, molecule_tag = "I10:" + cell_str, position = np.array([0.8161,0.9779,0.6928]) * scale_vect,hpath = ["I10"])
	I_atom11 = atom(I, molecule_tag = "I11:" + cell_str, position = np.array([0.8161,0.5221,0.6928]) * scale_vect,hpath = ["I11"])
	I_atom12 = atom(I,molecule_tag = "I12:" + cell_str, position = np.array([0.9623,0.75,0.0126]) * scale_vect,hpath = ["I12"])

	MA_molecule1 = create_MA_mol(molecule_tag = "MA1:" + cell_str, position = np.array([0,0,0]), hpath = ["MA1"])
	MA_molecule1.handles["C"].position = np.array([0.0324 +1,0.25,0.5622]) * scale_vect
	MA_molecule1.handles["H_C1"].position = np.array([0.1598 +1,0.25,0.5503]) * scale_vect
	MA_molecule1.handles["H_C2"].position = np.array([0.9948,0.3192,0.6201]) * scale_vect
	MA_molecule1.handles["H_C3"].position = np.array([0.9948,0.1808,0.6201]) * scale_vect
	MA_molecule1.handles["N"].position = np.array([0.9596,0.25,0.4151]) * scale_vect
	MA_molecule1.handles["H_N1"].position = np.array([0.9943,0.3146,0.357]) * scale_vect
	MA_molecule1.handles["H_N2"].position = np.array([0.9943,0.1854,0.357]) * scale_vect
	MA_molecule1.handles["H_N3"].position = np.array([0.8381,0.25,0.4214]) * scale_vect


	MA_molecule2 = create_MA_mol(molecule_tag = "MA2:" + cell_str, position = np.array([0,0,0]),hpath = ["MA2"])
	MA_molecule2.handles["C"].position = np.array([0.4676,0.75,0.0622 + 1]) * scale_vect
	MA_molecule2.handles["H_C1"].position = np.array([0.3402,0.75,0.0503 +1]) * scale_vect
	MA_molecule2.handles["H_C2"].position = np.array([0.5052,0.6808,0.1201 +1]) * scale_vect
	MA_molecule2.handles["H_C3"].position = np.array([0.5052,0.8192,0.1201+1]) * scale_vect
	MA_molecule2.handles["N"].position = np.array([0.5404,0.75,0.9151]) * scale_vect
	MA_molecule2.handles["H_N1"].position = np.array([0.6619,0.75,0.9214]) * scale_vect
	MA_molecule2.handles["H_N2"].position = np.array([0.5057,0.8146,0.857]) * scale_vect
	MA_molecule2.handles["H_N3"].position = np.array([0.5057,0.6854,0.857]) * scale_vect
	

	MA_molecule3 = create_MA_mol(molecule_tag = "MA3:" + cell_str, position = np.array([0,0,0]),hpath = ["MA3"])
	MA_molecule3.handles["C"].position = np.array([0.5324,0.25,0.9378]) * scale_vect
	MA_molecule3.handles["H_C1"].position = np.array([0.6598,0.25,0.9497]) * scale_vect
	MA_molecule3.handles["H_C2"].position = np.array([0.4948,0.1808,0.8799]) * scale_vect
	MA_molecule3.handles["H_C3"].position = np.array([0.4948,0.3192,0.8799]) * scale_vect
	MA_molecule3.handles["N"].position = np.array([0.4596,0.25,0.0849+1]) * scale_vect
	MA_molecule3.handles["H_N1"].position = np.array([0.3381,0.25,0.0786+1]) * scale_vect
	MA_molecule3.handles["H_N2"].position = np.array([0.4943,0.3146,0.143+1]) * scale_vect
	MA_molecule3.handles["H_N3"].position = np.array([0.4943,0.1854,0.143+1]) * scale_vect

	MA_molecule4 = create_MA_mol(molecule_tag = "MA4:" + cell_str, position = np.array([0,0,0]),hpath = ["MA4"])
	MA_molecule4.handles["C"].position = np.array([0.9676,0.75,0.4378]) * scale_vect
	MA_molecule4.handles["H_C1"].position = np.array([0.8402,0.75,0.4497]) * scale_vect
	MA_molecule4.handles["H_C2"].position = np.array([0.0052+1,0.6808,0.3799]) * scale_vect
	MA_molecule4.handles["H_C3"].position = np.array([0.0052+1,0.8192,0.3799]) * scale_vect
	MA_molecule4.handles["N"].position = np.array([0.0404+1,0.75,0.5849]) * scale_vect
	MA_molecule4.handles["H_N1"].position = np.array([0.1619+1,0.75,0.5786]) * scale_vect
	MA_molecule4.handles["H_N2"].position = np.array([0.0057+1,0.6854,0.643]) * scale_vect
	MA_molecule4.handles["H_N3"].position = np.array([0.0057+1, 0.8146, 0.643]) * scale_vect
	
	handles = {"Pb1": Pb_atom1,"Pb2": Pb_atom2, "Pb3":Pb_atom3, "Pb4":Pb_atom4, "I1":I_atom1,"I2":I_atom2,"I3":I_atom3,"I4":I_atom4,"I5":I_atom5,"I6":I_atom6,"I7":I_atom7,"I8":I_atom8,"I9":I_atom9,"I10":I_atom10,"I11":I_atom11,"I12":I_atom12,"MA1": MA_molecule1,"MA2": MA_molecule2,"MA3": MA_molecule3,"MA4": MA_molecule4}
	cell = unit_cell([Pb_atom1,Pb_atom2,Pb_atom3,Pb_atom4,I_atom1,I_atom2,I_atom3,I_atom4,I_atom5,I_atom6,I_atom7,I_atom8,I_atom9,I_atom10,I_atom11,I_atom12],[MA_molecule1,MA_molecule2,MA_molecule3,MA_molecule4],handles = handles,hpath = hpath)
	cell.translate(translation)
	return cell


def create_n_n_n_crystal(cell_func,Nx,Ny = None,Nz = None, lx = lx,ly = ly,lz = lz):
	atom_types = [H_N,H_C,N,C,Pb,I]
	bond_types = [C_H_bond,N_H_bond,C_N_bond]
	angle_types = [H_C_H_angle, H_N_H_angle, N_C_H_angle, H_N_C_angle]
	dihedral_types = [H_N_C_H_dihedral]
	pair_types = [buck,lj,lj_charmm]
	dielectric = 1.0 # or whatever it is
	if Ny == None: 
		Ny = Nx
	if Nz == None: 
		Nz = Nx
	xlo,ylo,zlo = 0.0,0.0,0.0
	xhi,yhi,zhi = lx*Nx,ly*Ny,lz*Nz
	crystal = structure(atom_types, bond_types, angle_types, dihedral_types,pair_types, xlo,xhi,ylo,yhi,zlo,zhi, dielectric = 1)
	for nx in range(Nx): 
		for ny in range(Ny): 
			for nz in range(Nz): 
				#crystal.add_cell(cell_func(nx,ny,nz,lx,ly,lz,hpath = [str(nx) + ":" + str(ny) + ":" + str(nz)]),handle = str(nx) + ":" + str(ny) + ":" + str(nz))
				crystal.add_cell(cell_func(nx,ny,nz,lx,ly,lz,hpath = [(nx,ny,nz)]),handle = (nx,ny,nz))
	return crystal


def create_n_n_n_tilted_crystal(Nx,cell_func = create_MAPbI3_cubic_cell, Ny = None,Nz = None, lx_ortho = lx,ly_ortho = ly,lz_ortho = lz):
	atom_types = [H_N,H_C,N,C,Pb,I]
	bond_types = [C_H_bond,N_H_bond,C_N_bond]
	angle_types = [H_C_H_angle, H_N_H_angle, N_C_H_angle, H_N_C_angle]
	dihedral_types = [H_N_C_H_dihedral]
	pair_types = [buck,lj,lj_charmm]
	if Ny == None: 
		Ny = Nx
	if Nz == None: 
		Nz = Nx
	xlo,ylo,zlo = 0.0,0.0,0.0
	lxz = np.sqrt(lx_ortho**2 + lz_ortho**2)/2
	costheta = (lx_ortho**2 - lz_ortho**2) / (lx_ortho**2 + lz_ortho**2)#(lx_ortho*lz_ortho)
	sintheta = (1 - costheta**2)**(1/2)
	xhi,yhi,zhi = lxz*Nx,ly_ortho*Ny/2,lxz * Nz * sintheta
	xz_cell = lxz * costheta
	crystal = structure(atom_types, bond_types, angle_types, dihedral_types,pair_types, xlo,xhi,ylo,yhi,zlo,zhi, dielectric = 1, xz_tilt = xz_cell * Nz)
	for nx in range(Nx): 
		for ny in range(Ny): 
			for nz in range(Nz): 
				#crystal.add_cell(cell_func(nx,ny,nz,lx_ortho,ly_ortho,lz_ortho,xz_cell,hpath = [str(nx) + ":" + str(ny) + ":" + str(nz)]),handle = str(nx) + ":" + str(ny) + ":" + str(nz))
				crystal.add_cell(cell_func(nx,ny,nz,lx_ortho,ly_ortho,lz_ortho,xz_cell if xz_cell else True,hpath = [(nx,ny,nz)]),handle = (nx,ny,nz))
	return crystal


#crystal = create_n_n_n_crystal(create_MAPbI3_cell,4, halve_lengths = True)
#crystal.write_structure_file(file_name = "MAPbI3_single_n4.dat", include_pairs = False)
MAPbI3_ortho_n1 = create_n_n_n_crystal(create_MAPbI3_supercell,1)
MAPbI3_ortho_n2 = create_n_n_n_crystal(create_MAPbI3_supercell,2)
MAPbI3_ortho_n4 = create_n_n_n_crystal(create_MAPbI3_supercell,4)
MAPbI3_ortho_n6 = create_n_n_n_crystal(create_MAPbI3_supercell,6)
MAPbI3_ortho_n757 = create_n_n_n_crystal(create_MAPbI3_supercell, Nx = 7, Ny = 5, Nz = 7)

l_cubic = (V/4)**(1/3)
MAPbI3_cubic_n4 = create_n_n_n_crystal(create_MAPbI3_cubic_cell,4,lx=l_cubic,ly=l_cubic,lz = l_cubic)
MAPbI3_cubic_n6 = create_n_n_n_crystal(create_MAPbI3_cubic_cell,6,lx = l_cubic, ly = l_cubic, lz = l_cubic)
MAPbI3_cubic_n8 = create_n_n_n_crystal(create_MAPbI3_cubic_cell,8, lx = l_cubic, ly = l_cubic, lz = l_cubic)

tilted_n2 = create_n_n_n_tilted_crystal(2)
#crystal_with_4_cell_symmetry.write_structure_file(file_name = "MAPbI3_n4.dat",include_pairs = False)
#crystal.write_pair_repulsions(file_name = "perovskite_pair_repulsions.dat")

def replicate_cubic(structure, N = 8, n = 2): 
	# primarily for the tilted case
	a = np.array((structure.xhi - structure.xlo,0,0))/N
	b = np.array((0,structure.yhi - structure.ylo,0))/N
	c = np.array([structure.xz, 0, structure.zhi - structure.zlo])/N
	atom_handles = ["Pb", "I1","I2", "I3", "H_C1","H_C2","H_C3", "H_N1", "H_N2", "H_N3", "C", "N"]
	atom_handles = [[i] for i in atom_handles[:4]] + [["MA",i] for i in atom_handles[4:]]
	for nx in range(N): 
		for ny in range(N): 
			for nz in range(N):
				#cell = [ str(nx) + ":" + str(ny) + ":" + str(nz)]
				#base_cell = [str(nx%n) + ":" + str(ny%n) + ":" + str(nz%n)]
				cell = (nx,ny,nz)
				base_cell = (nx%n, ny%n, nz%n)
				xtrans,ytrans, ztrans = n*(nx // n), n*(ny//n), n*(nz //n)
				for at in atom_handles: 
					handle = cell + at
					base_handle = base_cell + at 
					structure.handle_lookup(handle).position = structure.handle_lookup(base_handle).position + xtrans * a + ytrans * b + ztrans * c

