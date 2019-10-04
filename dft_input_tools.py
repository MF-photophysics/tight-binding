import numpy as np
import matplotlib.pyplot as plt 
import lammps_structure_tools as structure_tools
import lammps_analysis_tools as analysis_tools

lat_digits = 6

def write_all(dump, name, timestep, nkpts = 4, ibrav = 0, ecutwfc = 50,folder = "",write_js = True, N = [4,1/8], n = 64, t = ["8:00:00", "32:00:00"] ): 
	write_dft_input(dump, name + "_scf.in",timestep, title = name, ibrav = ibrav, calculation = "scf", nkpts = nkpts, ecutwfc = ecutwfc,folder = folder)
	write_dft_input(dump, name + "_nscf.in",timestep, title = name, ibrav = ibrav, calculation = "bands", nkpts = nkpts, ecutwfc = ecutwfc,folder = folder)
	write_wannier(dump, name + "_map.win", timestep, title = name , nkpts = nkpts)
	write_pw2wan(name + "_map.pw2wan", folder = folder,seedname = name + "_map", outdir = name + "_out")
	if write_js:
		write_dft_js(name + "_dft.sh", name + "_scf.in", name + "_nscf.in",N = N[0], n = n, log_prefix = name, queue = "regular", t = t[0])
		write_wannier_js(name + "_wan.sh", name + "_map.win", name + "_map.pw2wan", N = N[1], log_prefix = name, queue = "regular" if N[1] >= 1 else "shared", t = t[1])


def write_dft_input(dump,file,timestep,title = "", ecutwfc = 50, ibrav = 0,pseudo_dir = "/global/homes/d/d_abram/dft/pseudos",folder = "", calculation = "scf", nkpts = 4): 
	step = dump.timesteps[timestep]
	input_lines = []
	input_lines.append("&CONTROL")
	input_lines.append("calculation = '" + calculation + "'")
	#input_lines.append("MAP_" + title)
	input_lines.append("title = '" + title + "'" )
	input_lines.append("pseudo_dir = '" + str(pseudo_dir) + "'") #'../pseudos'")
	input_lines.append("etot_conv_thr = 1.0D-3")
	input_lines.append("outdir = '/global/cscratch1/sd/d_abram/" + folder +"/" +title + "_out'")
	input_lines.append("/")

	input_lines.append("&SYSTEM")
	if ibrav == -12: 
		A,B,C = step.box_vect[0], step.box_vect[1], (step.box_vect[2]**2 + step.xz**2)**(1/2)
		cosac = step.xz / A # = a dot c / (a*c) 
		input_lines.append("ibrav = " + str(ibrav))
		input_lines.append("A = " + str(A))
		input_lines.append("B = " + str(B))
		input_lines.append("C = " + str(C))
		input_lines.append("cosAC = " + str(cosac))
	elif ibrav == 8: 
		A,B,C = step.box_vect
		input_lines.append("ibrav = " + str(ibrav))
		input_lines.append("A = " + str(A))
		input_lines.append("B = " + str(B))
		input_lines.append("C = " + str(C))
	elif ibrav == 0: 
		input_lines.append("ibrav = " + str(ibrav))
		x,y,z = dump.timesteps[timestep].box_vect
		if hasattr(dump.timesteps[timestep],'xy'): 
			xy,xz,yz = dump.timesteps[timestep].xy,dump.timesteps[timestep].xz,dump.timesteps[timestep].yz
		else:
			xy,xz,yz = 0.0,0.0,0.0
		
	#input_lines.append("input_dft = \"PBE\"") can be read from pseudopotentials
	input_lines.append("ecutwfc = " + str(ecutwfc))
	natoms = len(dump.structure.atoms)
	input_lines.append("nat = " + str(natoms))
	input_lines.append("ntyp = 5")
	input_lines.append("nbnd = " + str(round(350 * (natoms / 48))))
	input_lines.append("lspinorb = .true.")
	input_lines.append("noncolin = .true.")
	input_lines.append("occupations = 'fixed'")
	input_lines.append("/")
	input_lines.append("")



	input_lines.append("&ELECTRONS")
	input_lines.append("mixing_mode = 'plain'")
	input_lines.append("conv_thr = 1.0D-4")
	input_lines.append("/")

	if ibrav == 0: 
		x,y,z,xy,xz,yz = [round(i,lat_digits) for i in [x,y,z,xy,xz,yz]]
		input_lines.append("CELL_PARAMETERS angstrom")
		input_lines.append(str(x) + " 0.0 0.0")
		input_lines.append(str(xy) + " " + str(y) + " 0.0")
		input_lines.append(str(xz) + " " + str(yz) + " " + str(z))
		#input_lines.append("/")

	input_lines.append("")
	input_lines.append("ATOMIC_SPECIES")
	input_lines.append(" H 1.008 h_srl_gga.upf")
	input_lines.append(" C 12.01 c_srl_gga.upf")
	input_lines.append(" N 14.01 n_srl_gga.upf")
	input_lines.append(" Pb 207.2 pb_frl_gga.upf")
	input_lines.append(" I 126.9 i_frl_gga.upf")
	input_lines.append("")
	types_to_species = {1: "H", 2: "H", 3: "N" , 4: "C", 5: "Pb", 6: "I"}

	input_lines.append("ATOMIC_POSITIONS angstrom") # include a way to do crystal coordinates too
	for a in dump.structure.atoms: 
		species = types_to_species[dump.structure.atom_types[a.atom_type]]
		x, y, z = dump.atom_pos(dump.structure.atoms[a], timestep)
		input_lines.append( " " + species + " " + str(x) + " " + str(y) + " " + str(z))
	input_lines.append("")


	if calculation == "scf":
		input_lines.append("K_POINTS automatic") # I need to figure out what setting for this
		input_lines.append(" " + str(nkpts) + " " + str(nkpts) + " " + str(nkpts) + " 0 0 0") # or 1 1 1 ? 
		input_lines.append("")
	elif calculation == "bands": 
		input_lines.append("K_POINTS (crystal)")
		input_lines.append(str(nkpts**3))
		for nx in range(nkpts): 
			for ny in range(nkpts):
				for nz in range(nkpts):
					kx,ky,kz = nx/nkpts, ny/nkpts, nz/nkpts
					input_lines.append(str(kx) + " " + str(ky) + " " + str(kz) + " 1.0")


	f = open(file, "w")
	f.write("\n".join(input_lines))# apparently .join is the most efficient way to write large strings


def write_wannier(dump, file, timestep, title = "", nkpts = 4): 
	input_lines = []
	natoms = len(dump.structure.atoms)
	num_wann = 152 * round(natoms/48)
	num_bands = 350 * round(natoms/48)
	input_lines.append("num_wann  = " + str(num_wann))
	input_lines.append("num_bands = " + str(num_bands))
	input_lines.append("num_iter = 0")
	input_lines.append("conv_window = 35")
	input_lines.append("conv_tol = 1.e-7")
	input_lines.append("spinors = .true.")
	input_lines.append("")
	input_lines.append("dis_win_max = 15.0")
	input_lines.append("dis_win_min = -10.0")
	input_lines.append("dis_froz_max = 9.2")
	input_lines.append("dis_froz_min = -5.0")
	input_lines.append("dis_num_iter = 29999")
	input_lines.append("dis_mix_ratio = 0.2")
	input_lines.append("dis_conv_tol = 4.0e-7")
	input_lines.append("dis_conv_window = 35")

	#input_lines.append("restart = plot")
	input_lines.append("write_hr = T")
	input_lines.append("begin unit_cell_cart")
	input_lines.append("ang")
	x,y,z = dump.timesteps[timestep].box_vect
	if hasattr(dump.timesteps[timestep],'xy'): 
		xy,xz,yz = dump.timesteps[timestep].xy,dump.timesteps[timestep].xz,dump.timesteps[timestep].yz
	else:
		xy,xz,yz = 0.0,0.0,0.0
	x,y,z,xy,xz,yz = [round(i,lat_digits) for i in [x,y,z,xy,xz,yz]]
	input_lines.append(str(x) + " 0.0 0.0")
	input_lines.append(str(xy) + " " + str(y) + " 0.0")
	input_lines.append(str(xz) + " " + str(yz) + " " + str(z))
	input_lines.append("end unit_cell_cart")

	input_lines.append("")
	input_lines.append("begin projections")
	input_lines.append("Pb:p")
	input_lines.append("I:p")
	input_lines.append("Pb:s")
	input_lines.append("C:p")
	input_lines.append("N:p")
	input_lines.append("end projections")
	input_lines.append("")

	input_lines.append("mp_grid : " + str(nkpts) + " " + str(nkpts) + " " + str(nkpts))
	input_lines.append("begin kpoints")
	for nx in range(nkpts): 
		for ny in range(nkpts):
			for nz in range(nkpts):
				kx,ky,kz = nx/nkpts, ny/nkpts, nz/nkpts
				input_lines.append(str(kx) + " " + str(ky) + " " + str(kz) + " 1.0")
	input_lines.append("end kpoints")

	# k points

	a_lat = x
	#a1,a2,a3 = np.array([a1,0,0]),np.array([0,a2,0]),np.array([0,0,a3])
	#kx,ky,kz = np.array([(2*np.pi*a_lat) / x,0,0]), np.array([0,2*np.pi*a_lat / y,0]), np.array([0,0,2*np.pi *a_lat/z])

	input_lines.append("")
	types_to_species = {1: "H", 2: "H", 3: "N" , 4: "C", 5: "Pb", 6: "I"}
	input_lines.append("begin atoms_cart")
	input_lines.append("ang")
	for a in dump.structure.atoms: 
		species = types_to_species[dump.structure.atom_types[a.atom_type]]
		x, y, z = dump.atom_pos(dump.structure.atoms[a], timestep)
		input_lines.append(species + " " + str(x) + " " + str(y) + " " + str(z))
	input_lines.append("end atoms_cart")
	f = open(file, "w")
	f.write("\n".join(input_lines))# apparently .join is the most efficient way to write large strings

def write_pw2wan(file, folder = "", prefix = "pwscf", seedname = "map", outdir = "scf_out"): 
	input_lines = []
	input_lines.append("&inputpp")
	input_lines.append("  outdir = '/global/cscratch1/sd/d_abram/" + folder + "/" + outdir + "'")
	input_lines.append("  prefix = '" + prefix + "'")
	input_lines.append("  seedname = '" + seedname + "'")
	input_lines.append("  spin_component = 'none'")
	input_lines.append("  write_mmn = .true.")
	input_lines.append("  write_amn = .true.")
	input_lines.append("/")
	f = open(file, "w")
	f.write("\n".join(input_lines))

def write_dft_js(file, scf_file, nscf_file, N = 4, n = 64, log_prefix = "",queue = "regular", t = "16:00:00",name = None): 
	input_lines = []
	input_lines.append("#!/bin/bash")
	input_lines.append("#SBATCH -N " + str(N))
	input_lines.append("#SBATCH -C haswell")
	input_lines.append("#SBATCH -q " + queue)
	if name is None: 
		name = file.split(".in")[0]
	input_lines.append("#SBATCH -J " + name)
	input_lines.append("#SBATCH -t " + t)
	input_lines.append("")
	input_lines.append("#OpenMP settings:")
	input_lines.append("#export OMP_NUM_THREADS=1")
	input_lines.append("#export OMP_PLACES=threads")
	input_lines.append("#export OMP_PROC_BIND=spread")
	input_lines.append("")
	input_lines.append("#run the application")
	input_lines.append("module load espresso")
	c = int(round(N*64/n))
	input_lines.append("srun -n " + str(n) + " -c " + str(c) + " --cpu_bind=cores " + "pw.x " + "-i " + str(scf_file) + " -npool 8 > " + log_prefix + ("_scf.log" if log_prefix else "scf.log"))
	input_lines.append("srun -n " + str(n) + " -c " + str(c) + " --cpu_bind=cores " + "pw.x " + "-i " + str(nscf_file) + " -npool 8 > " + log_prefix + ("_nscf.log" if log_prefix else "nscf.log"))
	f = open(file,'w')
	f.write("\n".join(input_lines))

def write_wannier_js(file, win_file, pw2wan_file, N = 1, log_prefix = "", queue = "regular", t = "16:00:00", name = None):
	input_lines = []
	input_lines.append("#!/bin/bash")
	if N >= 1:
		input_lines.append("#SBATCH -N " + str(N))
		n = 64
	else: 
		n = round(N * 32)
		input_lines.append("#SBATCH --ntasks=" + str(n))
		input_lines.append("#SBATCH --cpus-per-task=2")
		queue = "shared"
	input_lines.append("#SBATCH -C haswell")
	input_lines.append("#SBATCH -q " + queue)
	if name is None: 
		name = file.split(".in")[0]
	input_lines.append("#SBATCH -J " + name)
	input_lines.append("#SBATCH -t " + t)
	input_lines.append("")
	input_lines.append("#OpenMP settings:")
	input_lines.append("#export OMP_NUM_THREADS=1")
	input_lines.append("#export OMP_PLACES=threads")
	input_lines.append("#export OMP_PROC_BIND=spread")
	input_lines.append("")
	input_lines.append("#run the application")
	input_lines.append("module load espresso")
	input_lines.append("module load wannier90")
	input_lines.append("srun -n "  + str(1) + " wannier90.x -pp " + win_file + " > " + log_prefix + ("_wan_pre.log" if log_prefix else "wan_pre.log"))
	input_lines.append("srun -n "  + str(n) + " pw2wannier90.x < " + pw2wan_file + " > " + log_prefix + ("_pw2wan.log" if log_prefix else "pw2wan.log"))
	input_lines.append("srun -n "  + str(1) + " wannier90.x " + win_file + " > " + log_prefix + ("_wan.log" if log_prefix else "wan.log"))
	f = open(file, 'w')
	f.write("\n".join(input_lines))
