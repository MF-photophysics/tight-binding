import subprocess

input_file_1 = """clear
units real
dimension 3
boundary p p p
atom_style full
atom_modify map array
# define potential
bond_style harmonic
angle_style harmonic
dihedral_style charmm
kspace_style pppm 1.0e-5
kspace_modify diff ad
pair_style hybrid buck/coul/long 10 10 lj/cut/coul/long 10 10 lj/charmm/coul/long 10 10.1 10"""

potentials = """pair_coeff 5 5 buck/coul/long 70359906.629702 0.131258 0.0
pair_coeff 5 6 buck/coul/long 103496.13301 0.321737 0.0
pair_coeff 6 6 buck/coul/long 22793.338582 0.482217 696.949542
pair_coeff 3 5 buck/coul/long 32690390.937995 0.150947 0.0
pair_coeff 4 5 buck/coul/long 32690390.937995 0.150947 0.0
pair_coeff 3 6 buck/coul/long 112936.714213 0.342426 0.0
pair_coeff 4 6 buck/coul/long 112936.714213 0.342426 0.0
pair_coeff 1 5 lj/cut/coul/long 0.014 2.26454
pair_coeff 2 5 lj/cut/coul/long 0.014 2.70999
pair_coeff 1 6 lj/cut/coul/long 0.0574 2.75
pair_coeff 2 6 lj/cut/coul/long 0.0574 3.1
pair_coeff 3 3 lj/cut/coul/long 0.17 3.25
pair_coeff 1 3 lj/cut/coul/long 0.0517 2.1595
pair_coeff 3 4 lj/cut/coul/long 0.1364 3.3248
pair_coeff 2 3 lj/cut/coul/long 0.0517 2.605
pair_coeff 1 1 lj/cut/coul/long 0.0157 1.0691
pair_coeff 1 4 lj/cut/coul/long 0.0414 2.2344
pair_coeff 1 2 lj/charmm/coul/long 0.0157 1.5145
pair_coeff 4 4 lj/cut/coul/long 0.1094 3.3997
pair_coeff 2 4 lj/cut/coul/long 0.0414 2.6798
pair_coeff 2 2 lj/cut/coul/long 0.0157 1.96"""

zero_potentials = """pair_coeff 5 5 buck/coul/long 0.0 0.131258 0.0
pair_coeff 5 6 buck/coul/long 0.0 0.321737 0.0
pair_coeff 6 6 buck/coul/long 0.0 0.482217 0.0
pair_coeff 3 5 buck/coul/long 0.0 0.150947 0.0
pair_coeff 4 5 buck/coul/long 0.0 0.150947 0.0
pair_coeff 3 6 buck/coul/long 0.0 0.342426 0.0
pair_coeff 4 6 buck/coul/long 0.0 0.342426 0.0
pair_coeff 1 5 lj/cut/coul/long 0.0 2.26454
pair_coeff 2 5 lj/cut/coul/long 0.0 2.70999
pair_coeff 1 6 lj/cut/coul/long 0.0 2.75
pair_coeff 2 6 lj/cut/coul/long 0.0 3.1
pair_coeff 3 3 lj/cut/coul/long 0.0 3.25
pair_coeff 1 3 lj/cut/coul/long 0.0 2.1595
pair_coeff 3 4 lj/cut/coul/long 0.0 3.3248
pair_coeff 2 3 lj/cut/coul/long 0.0 2.605
pair_coeff 1 1 lj/cut/coul/long 0.0 1.0691
pair_coeff 1 4 lj/cut/coul/long 0.0 2.2344
pair_coeff 1 2 lj/charmm/coul/long 0.0 1.5145
pair_coeff 4 4 lj/cut/coul/long 0.0 3.3997
pair_coeff 2 4 lj/cut/coul/long 0.0 2.6798
pair_coeff 2 2 lj/cut/coul/long 0.0 1.96"""

timestep_run_input ="""fix 1 all nve
timestep 0.0
run 1
unfix 1"""

neighbor_thermo_compute = """neighbor 3.0 bin
neigh_modify delay 0 every 1 check yes page 500000 one 50000

thermo 100
thermo_style custom dt step pe temp lx ly lz press # pxx pyy pzz #c_eatoms

compute Eksp all pe/atom kspace
compute Epair all pe/atom pair
group PbI type 5 6"""

def restart_potential_run(restart_file,temp,tdamp = 50.0, pdamp = 500.0, ts_start = 0, timesteps = 60000, timestep = 0.5,anneal_time = 300000,pot_style = "oldfiles",  fix = "npt", name = "", mode = "restart",cool_anneal = 2000, relax = False, dumpevery = 60, dump = False,forces = True):
	if pot_style == "oldstrings": 
		real_pot = lambda : input_lines.append(potentials)
		zero_pot = lambda : input_lines.append(zero_potentials)
	elif pot_style == "newfiles": 
		real_pot = lambda : input_lines.append("include new_params.nonbonding")
		zero_pot = lambda : input_lines.append("include zero.nonbonding")
	elif pot_style == "oldfiles": 
		real_pot = lambda : input_lines.append("include params.nonbonding")
		zero_pot = lambda : input_lines.append("include zero.nonbonding")		
	elif pot_style == "NI_change": 
		real_pot = lambda : input_lines.append("include NI_change.nonbonding")
		zero_pot = lambda : input_lines.append("include zero.nonbonding")

	read_restart = "read_restart " + restart_file
	input_lines = []
	input_lines.append(input_file_1)
	if mode == "restart": 
		input_lines.append(read_restart)
	elif mode == "structure":
		input_lines.append("read_data " + restart_file)
		input_lines.append("velocity all create " + str(temp) + " 2093939 dist gaussian mom yes rot yes units box")
	input_lines.append("kspace_style pppm 1.0e-5")
	input_lines.append("reset_timestep " + str(ts_start))
	real_pot()
	input_lines.append(neighbor_thermo_compute)
	if relax == "minimize": 
		input_lines.append("minimize 1.0e-4 1.0e-8 1000 100000")
		input_lines.append("velocity all create " + str(temp) + " 2093939 dist gaussian mom yes rot yes units box")
	group = "group some id 1"
	annealdump = "dump annealdump some custom 10 " + name + "_anneal.dump id x y z type"
	input_lines.append("timestep " + str(timestep))
	input_lines.append(group)
	input_lines.append(annealdump)

	if cool_anneal: 
		max_temp = 300
		if temp < max_temp:
			cool_steps = int((max_temp - temp) * cool_anneal / timestep)
			group = "group some id 1"
			if fix == "npt":
				fixcool = "fix cool all npt temp " + str(max_temp) + " " + str(temp) + " " + str(tdamp) + " aniso 0.0 0.0 " + str(pdamp)
			elif fix == "nvt": 
				fixcool = "fix cool all nvt temp " + str(max_temp) + " " + str(temp) + " " + str(tdamp) 
			elif fix == "npt/iso": 
				fixcool = "fix cool all npt temp " + str(max_temp) + " " + str(temp) + " " + str(tdamp) + " iso 0.0 0.0 " + str(pdamp)
			input_lines.append(fixcool)
			input_lines.append("run " + str(cool_steps))
			input_lines.append("unfix cool")

	if anneal_time > 0: 
		anneal_steps = int(anneal_time/timestep)
		group = "group some id 1"
		annealdump = "dump annealdump some custom 10 " + name + "_anneal.dump id x y z type"
		if fix == "npt":
			fixanneal = "fix anneal all npt temp " + str(temp) + " " + str(temp) + " " + str(tdamp) + " aniso 0.0 0.0 " + str(pdamp)
		elif fix == "nvt": 
			fixanneal = "fix anneal all nvt temp " + str(temp) + " " + str(temp) + " " + str(tdamp) 
		elif fix == "npt/iso": 
			fixanneal = "fix anneal all npt temp " + str(temp) + " " + str(temp) + " " + str(tdamp) + " iso 0.0 0.0 " + str(pdamp)

		input_lines.append(fixanneal)
		input_lines.append("run " + str(anneal_steps))
		input_lines.append("unfix anneal")
		input_lines.append("reset_timestep " + str(ts_start))
	input_lines.append("undump annealdump")
	if forces: 
		rundump = "dump rundump all custom " + str(dumpevery) + " " + name + "_run.dump id x y z c_Epair fx fy fz type "
		couldump = "dump couldump PbI custom " + str(timesteps*2) + " " + name + "_coul.dump id x y z c_Epair fx fy fz type "
	else:
		rundump = "dump rundump all custom " + str(dumpevery) + " " + name + "_run.dump id x y z c_Epair type"
		couldump = "dump couldump PbI custom " + str(timesteps*2) + " " + name + "_coul.dump id x y z c_Epair type"
	# ^ will dump when the modify is used, not otherwise, except at beginning. 
	input_lines.append(rundump)
	input_lines.append(couldump)
	if dump:
		input_lines.append("read_dump " + str(dump) + " 0 x y z box yes")
	for i in range(timesteps//dumpevery):
		if fix == "npt":
			fixnpt = "fix 1 all npt temp " + str(temp) + " " + str(temp) + " " + str(tdamp) + " aniso 0.0 0.0 " + str(pdamp)
		elif fix == "nvt": 
			fixnpt = "fix 1 all nvt temp " + str(temp) + " " + str(temp) + " " + str(tdamp) 
		elif fix == "npt/iso": 
			fixnpt = "fix 1 all npt temp " + str(temp) + " " + str(temp) + " " + str(tdamp) + " iso 0.0 0.0 " + str(pdamp)
		input_lines.append(fixnpt) #may or may not actually be npt
		input_lines.append("timestep " + str(timestep))
		input_lines.append("run " + str(dumpevery))
		input_lines.append("unfix 1")

		zero_pot()
		fixnve = "fix 2 all nve"
		input_lines.append(fixnve)
		input_lines.append("timestep 0.0")
		input_lines.append("dump_modify couldump every 60")
		input_lines.append("reset_timestep " + str(ts_start + (i+1)*dumpevery -1))
		input_lines.append("run 1")
		input_lines.append("unfix 2")
		input_lines.append("dump_modify couldump every " + str(2 * timesteps))
		real_pot()
	input_lines.append("timestep " + str(timestep))
	input_lines.append("write_restart " + restart_file + "_after_coul")
	file = name + ".in"
	f = open(file, 'w')
	f.write("\n".join(input_lines))




def dump_to_potentials(dump_file,structure_file, timesteps, name = "", cores = 4,log = None, time = "00:01:00", queue = "regular", n = 1, C = "haswell", interactive = False): 
	input_file_1 = """clear
	units real
	dimension 3
	boundary p p p
	atom_style full
	atom_modify map array
	# define potential
	bond_style harmonic
	angle_style harmonic
	dihedral_style charmm
	kspace_style pppm 1.0e-5
	kspace_modify diff ad
	pair_style hybrid buck/coul/long 10 10 lj/cut/coul/long 10 10 lj/charmm/coul/long 10 10.1 10
	read_data """
	# add in structure file 

	input_file_2 = """pair_coeff 5 5 buck/coul/long 0.0 0.131258 0.0
	pair_coeff 5 6 buck/coul/long 0.0 0.321737 0.0
	pair_coeff 6 6 buck/coul/long 0.0 0.482217 0.0
	pair_coeff 3 5 buck/coul/long 0.0 0.150947 0.0
	pair_coeff 4 5 buck/coul/long 0.0 0.150947 0.0
	pair_coeff 3 6 buck/coul/long 0.0 0.342426 0.0
	pair_coeff 4 6 buck/coul/long 0.0 0.342426 0.0
	pair_coeff 1 5 lj/cut/coul/long 0.0 2.26454
	pair_coeff 2 5 lj/cut/coul/long 0.0 2.70999
	pair_coeff 1 6 lj/cut/coul/long 0.0 2.75
	pair_coeff 2 6 lj/cut/coul/long 0.0 3.1
	pair_coeff 3 3 lj/cut/coul/long 0.0 3.25
	pair_coeff 1 3 lj/cut/coul/long 0.0 2.1595
	pair_coeff 3 4 lj/cut/coul/long 0.0 3.3248
	pair_coeff 2 3 lj/cut/coul/long 0.0 2.605
	pair_coeff 1 1 lj/cut/coul/long 0.0 1.0691
	pair_coeff 1 4 lj/cut/coul/long 0.0 2.2344
	pair_coeff 1 2 lj/charmm/coul/long 0.0 1.5145
	pair_coeff 4 4 lj/cut/coul/long 0.0 3.3997
	pair_coeff 2 4 lj/cut/coul/long 0.0 2.6798
	pair_coeff 2 2 lj/cut/coul/long 0.0 1.96
	
	neighbor 3.0 bin
	neigh_modify delay 0 every 1 check yes page 500000 one 50000

	thermo 100
	thermo_style custom dt step pe temp lx ly lz press # pxx pyy pzz #c_eatoms

	compute Eksp all pe/atom kspace
	compute Ecoul all pe/atom pair
	group PbI type 5 6
	dump position_dump PbI custom 60 """ + name + """_coul.dump id x y z c_Ecoul
	timestep 0.0"""

	timestep_run_input ="""fix 1 all nve
	run 1
	unfix 1"""
	
	input_file = input_file_1 + structure_file + "\n" + input_file_2

	for step in timesteps: 
		input_file += "\n" + "read_dump " + dump_file + " " + str(step) + " x y z box yes" + "\n" + timestep_run_input

	input_name = "TEMP_" + name + "_coul" + ".in"
	f = open(input_name, 'w')
	f.write(input_file)
	#timestep_file.append("TEMP_coul_ts" + str(step) + "_" + name + ".dump")
	#load = subprocess.Popen("module load lammps/2017.03.31-hsw")
	command = ["srun", "-n", str(n), "-c", str(cores), "-C", C, "-q", queue, "-t",time , "lmp_cori", "-in", input_name] 
	if log: 
		command += [">>", log]
	call = subprocess.Popen(command)

def merge_timestep_files(timesteps, filename): 
	timestep_files = []
	f = open(filename,"w")
	for step in timesteps: 
		step_file = "TEMP_coul_ts" + str(step) + "_" + name + ".dump"
		command = ["cat", stepfile, ">>", filename]
		call = subprocess.Popen(command)