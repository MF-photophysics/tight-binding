from cell_geometry import *

def sighop(x,params):
	return params[0]*np.exp(-params[1]*x[0] - params[2]*x[1]**2) + params[3]

def pihop(x,params): 
	return params[0]*np.exp(-params[1]*x[0] - params[2]*x[1] - params[3]*x[2]) + params[4]

def jnchop(x,params):
	return params[0]*x[0]*np.exp(-params[1]*x[1] - params[2]*x[2])

def onsite(x,params,T): 
	# x[0]/params[0] are coulomb potential, x[1]/params[1] are hessian, params[2] is temperature coefficient, params[3] is 
	return params[0]*x[0] + params[1]*x[1] + params[2]*T + params[3]

# fitted parameters
spsig_20ts_150_300_params = [72.5498,1.29895,0.268934,-0.143073]
ppsig_20ts_150_300_params = [12.1361, 0.420229, 0.118432,-1.49759]
pppi_20ts_150_300_params = [23.5502, 1.17949, 0.394076, 0.127548, -0.129589]

pbs_20ts_150_300_params = [-1.59578, 0.00164292, -0.00114754, -5.18705]
pbp_20ts_150_300_params = [-0.692346, -0.38538, -0.00155891, 4.9896]
ip_20ts_150_300_params = [-1.48536, 0.0626952, -0.00117001, 4.75789]

# using volume (of 2x2x2, can change this)
pbs_20ts_150_220_300_params = [-1.56582, 0.0026258, -0.007925, 2.40268]
pbp_20ts_150_220_300_params = [-0.66472, -0.368484, -0.0108627, 15.3676]
ip_20ts_150_220_300_params = [-1.5067, 0.0617038, -0.00811448, 12.5022]

spsighop = lambda x, T: sighop(x,spsig_20ts_150_300_params)
ppsighop = lambda x, T: sighop(x,ppsig_20ts_150_300_params)
pppihop = lambda x, T: pihop(x,pppi_20ts_150_300_params)

pbs = lambda x, T: onsite(x,pbs_20ts_150_300_params, T)
pbp = lambda x, T: onsite(x,pbp_20ts_150_300_params, T)
ip = lambda x, T: onsite(x,ip_20ts_150_300_params, T)


psets1 = {"pbs": pbs_onsites, "pbp": pbp_onsites, "ip": ip_onsites, "spsig": spsig, "ppsig":ppsig, "pppi":pppi}
featsets1 = {"pbs": onsite_no_features, "pbp": onsite_no_features, "Ip" : onsite_no_features, "spsig":bond_min_features, "ppsig": bond_min_features, "pppi":bond_min_features}
model1 = {"pbs": pbs, "pbp": pbp, "ip": ip, "spsig" : spsighop, "ppsig": ppsighop, "pppi": pppihop}

def full_hopping(dump, timesteps, n_vec, T, filename, coul_dump = False,model = model1):
	features = write_param_file(dump, timesteps, n_vec, ortho_cell = False, N_ortho = None, param_sets = psets1, feature_sets = featsets1, coul = True,vdw = False, coul_dump = coul_dump,hessian = True, textfile = False)
	pbs = [model["pbs"](i[4:],T) for i in features["pbs"]]
	pbp = [model["pbp"](i[4:],T) for i in features["pbp"]]
	ip = [model["ip"](i[4:],T) for i in features["ip"]]
	spsig = [model["spsig"](i[4:],T) for i in features["spsig"]]
	ppsig = [model["ppsig"](i[4:],T) for i in features["ppsig"]]
	pppi = [model["pppi"](i[4:],T) for i in features["pppi"]]
	#for t in timesteps:
	#	for ix in range(n_vec[0]):
	#		for iy in range(n_vec)
	lines = []
	for i in range(len(timesteps) * n_vec[0]*n_vec[1]*n_vec[2]):
		lines.append(" ".join([str(i) for i in [pbs[i], pbp[3*i], pbp[3*i + 1], pbp[3*i + 2]]]))
		lines.append(" ".join([str(i) for i in [ip[9*i], ip[9*i], ip[9*i + 1], ip[9*i + 2]]]))
		lines.append(" ".join([str(i) for i in [ip[9*i + 3], ip[9*i + 3], ip[9*i + 4], ip[9*i + 5]]]))
		lines.append(" ".join([str(i) for i in [ip[9*i + 6], ip[9*i + 6], ip[9*i + 7], ip[9*i + 8]]]))
		lines.append(" ".join([str(i) for i in [spsig[6*i], ppsig[6*i], pppi[12*i], pppi[12*i + 1]]]))
		lines.append(" ".join([str(i) for i in [spsig[6*i + 1], ppsig[6*i + 1], pppi[12*i + 2], pppi[12*i + 3]]]))
		lines.append(" ".join([str(i) for i in [spsig[6*i + 2], ppsig[6*i + 2], pppi[12*i + 4], pppi[12*i + 5]]]))
		lines.append(" ".join([str(i) for i in [spsig[6*i + 3], ppsig[6*i + 3], pppi[12*i + 6], pppi[12*i + 7]]]))
		lines.append(" ".join([str(i) for i in [spsig[6*i + 4], ppsig[6*i + 4], pppi[12*i + 8], pppi[12*i + 9]]]))
		lines.append(" ".join([str(i) for i in [spsig[6*i + 5], ppsig[6*i + 5], pppi[12*i + 10], pppi[12*i + 11]]]))
	f = open(filename, "w")
	f.write("\n".join(lines))