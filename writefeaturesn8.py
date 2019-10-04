from cell_geometry import *
for t in [100,150,160,220,300,120,140,180,200,260]:
	d = analysis_tools.dump("T" + str(t) + "_exp8_run.dump", structure_tools.MAPbI3_cubic_n8, tilting = True, potentials = True)
	dc = analysis_tools.dump("T" + str(t) + "_exp8_coul.dump", structure_tools.MAPbI3_cubic_n8, tilting = True, potentials = True)
	write_param_file(d, [i for i in d.timesteps.keys() if i > 0], np.array([8,8,8]), param_sets = {"T" + str(t) + "_I-p_feat.dat":ip_onsites, "T" + str(t) + "_Pb-p_feat.dat":pbp_onsites, "T" + str(t) + "_Pb-s_feat.dat" : pbs_onsites}, feature_sets = {"T" + str(t) + "_I-p_feat.dat": onsite_no_features, "T" + str(t) + "_Pb-p_feat.dat": onsite_no_features, "T" + str(t) + "_Pb-s_feat.dat":onsite_no_features}, coul = True, vdw = True, coul_dump = dc, hessian = True)