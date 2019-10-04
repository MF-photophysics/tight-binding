import numpy as np
import matplotlib.pyplot as plt 
import brokenaxes as bax

#plt.rc("font", **{"family" : "normal", "weight" : "bold", "size" : 12})
n8_boxes = {100: [49.4986,50.4537,49.47606], 120:[49.55618, 50.515568,49.533306], 140:[49.53516, 50.600488, 49.511046], 150: [49.592049, 50.6261646, 49.5719502], 160: [49.669255, 50.767169, 49.669255],180:[49.7314733, 50.7685564, 49.7314744], 200: [49.790474, 50.77616, 49.7904746], 220: [49.865583, 50.771677, 49.865583], 260: [49.9888239, 50.74272048, 48.9888239], 300:[50.1278119, 50.626237, 50.1278119]}

def plot_gaps(allgaps, timestep = 30, display = True, color = "b", title = "Bandgap Over Time"): 
	gaps = [float(gap) for gap in open(allgaps,'r').read().split("\n")[:-1]]
	plt.plot([timestep * i for i in range(len(gaps))], gaps, color)
	plt.title(title)
	plt.ylabel("Bandgap (eV)")
	plt.xlabel("time (fs)")
	if display: 
		plt.show()

def read_bands(bands_file): 
	bands = {}
	lines = [[round(float(i), 15) for i in  line.split(" ")] for line in open(bands_file, "r").read().split("\n")[:-1]]
	for line in lines: 
		kpt = tuple(line[:3])
		try:
			bands[kpt].append(line[3])
		except KeyError: 
			bands[kpt] = [line[3]]
	return bands

def axsq_fit(points):
	A = np.array([[p[0]**2 for p in points]]).T
	b = np.array([[p[1] for p in points]]).T
	return float(np.linalg.inv(A.T @ A) @ A.T @ b)

def conduction_band_min(bands,E_fermi = 4.7,kpoint = (0,0,0)):
	kenergies = bands[kpoint]
	for n in range(len(kenergies)):
		if kenergies[n] > E_fermi and kenergies[n-1] < E_fermi: 
			return kenergies[n]

def valence_band_max(bands,E_fermi = 4.7,kpoint = (0,0,0)):
	kenergies = bands[kpoint]
	for n in range(len(kenergies)):
		if kenergies[n] > E_fermi and kenergies[n-1] < E_fermi: 
			return kenergies[n-1]

def electron_mass(bands,E_fermi = 4.7,box = np.array([50,50,50])):
	gamma = conduction_band_min(bands,E_fermi,(0,0,0))
	conv  = 7.61996563 # (hbar^2 / angstrom^2 * eV * electron mass )
	kx0,ky0,kz0 = 2*np.pi/box[0],2*np.pi/box[1], 2*np.pi / box[2]
	mx = conv * kx0**2/ (2*(axsq_fit([(0.1,conduction_band_min(bands,E_fermi, (0.1,0,0)) - gamma),(0.2,conduction_band_min(bands,E_fermi, (0.2,0,0)) - gamma)])))
	my = conv * ky0**2/ (2*(axsq_fit([(0.1,conduction_band_min(bands,E_fermi, (0,0.1,0)) - gamma),(0.2,conduction_band_min(bands,E_fermi, (0,0.2,0)) - gamma)])))
	mz = conv * kz0**2/ (2*(axsq_fit([(0.1,conduction_band_min(bands,E_fermi, (0,0,0.1)) - gamma),(0.2,conduction_band_min(bands,E_fermi, (0,0,0.2)) - gamma)])))
	return mx,my,mz

def hole_mass(bands,E_fermi = 4.7,box = np.array([50,50,50])):
	gamma = valence_band_max(bands,E_fermi,(0,0,0))
	conv  = 7.61996563 # (hbar^2 / angstrom^2 * eV * electron mass )
	kx0,ky0,kz0 = 2*np.pi/box[0],2*np.pi/box[1], 2*np.pi / box[2]
	mx = -conv*kx0**2/ (2*(axsq_fit([(0.1,valence_band_max(bands,E_fermi, (0.1,0,0))-gamma),(0.2,valence_band_max(bands,E_fermi, (0.2,0,0)) - gamma)])))
	my = -conv*ky0**2/ (2*(axsq_fit([(0.1,valence_band_max(bands,E_fermi, (0,0.1,0))-gamma),(0.2,valence_band_max(bands,E_fermi, (0,0.2,0)) - gamma)])))
	mz = -conv*kz0**2/ (2*(axsq_fit([(0.1,valence_band_max(bands,E_fermi, (0,0,0.1))-gamma),(0.2,valence_band_max(bands,E_fermi, (0,0,0.2)) - gamma)])))
	return mx,my,mz


def plot_bands(bands, kpt_to_index = lambda kpt: np.array(kpt) @ np.array([-1, 0,1]), marker = "bo" , title = "Band Structure" , display = True, output_structure = False, E_fermi = 4.7):
	kpt_indices = []
	energies = []
	energies_per_kpt = len(bands[(0,0,0)])
	kpts = len(bands.keys())
	for i in range(energies_per_kpt): # to get states per kpt, assuming gamma point is included
		for kpt in bands.keys(): 
			kpt_indices.append(kpt_to_index(kpt))
			energies.append(bands[kpt][i])
	i = 0
	#while i < len(kpt_indices): 
	#	xvals,yvals = [],[]	
	plt.plot(kpt_indices, energies,marker)
	plt.title(title)
	plt.ylabel("Energy (eV)")
	plt.xlabel("k point index")
	if display: 
		plt.show()
	if output_structure: 
		return kpt_indices, energies

def plot_function(files, func, timestep = 30, title = "",colors = ["b","g","r", "k"],ylabel = "", axes = "auto",display = True,valsonly = False,xlabel = "time (fs)"): 
	bands = [read_bands(file) for file in files]
	tvals = [timestep * i for i in range(len(files))]
	fvals = [func(b) for b in bands]
	if valsonly: 
		return fvals
	
	if axes == "auto":
		axes = plt
	try: 
		mean,sd = [],[]
		for i in range(len(fvals[0])):
			axes.plot(tvals,[val[i] for val in fvals], colors[i])
			mean.append(sum([x[i] for x in fvals]) / len(fvals))
			sd.append((sum([(x[i]-mean[i])**2 for x in fvals]) / len(fvals))**(1/2))
	except TypeError:
		axes.plot(tvals, fvals,colors[0])
		mean = sum([x for x in fvals]) / len(fvals)
		sd = (sum([(x-mean)**2 for x in fvals]) / len(fvals))**(1/2)
	if axes == plt:
		plt.title(title, fontsize = 22)
		plt.xlabel(xlabel,fontsize = 18)
		plt.ylabel(ylabel,fontsize = 18)
	else:
		axes.set_title(title, fontsize = 22)
		axes.set_xlabel(xlabel, fontsize = 18)
		axes.set_ylabel(ylabel,fontsize = 18)
	if display: 
		#plt.tight_layout()
		plt.show()
	return mean,sd

def band_edge_plot(files, margin = 0.03, title = "VBM and CBM"):
	mean = lambda x : sum(x) / len(x)
	vbm = mean(plot_function(files, valence_band_max, valsonly = True))
	cbm = mean(plot_function(files, conduction_band_min, valsonly = True))
	baxes = bax.brokenaxes(ylims=((vbm - margin, vbm + margin),(cbm - margin, cbm + margin)), hspace = 0.05)
	plot_function(files, valence_band_max, title = title, timestep = 0.03 ,xlabel = "time (ps)",ylabel = "Energy (eV)",display = False,axes = baxes,colors = ["b"])
	plot_function(files, conduction_band_min, title = title, timestep = 0.03 ,xlabel = "time (ps)",ylabel = "Energy (eV)",display = True,axes = baxes,colors = ["g"])
	
temps = {t:["T" +str(t) + "work/bands" + str(n) + ".dat" for n in range(1,1001)] for t in [100,120,140,150,160,180,200,220,260,300]}

def mass_temp_plot(temp_files, directions = ["x","y","z"],func = hole_mass,title = "mass", display = True):
	mean = lambda x : sum(x) / len(x)
	sd = lambda x, mean : (sum([(i-mean)**2 for i in x])/ len(x))**(1/2)
	xvals, yvals, zvals = {},{},{}
	temps = list(temp_files.keys())
	xaves,yaves,zaves = [],[],[]
	xsds,ysds,zsds = [],[],[]
	for t in temp_files:
		vals = plot_function(temp_files[t], func, valsonly = True)
		xvals[t], yvals[t], zvals[t] = [v[0] for v in vals], [v[1] for v in vals], [v[2] for v in vals] 
	
	if "x" in directions:
		xaves = {t:mean(xvals[t]) for t in temps}
		xsds = {t:sd(xvals[t], xaves[t]) for t in temps}
		plt.plot(list(xaves.keys()), list(xaves.values()), 'b')
	if "y" in directions: 
		yaves = {t:mean(yvals[t]) for t in temps}
		ysds = {t:sd(yvals[t], yaves[t]) for t in temps}
		plt.plot(list(yaves.keys()), list(yaves.values()), 'g')
	if "z" in directions: 
		zaves = {t:mean(zvals[t]) for t in temps}
		zsds = {t:sd(zvals[t], zaves[t]) for t in temps}
		plt.plot(list(zaves.keys()), list(zaves.values()), 'r')
	plt.xlabel("temperature (K)")
	plt.ylabel("hole mass" if func is hole_mass else "electron mass")
	plt.title(title)
	if display: 
		plt.show()

def write_data(files, outfile, functions = [lambda bands: conduction_band_min(bands) - valence_band_max(bands), valence_band_max, conduction_band_min, hole_mass, electron_mass]):
	data = []
	for f in functions:
		output = plot_function(files, f,valsonly = True)
		try: 
			for i in range(len(output[0])):
				dat = [out[i] for out in output]
				data.append(dat)
		except TypeError:
			data.append(output)
	outlines = []
	for i in range(len(files)):
		line = " ".join([str(data[c][i]) for c in range(len(data))])
		outlines.append(line)
	f = open(outfile, "w")
	f.write("\n".join(outlines))
