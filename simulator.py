import meanfield_solve_ivp_w_vrna
import numpy
import scipy
import sys

def load():
	# load base parameter dictionary
	pdic = dict(
		p = 10**0.11*1.781*2,  # infectious virus production rate
		beta = 10**-5.3/1.781, # viral entry rate 
		gamma = 0.5, # probability of a successful cell infection post viral entry
		prna = 10**3.03, # viral RNA production rate
		tauE = 7, # average length of eclipse phase
		tauI = 41, # average length of infectious phase
		nE = 60, # number of eclipse phase compartments
		nI = 60, # number of infectious phase compartments
		c = 0.0573, # rate of virus loss of infectivity
		crna = 0.001, # viral RNA degradation rate
		Nx = 1.9e6, # number of cells in infection experiment
		S = 10., # volume of supernatant in infection experiments
		Nxassay = 1e5, # number of cells in ED assay well
		Vinoc = 0.05, # volume of inoculum placed in ED assay well
		Vvir = 5.236e-16, # volume of one virion
		MOIsc = 3., # multiplicity of infection for single-cycle assay
		V0pr = 10**6.00, # SIN/mL concentration post-rinse (single-cycle)
		V0mc = 10**1.79, # SIN/mL concentration initially (multiple-cycle)
		V0scrna = 10**9.97, # RNA/mL cocentration initially (single-cycle)
		V0prrna = 10**6.76, # RNA/mL concentration post-rinse (single-cycle)
		V0mcrna = 10**4.02, # RNA/mL concentration initially (multiple-cycle)
		sigma_MC_RNA = 0.258, # characterizes variability in multiple-cycle total virus data
		sigma_SC_RNA = 0.237, # characterizes variability in single-cycle total virus data 
		frinse = 1-0.5/10, # rinse factor
		lim = 1.4940 # log10(SIN/mL) ED assay limit of detection 
	)
	
	# load sampling times and dilutions
	base = 'psimon-sH1N1/sH1N1_'
	dat_MC_ED = numpy.loadtxt(base+'MC_ED.dat').T
	dat_SC_ED = numpy.loadtxt(base+'SC_ED.dat').T
	t_samples = dict(
		MC = dat_MC_ED[0][0::3],
		SC = numpy.concatenate(([0],dat_SC_ED[0][0::3]))
	)
	D = dict(
		MC = dat_MC_ED[1:9].T[0::3],
		SC = dat_SC_ED[1:9].T[0::3]
	)
	
	return pdic, t_samples, D

def get_rV(pdic):
	cogr = (1+pdic['c']/(pdic['beta']*pdic['Nxassay']/pdic['Vinoc']))/pdic['gamma']
	B = pdic['p']*pdic['tauI']; nI = pdic['nI']
	f = lambda rV: 1-cogr*(1-rV)-(B*(1-rV)/nI+1)**-nI
	x = scipy.optimize.fsolve(func=f,x0=1-1/cogr)[0]
	if numpy.abs(f(x)) > 1e-5:
		return None	
	return x

def run(pdic,t_samples,D):
	
	# check parameters
	if pdic['prna'] < 1e-3 or \
		pdic['prna'] > 1e12 or \
		pdic['p'] < 1e-3 or \
		pdic['p'] > pdic['prna'] or \
		pdic['gamma'] < 1e-6 or \
		pdic['gamma'] > 1 or \
		pdic['beta'] < 1e-10 or \
		pdic['beta'] > 1e-2 or \
		pdic['tauE'] < 0 or \
		pdic['tauE'] > 720 or \
		pdic['tauI'] < 0 or \
		pdic['tauI'] > 720:
		return 4*[None]	
	
	# obtain extinction probability 
	rV = get_rV(pdic)
	if rV is None:
		return 4*[None]	
	
	# check at t=0h, # infectious virions <= # viral RNA copies
	if rV > 1-min(pdic['V0mc']/pdic['V0mcrna'],pdic['V0pr']/pdic['V0prrna']):
		return 4*[None]

	dat = []	
	for lab in ['MC','SC']:
		# run infection experiment
		_,_,_,V,Vrna = meanfield_solve_ivp_w_vrna.solver(pdic,t_samples[lab],lab,rV)
		if V is None or min(V) < 0 or min(Vrna) < 0:
			return 4*[None]
	
		# simulate ED assay experiments
		sim_dat_ED = []
		for v,d,t in zip(V,D[lab],t_samples[lab][lab=='SC':]):
			V0 = numpy.random.binomial(
				n = (pdic['Vinoc']*10**d/pdic['Vvir']).astype(int),
				p = v/pdic['S']*pdic['Vvir'],
				size = [4,8]
			)
			ninf = numpy.sum(numpy.random.random((4,8)) >= rV**V0,axis=0)
			sim_dat_ED.append(numpy.concatenate(([t],d,ninf)))
		sim_dat_ED = numpy.array(sim_dat_ED)
		dat.append(sim_dat_ED)
	
		# adding RNA measurement noise
		Crna = Vrna/pdic['S']
		Crna *= 10.**(numpy.random.normal(0,pdic['sigma_%s_RNA'%lab],len(Vrna)))
		
		# formatting simulated RNA data
		expon = numpy.floor(numpy.log10(Crna))
		coeff = numpy.around(Crna/10**expon,2)
		Vrna = coeff*10**expon
		sim_dat_RNA = numpy.concatenate((t_samples[lab][lab=='SC':],Crna)).reshape(2,-1).T
		dat.append(sim_dat_RNA)
	
	# return results
	sim_dat_MC_ED, sim_dat_MC_RNA, sim_dat_SC_ED, sim_dat_SC_RNA = dat
	return sim_dat_MC_ED, sim_dat_MC_RNA, sim_dat_SC_ED, sim_dat_SC_RNA
