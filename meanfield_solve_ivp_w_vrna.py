from scipy.integrate import solve_ivp
import numpy
import scipy
import sys

def solver(pdic,t_samples,lab,rV):
	
	# function to find initial amount of infectious virus at t=-1h based on MOI (single-cycle assay)
	def find_V0(MOI):
		def func(V0):
			def ode(t,y):
				(T,V) = (y[0],y[-1])
				bTV = pdic['beta']*T*V/pdic['S']
				dT = -pdic['gamma']*bTV
				dV = -pdic['c']*V-bTV
				return numpy.hstack((dT,dV))
			sol = solve_ivp(ode,[0,1],[pdic['Nx'],V0[0]],method='BDF',t_eval=[1])
			return sol.y[0][-1]/pdic['Nx']-numpy.exp(-MOI)
		x0 = (pdic['c']+pdic['beta']*pdic['Nx']/pdic['S'])*MOI \
			/(pdic['gamma']*pdic['beta']/pdic['S']*(1-numpy.exp(-(pdic['c']+pdic['beta']*pdic['Nx']/pdic['S'])*1.0)))
		x = scipy.optimize.fsolve(func=func,x0=x0)
		if numpy.abs(func(x)) > 1e-4:
			return None
		return x[0]
	
	# MFM
	def ODE(t,y):
		(T,E,I,V,Vrna) = (y[0],y[1:1+pdic['nE']],y[-2-pdic['nI']:-2],y[-2],y[-1])
		# Pre-calculation of shared terms
		sI = numpy.sum(I)
		bTV = pdic['beta']*T*V/pdic['S']
		kE = (pdic['nE']/pdic['tauE'])*E
		dI = (pdic['nI']/pdic['tauI'])*I
		# ODEs
		dT = -pdic['gamma']*bTV
		dE1 = pdic['gamma']*bTV-kE[0]
		dEi = -numpy.diff(kE)
		dI1 = kE[-1]-dI[0]
		dIi = -numpy.diff(dI)
		dV = pdic['p']*sI-pdic['c']*V-bTV
		dVrna = pdic['prna']*sI-pdic['crna']*Vrna-bTV
		return numpy.hstack((dT,dE1,dEi,dI1,dIi,dV,dVrna))
	
	# define some empty lists
	t = []; V = []; Vrna = []; V_samples = []; Vrna_samples = []
	
	if lab == 'SC': # single-cycle assay
		# pre-rinse
		y0 = numpy.zeros(pdic['nE']+pdic['nI']+3)
		y0[0] = pdic['Nx']
		y0[-2] = find_V0(pdic['MOIsc'])
		if y0[-2] is None or numpy.isnan(y0[-2]):
			return 5*[None]
		y0[-1] = pdic['V0scrna']*pdic['S']
		if (y0[-2] > y0[-1]): 
			return 5*[None]
		sol = solve_ivp(ODE,[0,1],y0,method='BDF',t_eval=[1])
		y0 = sol.y[:,-1]
		# post-rinse
		y0[-2] = pdic['V0pr']*pdic['S']/(1-rV)
		y0[-1] = pdic['V0prrna']*pdic['S']
		if (y0[-2] > y0[-1]): 
			return 5*[None]
		for i in range(len(t_samples)-1):
			t.append(numpy.linspace(t_samples[i],t_samples[i+1],round(300/len(t_samples))))
			sol = solve_ivp(ODE,[t_samples[i],t_samples[i+1]],y0,method='BDF',t_eval=t[i])
			V.append(sol.y[-2]); Vrna.append(sol.y[-1])
			V_samples.append(V[i][-1]); Vrna_samples.append(Vrna[i][-1])
			y0 = sol.y[:,-1]
			y0[-2] *= pdic['frinse']
			y0[-1] *= pdic['frinse']
	
	else: # multiple-cylce assay
		y0 = numpy.zeros(pdic['nE']+pdic['nI']+3)
		y0[0] = pdic['Nx']
		y0[-2] = pdic['V0mc']*pdic['S']/(1-rV)
		y0[-1] = pdic['V0mcrna']*pdic['S']
		if (y0[-2] > y0[-1]): 
			return 5*[None]	
		V_samples.append(y0[-2]); Vrna_samples.append(y0[-1])
		for i in range(len(t_samples)-1):
			t.append(numpy.linspace(t_samples[i],t_samples[i+1],round(300/len(t_samples))))
			sol = solve_ivp(ODE,[t_samples[i],t_samples[i+1]],y0,method='BDF',t_eval=t[i])
			V.append(sol.y[-2]); Vrna.append(sol.y[-1])
			V_samples.append(V[i][-1]); Vrna_samples.append(Vrna[i][-1])
			y0 = sol.y[:,-1]
			y0[-2] *= pdic['frinse']
			y0[-1] *= pdic['frinse']
	
	# return results
	[t, V, Vrna, Vsamples, Vrna_samples] = [
		numpy.concatenate(t),
		numpy.concatenate(V),
		numpy.concatenate(Vrna),
		numpy.array(V_samples),
		numpy.array(Vrna_samples)
	]
	return t, V, Vrna, Vsamples, Vrna_samples
