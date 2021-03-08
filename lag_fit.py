import pymultinest
import numpy as np
import pandas as pd
from ..spectrum.analysis import Analyzer
from ..spectrum.PlotMarginalModes import PlotMarginalModes
import matplotlib.pyplot as plt

class Lag_fit(object):
	
	def __init__(self,model,priorlist,parameters,result = None, lag_data = None):
		
		if result is not None:
			band = result['band']
			index_,lag,lag_errl,lag_errh = result['lag']
			index_ = index_.astype(np.int)
			band = band[index_]
			self.energyc =  np.sqrt(band[:,0]*band[:,1])
			self.lag = lag
			self.lag_errl = lag_errl
			self.lag_errh = lag_errh
		elif lag_data is not None :
			if isinstance(lag_data,pd.DataFrame):
				band_l = lag_data['band_l'].values
				band_h = lag_data['band_h'].values
				self.lag = lag_data['lag'].values
				self.lag_errl = lag_data['lag_errl'].values
				self.lag_errh = lag_data['lag_errh'].values
				self.energyc = np.sqrt(band_l*band_h)
			else:
				self.energyc = np.array(lag_data[0])
				self.lag = np.array(lag_data[1])
				self.lag_errl = np.array(lag_data[2])
				self.lag_errh = np.array(lag_data[3])
			
		self.model = model
		self.prior_list = priorlist
		self.parameters = parameters
		self.n_params = len(parameters)
		

	def prior(self,cube,ndim,nparams):
		for i in range(nparams):
			cube[i] = self.prior_list[i].get_value(cube[i])

	def log_like(self,cube,ndim,nparams):
		
		model_dt = self.model(self.energyc,cube)
		diff = self.lag - model_dt
		index = np.where(diff<=0)[0]
		err = self.lag_errl
		err[index] = self.lag_errh[index]
		loglik = -0.5*((diff/err)**2).sum()
		
		return loglik

	def run(self,outputfiles_basename,resume = False, verbose = True):
		
		pymultinest.run(self.log_like,self.prior,self.n_params,outputfiles_basename = outputfiles_basename,resume = resume, verbose = verbose)
		a1 = Analyzer(outputfiles_basename=outputfiles_basename, n_params = self.n_params)
		data = [self.energyc,self.lag,self.lag_errl,self.lag_errh]
		return Fit_plot(a1,self.parameters,data,self.model)

	
class Fit_plot(object):
	
	def __init__(self,a,parameters,data,model):
		self.result = a
		self.n_params = a.n_params
		self.parameters = parameters
		self.energyc = data[0]
		self.lag = data[1]
		self.lag_errl = data[2]
		self.lag_errh = data[3]
		self.model = model
	
	def get_analysis(self):
		return self.result
	def get_best_fit(self):
		c = self.result.get_stats()['modes'][0]
		maximum = np.array(c['maximum'])
		mean = np.array(c['mean'])
		sigma = np.array(c['sigma'])
		dx = mean - maximum
		sigmal = sigma - dx
		sigmah = sigma + dx
		return maximum,sigmal,sigmah
	
	def plot_corner(self):
		c = self.result.get_stats()['modes'][0]
		maximum = np.array(c['maximum'])
		mean = np.array(c['mean'])
		sigma = np.array(c['sigma'])
		dx = mean - maximum
		p = PlotMarginalModes(self.result)
		sigmal = sigma - dx
		sigmah = sigma + dx
		uplims = np.zeros(len(sigmal),dtype=bool)
		lolims = np.zeros(len(sigmal),dtype=bool)
		lolims[sigmal<=0]=True
		uplims[sigmah<=0]=True
		sigmal[sigmal < 0] = 0
		sigmah[sigmah < 0] = 0
		
		plt.figure(figsize=(5*self.n_params,5*self.n_params))
		for i in range(self.n_params):
			for j in range(i,self.n_params):
				if(i == 0):
					plt.subplot(self.n_params, self.n_params,i*self.n_params+j+1)
					plt.title(self.parameters[j],size = 30)
					plt.tick_params(labelsize = 15)
					p.plot_marginal(j, with_ellipses = True, with_points = False, grid_points=50)
					plt.errorbar(maximum[j],0,xerr=[[sigmal[j]],[sigmah[j]]],xlolims=lolims[j],xuplims=uplims[j],
					             markersize=10,elinewidth = 3,
					             fmt = 'o',ecolor = 'r',capthick=5,color = 'r')
					if(j == 0):
						plt.ylabel("Probability",size = 30)
						plt.xlabel(self.parameters[j],size = 30)
				else:
					plt.subplot(self.n_params, self.n_params,i*self.n_params+j+1)
					plt.tick_params(labelsize = 15)
					p.plot_conditional(j, i-1, with_ellipses = False, with_points = False, grid_points=30)
					plt.errorbar(maximum[j], maximum[i-1], xerr=[[sigmal[j]],[sigmah[j]]], yerr=[[sigmal[i-1]],[sigmah[i-1]]],
					             xlolims=lolims[j],xuplims=uplims[j],uplims=uplims[i-1],lolims=lolims[i-1],
					             markersize=10,elinewidth = 3,
					             fmt='o', ecolor='r', capthick=5, color='r')
					if(j == i):
						plt.xlabel(self.parameters[j],size = 30)
						plt.ylabel(self.parameters[i-1],size = 30)
		
	def plot_fit(self,ax = None):
		
		a1 = self.result
		for para in a1.get_equal_weighted_posterior()[::100,:-1]:
			if ax is not None:
				ax.plot(self.energyc,self.model(self.energyc,para),'-', color='k', alpha=0.3)
			else:
				plt.plot(self.energyc,self.model(self.energyc,para),'-', color='k', alpha=0.3)
		
		best_value = a1.get_best_fit()['parameters']
		if ax is not None:
			ax.plot(self.energyc,self.model(self.energyc,best_value),'-', color='r',label = 'best model')
		else:
			plt.plot(self.energyc,self.model(self.energyc,best_value),'-', color='r',label = 'best model')
		
	def plot_different(self,ax = None):
		
		a1 = self.result
		best_value = a1.get_best_fit()['parameters']
		model_lag = self.model(self.energyc,best_value)
		if ax is not None:
			ax.errorbar(self.energyc[1:],(self.lag-model_lag)[1:],yerr= [self.lag_errl[1:],self.lag_errh[1:]],elinewidth=1,capsize=2,fmt = '.',label = 'lag data')
			ax.axhline(y=0,label = 'model',linestyle = '-.')
		else:
			plt.errorbar(self.energyc[1:],(self.lag-model_lag)[1:],yerr= [self.lag_errl[1:],self.lag_errh[1:]],elinewidth=1,capsize=2,fmt = '.',label = 'lag data')
			plt.axhline(y=0,label = 'model',linestyle = '-.')
		
		