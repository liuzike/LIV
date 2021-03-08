
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np


class Lag_plot(object):
	
	def __init__(self,results):
		
		self.results = results
		
	def plot_lightcurve(self,sigma = 1,fig = None):
		
		lc = self.results['lc']
		band = self.results['band']
		wind = self.results['wind']
		time = lc['time']
		rate_list = lc['rate']
		n = len(rate_list)
		
		if fig is None:
			fig = plt.figure(constrained_layout=True,figsize = (5,2*n))
		gs = GridSpec(n, 1, figure=fig)
		for index_,(ratei,bsi,sigmai) in enumerate(rate_list):
			bandi = band[index_]
			labeli = '%.1f-%.1f keV' % (bandi[0],bandi[-1])
			ax = fig.add_subplot(gs[-index_-1])
			ax.plot(time,ratei,'-',color = 'k',label = labeli)
			ax.plot(time,bsi,'-',color = 'r')
			ax.plot(time,bsi+sigma*sigmai,'-',color = 'g')
			if wind is not None:
				ax.set_xlim(wind[0],wind[-1])
			else:
				ax.set_xlim(time.min(),time.max())
			ax.set_ylabel('Rate')
			ax.legend()
			if index_ != 0 :
				ax.tick_params(labelbottom = False)
			else:
				ax.set_xlabel('Time')
		
		return fig
	
	def plot_lag(self,ax = None):
		band = self.results['band']
		index_,lag,lag_errl,lag_errh = self.results['lag']
		#index_ = np.trunc(index_)
		band_l, band_h = (band[index_.astype(np.int)]).T
		energyc = np.sqrt(band_l*band_h)
		if ax is None:
			fig = plt.figure(constrained_layout=True)
			ax = fig.add_subplot(1,1,1)
		ax.errorbar(energyc[1:],lag[1:],yerr = [lag_errl[1:],lag_errh[1:]],elinewidth=1,capsize=2,label = 'lag data',fmt = '.')
		if ax is None:
			return fig
		else:
			return ax







