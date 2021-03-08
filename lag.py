
import numpy as np
import matplotlib.pyplot as plt
from Data_analysis.Baseline import TD_bs,WhittakerSmooth
from scipy.interpolate import interp1d
import os



def get_lag(data,band,bins,wind=None,sigma = 4,mcmc_num = 500,plot_savedir = None,Positive_lag_positive_value = True):
	
	binsize = bins[1]-bins[0]
	t_c = 0.5*(bins[1:]+bins[:-1])
	if wind is not None:
		wind_index = np.where((t_c>=wind[0])&(t_c<=wind[-1]))[0]
	else:
		wind_index = np.arange(0,len(t_c),dtype = int)
	t,energy = data
	event_band = get_event_with_band(t,energy,band)
	cs0 = None
	cs_err0 = None
	lag_errh2 = 0
	lag_errl2 = 0
	lag_all = 0
	return_list = []
	lc_list = []
	for index_,(t_i,e_i) in enumerate(event_band):
		
		if cs0 is None or cs_err0 is None:
			num_i = np.histogram(t_i,bins=bins)[0]
			num_err_i = np.sqrt(num_i)
			rate_i = num_i/binsize
			rate_err_i = num_err_i/binsize
			csi,bsi,sigma_i = TD_bs(t_c,rate_i,sigma=True)
			lc_list.append([rate_i,bsi,sigma_i])
			SNRi = csi/sigma_i
			if SNRi.max() > sigma:
				cs0 = csi
				cs_err0 = rate_err_i
				return_list.append([index_,0,1,1])
		else:
			num_i = np.histogram(t_i,bins=bins)[0]
			num_err_i = np.sqrt(num_i)
			rate_i = num_i/binsize
			rate_err_i = num_err_i/binsize
			csi,bsi,sigma_i = TD_bs(t_c,rate_i,sigma=True)
			lc_list.append([rate_i,bsi,sigma_i])
			SNRi = csi/sigma_i
			if SNRi.max() > sigma:
				if plot_savedir is not None:
					if os.path.exists(plot_savedir) == False:
						os.makedirs(plot_savedir)
					savename = plot_savedir + 'A_nccf_' + str(index_) + '.png'
				else:
					savename = None
				lag,lag_errl,lag_errh = get_one_lag(cs0[wind_index],csi[wind_index],cs_err0[wind_index],rate_err_i[wind_index],t_c[wind_index],mcmc_num=mcmc_num,save = savename)
				lag_all = lag_all + lag
				lag_errl2 = lag_errl2+lag_errl**2
				lag_errh2 = lag_errh2 + lag_errh**2
				return_list.append([index_,lag_all,np.sqrt(lag_errl2),np.sqrt(lag_errh2)])
				cs0 = csi
				cs_err0 = rate_err_i
	lag0 = np.array(return_list).T
	if Positive_lag_positive_value:
		lag = lag0
	else:
		lag = np.vstack([lag0[0], -1 * lag0[1], lag0[3], lag0[2]])
		
	return {
		'band':band,
		'wind':wind,
		'lag':lag,
		'lc':{
			'time':t_c,
			'rate':lc_list
		}
	}
	

def get_one_lag(cs1,cs2,cs1_err,cs2_err,time,mcmc_num = 1000,save = None):
	
	n = len(time)
	dt = np.abs(time[1]-time[0])
	lag_t = np.arange(-n+1,n,1)*dt
	nccf = np.correlate(cs1 / cs1.max(), cs2 / cs2.max(), 'full')
	lag_time = get_nccf_peak_time(lag_t,nccf,save = save)
	lag_list = []
	for i in range(mcmc_num):
		rand1 = np.random.randn(n)
		cs1_sim = cs1 + cs1_err*rand1
		rand2 = np.random.randn(n)
		cs2_sim = cs2 + cs2_err*rand2
		nccf_sim = np.correlate(cs1_sim/cs1_sim.max(),cs2_sim/cs2_sim.max(),'full')
		lag_sim = get_nccf_peak_time(lag_t,nccf_sim)
		lag_list.append(lag_sim)
	lag_list = np.array(lag_list)
	lag_mean = lag_list.mean()
	lag_sigma = lag_list.std(ddof = 1)
	dx = lag_mean-lag_time
	lag_errl = lag_sigma-dx
	lag_errh = lag_sigma+dx
	
	return lag_time,lag_errl,lag_errh
	
	
def get_nccf_peak_time(lag_t,nccf,precision = 0.001,save = None):
	
	dt = np.abs(lag_t[1]-lag_t[0])
	w = np.ones(len(nccf))
	smoo_nccf = WhittakerSmooth(nccf,w,0.15/dt**1.5)
	intef = interp1d(lag_t,smoo_nccf,kind = 'quadratic')
	new_lat_t = np.arange(lag_t[1],lag_t[-2]+precision,precision)
	new_nccf = intef(new_lat_t)
	index_ = np.argmax(new_nccf)
	t_max = new_lat_t[index_]
	if save is not None:
		fig = plt.figure(constrained_layout=True)
		ax = fig.add_subplot(1,1,1)
		ax.plot(lag_t,nccf,'.',color = 'k',label = 'nccf')
		ax.plot(new_lat_t,new_nccf,'-',color='#f47920',label = 'fit')
		ax.axvline(x =t_max,color = 'r' )
		ax.set_xlim(t_max-5,t_max+5)
		ax.set_xlabel('lag time (s)')
		ax.set_ylabel('nccf')
		ax.legend()
		fig.savefig(save)
		plt.close(fig)
	return t_max

def get_event_with_band(time,energy,band):
	
	return_list = []
	
	for (el,eh) in band:
		index = np.where((energy>=el)&(energy<=eh))[0]
		t_i = time[index]
		e_i = energy[index]
		return_list.append([t_i,e_i])
	
	return return_list

