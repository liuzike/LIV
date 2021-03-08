
import numpy as np


def get_band(band,numb,ovelap=0.5,scale='log'):
	
	if scale == 'log':
		return get_band_in_log(band,numb,ovelap)
	else:
		return get_band_in_unif(band,numb,ovelap)


def get_band_in_log(band,numb,ovelap=0.5):
	if ovelap >= 0.95:
		ovelap = 0.95
	non_ovelap = 1-ovelap
	log_band = np.log10(band)
	len_ = np.abs(log_band[-1]-log_band[0])
	x = len_/(non_ovelap*(numb-1)+1)
	log_el_edges = np.linspace(log_band[0], log_band[-1] - x, numb)
	log_eh_edges = log_el_edges + x
	return np.vstack([10**log_el_edges,10**log_eh_edges]).T
	
	
def get_band_in_unif(band,numb,ovelap=0.5):
	
	if ovelap >= 0.95:
		ovelap = 0.95
	non_ovelap = 1-ovelap
	len_ = np.abs(band[-1]-band[0])
	x = len_/(non_ovelap*(numb-1)+1)
	el_edges = np.linspace(band[0],band[-1]-x,numb)
	eh_edges = el_edges+x
	
	return np.vstack([el_edges,eh_edges]).T






















