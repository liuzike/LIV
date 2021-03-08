
import pandas as pd
import numpy as np


class Lag_save(object):
	
	def __init__(self,result):
		self.result = result
	
	def save_lag_csv(self,savename):
		
		index_,lag,lag_errl,lag_errh = self.result['lag']
		band_l,band_h = np.array(self.result['band'][index_.astype(np.int)]).T
		
		c = {'band_l':band_l,
		     'band_h':band_h,
		     'lag':lag,
		     'lag_errl':lag_errl,
		     'lag_errh':lag_errh}
		hl = pd.DataFrame(c)
		hl.to_csv(savename,index=False)
		





