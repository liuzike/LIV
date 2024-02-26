import os
import sys
import json
import time
import operator
import matplotlib
import numpy as np
matplotlib.use('Agg')
from tools import savetxt
from os.path import isdir
from tools import loadtxt
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import rcParams
from astropy.modeling import models
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip, mad_std


class CCF:

    def __init__(self, dt, x, y, x_err, y_err, save, save_dir, suffix):
        self.dt = dt
        self.x = np.array(x)
        self.y = np.array(y)
        self.x_err = np.array(x_err)
        self.y_err = np.array(y_err)
        self.save = save
        self.save_dir = save_dir
        self.suffix = suffix
        self.N = len(self.x)


    @staticmethod
    def gmo(x, amp, mean, std):
        g = models.Gaussian1D(amplitude=amp, mean=mean, stddev=std)
        return g(x)


    @staticmethod
    def ccf_band(d, x, y):
        x = np.array(x)
        y = np.array(y)
        N = len(x)
        sx = np.sum(x ** 2)
        sy = np.sum(y ** 2)
        nor = np.sqrt(sx * sy)
        low = max(1, 1 - d)
        upp = min(N, N - d)
        sxy = np.sum([x[i] * y[i+d] for i in range(low - 1, upp, 1)])
        ccf = sxy / nor
        return ccf


    @staticmethod
    def ccf_series(dt, x, y):
        x = np.array(x)
        y = np.array(y)
        N = len(x)
        tau = [dt * d for d in range(-N + 1, N, 1)]
        ccfs = [CCF.ccf_band(d, x, y) for d in range(-N + 1, N, 1)]
        tau, ccfs = tau[int(N / 4):int(7 * N / 4)], ccfs[int(N / 4):int(7 * N / 4)]
        return np.array(tau), np.array(ccfs)


    @staticmethod
    def correlate(dt, x, y):
        x = np.array(x)
        y = np.array(y)
        N = len(x)
        sx = np.sum(x ** 2)
        sy = np.sum(y ** 2)
        nor = np.sqrt(sx * sy)
        tau = [dt * d for d in range(-N + 1, N, 1)]
        ccfs = np.correlate(y, x, mode='full') / nor
        tau, ccfs = tau[int(N / 4):int(7 * N / 4)], ccfs[int(N / 4):int(7 * N / 4)]
        return np.array(tau), np.array(ccfs)


    def mc_simulation(self, nmc):
        self.nmc = nmc
        self.xs = [self.x]
        self.ys = [self.y]
        for m in range(self.nmc):
            self.xs = self.xs + [self.x + self.x_err * np.random.randn(self.N)]
            self.ys = self.ys + [self.y + self.y_err * np.random.randn(self.N)]
        if self.save:
            np.savetxt(self.save_dir + '/mc_xs%s.txt'%self.suffix, self.xs)
            np.savetxt(self.save_dir + '/mc_ys%s.txt'%self.suffix, self.ys)


    def get_lag(self):
        mc_ccfs, mc_lags = [], []
        for i, (x, y) in enumerate(zip(self.xs, self.ys)):
            tau, ccfs = CCF.correlate(self.dt, x, y)
            peak_idx = np.argmax(ccfs)
            lag = tau[peak_idx]
            mc_ccfs.append(ccfs)
            mc_lags.append(lag)

        tau, _ = CCF.correlate(self.dt, self.xs[0], self.ys[0])
        true_lag = mc_lags[0]

        mask = sigma_clip(mc_lags, sigma=5, maxiters=5, stdfunc=mad_std).mask
        not_mask = list(map(operator.not_, mask))
        mc_lags = np.array(mc_lags)[not_mask]

        low_lag, upp_lag = np.percentile(mc_lags, [16, 84])
        true_lag_err = np.diff([low_lag, true_lag, upp_lag])
        lagpm = [true_lag, true_lag_err[0], true_lag_err[1]]

        if self.save:
            np.savetxt(self.save_dir + '/delay_times%s.txt'%self.suffix, np.array(tau)+0)
            np.savetxt(self.save_dir + '/mc_ccfs%s.txt'%self.suffix, np.array(mc_ccfs)+0)
            np.savetxt(self.save_dir + '/mc_lags%s.txt'%self.suffix, np.array(mc_lags)+0)
            with open(self.save_dir + '/lagpm%s.json'%self.suffix, 'w') as f_obj:
                json.dump(lagpm, f_obj, indent=4)
            print('\n-------------------------------------------------')
            print('lag calculator:')
            with open(self.save_dir + '/lagpm%s.json'%self.suffix) as f_obj:
                for line in f_obj:
                    print(line.rstrip())
            print('-------------------------------------------------')

            rcParams['mathtext.fontset'] = 'stix'
            rcParams['font.family'] = 'Times New Roman'
            rcParams['font.serif'] = 'Times New Roman'
            rcParams['font.size'] = '16'

            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.scatter(tau, mc_ccfs[0], marker='+', color='black', s=20, linewidths=0.5, alpha=1.0)
            for i in np.arange(1, len(mc_ccfs), 100):
                ax.scatter(tau, mc_ccfs[i], marker='+', color='blue', s=10, linewidths=0.5, alpha=0.5)
            ax.set_xlabel('Time Delay (sec)')
            ax.set_ylabel('CCF Value')
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/tau_ccfs%s.pdf'%self.suffix, bbox_inches='tight', pad_inches=0.1, dpi=100)
            print('\n+++++tau_ccfs%s.pdf+++++'%self.suffix)
            plt.close(fig)

            bins = np.linspace(min(mc_lags), max(mc_lags), 30)
            histvalue, histbin = np.histogram(mc_lags, bins=bins)
            lags = (histbin[1:] + histbin[:-1]) / 2
            counts = histvalue
            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.plot(lags, counts, ds='steps-mid', color='blue')
            ax.axvline(true_lag, c='grey')
            ax.axvline(low_lag, c='grey', ls='--')
            ax.axvline(upp_lag, c='grey', ls='--')
            ax.set_xlabel('lags (sec)')
            ax.set_ylabel('counts')
            ax.set_title(r'$lag=%.2f_{-%.2f}^{+%.2f}\,sec$'%(true_lag, true_lag_err[0], true_lag_err[1]))
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/lags_pdf%s.pdf'%self.suffix, bbox_inches='tight', pad_inches=0.1, dpi=100)
            print('\n+++++lags_pdf%s.pdf+++++'%self.suffix)
            plt.close(fig)

        return lagpm


    def fit_tau_ccf_gau(self, nf):
        self.nf = nf
        tau0, ccfs0 = CCF.correlate(self.dt, self.xs[0], self.ys[0])
        mean0 = np.average(tau0, weights=ccfs0)
        std0 = np.sqrt(np.average((tau0 - mean0) ** 2, weights=ccfs0))
        g0 = models.Gaussian1D(amplitude=1, mean=mean0, stddev=std0)
        amp0 = np.max(ccfs0) / g0(mean0)

        try:
            popt, pcov = curve_fit(CCF.gmo, tau0, ccfs0, p0=[amp0, mean0, std0], maxfev=5000)
        except RuntimeError:
            print('\n-------------------------------------------------')
            print('RuntimeError: curve_fit error for 1-th pair light curves!')
            print('-------------------------------------------------')
            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.scatter(tau0, ccfs0, marker='+', color='bl', s=20, linewidths=0.5, alpha=1.0)
            ax.set_xlabel('Time Delay (sec)')
            ax.set_ylabel('CCF Value')
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/1-th_tau_ccfs%s.pdf'%self.suffix,
                        bbox_inches='tight', pad_inches=0.1, dpi=100)
            plt.close(fig)
            sys.exit()

        amp, mean, std = popt
        fwhm = 2 * std * np.sqrt(2 * np.log(2))
        range = [max(min(tau0), mean-fwhm/self.nf), min(max(tau0), mean+fwhm/self.nf)]
        idx = (tau0 >= range[0]) * (tau0 <= range[1])

        mc_ccfs, mc_gmos, mc_lags = [], [], []
        for i, (x, y) in enumerate(zip(self.xs, self.ys)):
            tau, ccfs = CCF.correlate(self.dt, x, y)
            tau, ccfs = tau[idx], ccfs[idx]
            mean = np.average(tau, weights=ccfs)
            std = np.sqrt(np.average((tau - mean) ** 2, weights=ccfs))
            g = models.Gaussian1D(amplitude=1, mean=mean, stddev=std)
            amp = np.max(ccfs) / g(mean)

            try:
                popt, pcov = curve_fit(CCF.gmo, tau, ccfs, p0=[amp, mean, std], maxfev=5000)
            except RuntimeError:
                print('\n-------------------------------------------------')
                print('RuntimeError: curve_fit error for %d-th pair light curves!'%(i+1))
                print('-------------------------------------------------')

            tau_interp = np.arange(min(tau), max(tau), 1e-4)
            ccfs_interp = CCF.gmo(tau_interp, *popt)
            peak_idx = np.argmax(ccfs_interp)
            lag = tau_interp[peak_idx]
            mc_ccfs.append(ccfs)
            mc_gmos.append(popt)
            mc_lags.append(lag)

        tau, _ = CCF.correlate(self.dt, self.xs[0], self.ys[0])
        tau = tau[idx]
        true_lag = mc_lags[0]

        mask = sigma_clip(mc_lags, sigma=5, maxiters=5, stdfunc=mad_std).mask
        not_mask = list(map(operator.not_, mask))
        mc_lags = np.array(mc_lags)[not_mask]

        low_lag, upp_lag = np.percentile(mc_lags, [16, 84])
        true_lag_err = np.diff([low_lag, true_lag, upp_lag])
        lagpm = [true_lag, true_lag_err[0], true_lag_err[1]]

        if self.save:
            np.savetxt(self.save_dir + '/delay_times%s.txt'%self.suffix, np.array(tau)+0)
            np.savetxt(self.save_dir + '/mc_ccfs%s.txt'%self.suffix, np.array(mc_ccfs)+0)
            np.savetxt(self.save_dir + '/mc_gmos%s.txt'%self.suffix, np.array(mc_gmos)+0)
            np.savetxt(self.save_dir + '/mc_lags%s.txt'%self.suffix, np.array(mc_lags)+0)
            with open(self.save_dir + '/lagpm%s.json'%self.suffix, 'w') as f_obj:
                json.dump(lagpm, f_obj, indent=4)
            print('\n-------------------------------------------------')
            print('lag calculator:')
            with open(self.save_dir + '/lagpm%s.json'%self.suffix) as f_obj:
                for line in f_obj:
                    print(line.rstrip())
            print('-------------------------------------------------')

            rcParams['mathtext.fontset'] = 'stix'
            rcParams['font.family'] = 'Times New Roman'
            rcParams['font.serif'] = 'Times New Roman'
            rcParams['font.size'] = '16'

            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.scatter(tau, mc_ccfs[0], marker='+', color='black', s=20, linewidths=0.5, alpha=1.0)
            ax.plot(tau, CCF.gmo(tau, *mc_gmos[0]), c='black', lw=0.5, alpha=1.0)
            for i in np.arange(1, len(mc_ccfs), 100):
                ax.scatter(tau, mc_ccfs[i], marker='+', color='blue', s=10, linewidths=0.5, alpha=0.5)
                ax.plot(tau, CCF.gmo(tau, *mc_gmos[i]), c='blue', lw=0.2, alpha=0.5)
            ax.set_xlabel('Time Delay (sec)')
            ax.set_ylabel('CCF Value')
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/tau_ccfs%s.pdf'%self.suffix, bbox_inches='tight', pad_inches=0.1, dpi=100)
            print('\n+++++tau_ccfs%s.pdf+++++'%self.suffix)
            plt.close(fig)

            bins = np.linspace(min(mc_lags), max(mc_lags), 30)
            histvalue, histbin = np.histogram(mc_lags, bins=bins)
            lags = (histbin[1:] + histbin[:-1]) / 2
            counts = histvalue
            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.plot(lags, counts, ds='steps-mid', color='blue')
            ax.axvline(true_lag, c='grey')
            ax.axvline(low_lag, c='grey', ls='--')
            ax.axvline(upp_lag, c='grey', ls='--')
            ax.set_xlabel('lags (sec)')
            ax.set_ylabel('counts')
            ax.set_title(r'$lag=%.2f_{-%.2f}^{+%.2f}\,sec$'%(true_lag, true_lag_err[0], true_lag_err[1]))
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/lags_pdf%s.pdf'%self.suffix, bbox_inches='tight', pad_inches=0.1, dpi=100)
            print('\n+++++lags_pdf%s.pdf+++++'%self.suffix)
            plt.close(fig)

        return lagpm


    def fit_tau_ccf(self, deg=4):
        tau0, ccfs0 = CCF.correlate(self.dt, self.xs[0], self.ys[0])
        peak_idx0 = np.argmax(ccfs0)
        idx = np.arange(peak_idx0 - 5 if peak_idx0 >= 5 else 0, peak_idx0 + 6, 1)

        mc_ccfs, mc_poly, mc_lags = [], [], []
        for i, (x, y) in enumerate(zip(self.xs, self.ys)):
            tau, ccfs = CCF.correlate(self.dt, x, y)
            tau, ccfs = tau[idx], ccfs[idx]

            tau_interp = np.arange(min(tau), max(tau), 1e-4)
            poly_fit = np.polyfit(tau, ccfs, deg=deg)
            ccfs_interp = np.polyval(poly_fit, tau_interp)
            peak_idx = np.argmax(ccfs_interp)
            lag = tau_interp[peak_idx]
            mc_ccfs.append(ccfs)
            mc_poly.append(poly_fit)
            mc_lags.append(lag)

        tau, _ = CCF.correlate(self.dt, self.xs[0], self.ys[0])
        tau = tau[idx]
        tau_interp = np.arange(min(tau), max(tau), 1e-4)
        true_lag = mc_lags[0]

        mask = sigma_clip(mc_lags, sigma=5, maxiters=5, stdfunc=mad_std).mask
        not_mask = list(map(operator.not_, mask))
        mc_lags = np.array(mc_lags)[not_mask]

        low_lag, upp_lag = np.percentile(mc_lags, [16, 84])
        true_lag_err = np.diff([low_lag, true_lag, upp_lag])
        lagpm = [true_lag, true_lag_err[0], true_lag_err[1]]

        if self.save:
            np.savetxt(self.save_dir + '/delay_times%s.txt'%self.suffix, np.array(tau)+0)
            np.savetxt(self.save_dir + '/mc_ccfs%s.txt'%self.suffix, np.array(mc_ccfs)+0)
            np.savetxt(self.save_dir + '/mc_lags%s.txt'%self.suffix, np.array(mc_lags)+0)
            with open(self.save_dir + '/lagpm%s.json'%self.suffix, 'w') as f_obj:
                json.dump(lagpm, f_obj, indent=4)
            print('\n-------------------------------------------------')
            print('lag calculator:')
            with open(self.save_dir + '/lagpm%s.json'%self.suffix) as f_obj:
                for line in f_obj:
                    print(line.rstrip())
            print('-------------------------------------------------')

            rcParams['mathtext.fontset'] = 'stix'
            rcParams['font.family'] = 'Times New Roman'
            rcParams['font.serif'] = 'Times New Roman'
            rcParams['font.size'] = '16'

            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.scatter(tau, mc_ccfs[0], marker='+', color='black', s=20, linewidths=0.5, alpha=1.0)
            ax.plot(tau_interp, np.polyval(mc_poly[0], tau_interp), c='black', lw=0.5, alpha=1.0)
            for i in np.arange(1, len(mc_ccfs), 100):
                ax.scatter(tau, mc_ccfs[i], marker='+', color='blue', s=10, linewidths=0.5, alpha=0.5)
                ax.plot(tau_interp, np.polyval(mc_poly[i], tau_interp), c='blue', lw=0.2, alpha=0.5)
            ax.set_xlabel('Time Delay (sec)')
            ax.set_ylabel('CCF Value')
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/tau_ccfs%s.pdf'%self.suffix, bbox_inches='tight', pad_inches=0.1, dpi=100)
            print('\n+++++tau_ccfs%s.pdf+++++'%self.suffix)
            plt.close(fig)

            bins = np.linspace(min(mc_lags), max(mc_lags), 30)
            histvalue, histbin = np.histogram(mc_lags, bins=bins)
            lags = (histbin[1:] + histbin[:-1]) / 2
            counts = histvalue
            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(111)
            ax.plot(lags, counts, ds='steps-mid', color='blue')
            ax.axvline(true_lag, c='grey')
            ax.axvline(low_lag, c='grey', ls='--')
            ax.axvline(upp_lag, c='grey', ls='--')
            ax.set_xlabel('lags (sec)')
            ax.set_ylabel('counts')
            ax.set_title(r'$lag=%.2f_{-%.2f}^{+%.2f}\,sec$'%(true_lag, true_lag_err[0], true_lag_err[1]))
            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', direction='in')
            ax.tick_params(which='major', width=1.5, length=7)
            ax.tick_params(which='minor', width=1.5, length=4)
            plt.savefig(self.save_dir + '/lags_pdf%s.pdf'%self.suffix, bbox_inches='tight', pad_inches=0.1, dpi=100)
            print('\n+++++lags_pdf%s.pdf+++++'%self.suffix)
            plt.close(fig)

        return lagpm


from_dir, nmc, t1, t2 = sys.argv[1:]
nmc = int(nmc)
t1, t2 = float(t1), float(t2)

e1e2_txt = from_dir + '/e1e2.txt'
lc_list_txt = from_dir + '/lc_list.txt'

e1e2 = np.array(loadtxt(file=e1e2_txt))
ergs = [np.sqrt(e1e2[i][0] * e1e2[i][1]) for i in range(0, len(e1e2))]
ergs_err = [np.array(ergs) - e1e2[:, 0], e1e2[:, 1] - np.array(ergs)]
lc_list = loadtxt(file=lc_list_txt, trans=True)[0]

print('\n-------------------------------------------------')
print('Begin to calculate Lyman lags!')
print('-------------------------------------------------')
save_dir = from_dir + '/Lyman_Lag_Res'
if not isdir(save_dir): os.mkdir(save_dir)

lagpms, lags, lags_err = {}, [0], [[0], [0]]
lc_table = ascii.read(lc_list[0])
time0 = lc_table.field(0)
dt = time0[1] - time0[0]
idx = (time0 >= t1) & (time0 <= t2)
cts0 = (lc_table.field(1) * dt)[idx]

for i in range(1, len(e1e2)):
    pvp_dir = save_dir + '/%02dv%02d'%(i+1, 1)
    if not isdir(pvp_dir): os.mkdir(pvp_dir)

    print('\n-------------------------------------------------')
    print('doing %02d vs. %02d!'%(i+1, 1))
    print('-------------------------------------------------')

    lc_table = ascii.read(lc_list[i])
    cts = (lc_table.field(1) * dt)[idx]

    ccf = CCF(dt, cts, cts0, np.sqrt(cts), np.sqrt(cts0), save=True, save_dir=pvp_dir, suffix='')
    ccf.mc_simulation(nmc)
    lagpm = ccf.fit_tau_ccf()
    lags.append(lagpm[0])
    lags_err[0].append(lagpm[1])
    lags_err[1].append(lagpm[2])
    lagpms['lag-%02d%02d'%(i+1, 1)] = lagpm

    sys.stdout.flush()
    time.sleep(5)

with open(save_dir + '/Erg_Lags.json', 'w') as f_obj:
    json.dump(lagpms, f_obj, indent=4)
savetxt(file=save_dir + '/Erg_Lags.txt', data=[ergs] + ergs_err + [lags] + lags_err, trans=True)

print('\n-------------------------------------------------')
print('Erg vs lags:')
with open(save_dir + '/Erg_Lags.json') as f_obj:
    for line in f_obj:
        print(line.rstrip())
print('\n-------------------------------------------------')

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111)
ax.errorbar(ergs, lags, xerr=ergs_err, yerr=lags_err, ms=2, fmt='ob', elinewidth=1, capsize=1, alpha=1.0)
ax.axhline(0, c='grey', ls='--', lw=0.5)
ax.set_xscale('log')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Lag (s)')
ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in')
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.savefig(save_dir + '/Erg_Lags.pdf', bbox_inches='tight', pad_inches=0.1, dpi=100)
print('\n+++++Erg_Lags.pdf+++++')
plt.close(fig)

print('\n-------------------------------------------------')
print('Calculated Lyman lags!')
print('-------------------------------------------------')

print('\n-------------------------------------------------')
print('Begin to calculate Relay lags!')
print('-------------------------------------------------')
save_dir = from_dir + '/Relay_Lag_Res'
if not isdir(save_dir): os.mkdir(save_dir)

lc_table0 = ascii.read(lc_list[0])
time0 = lc_table0.field(0)
dt = time0[1] - time0[0]
idx = (time0 >= t1) & (time0 <= t2)

lagpms, lags, lags_err = {}, [0], [[0], [0]]
for i in range(1, len(e1e2)):
    pvp_dir = save_dir + '/%02dv%02d'%(i+1, i)
    if not isdir(pvp_dir): os.mkdir(pvp_dir)

    print('\n-------------------------------------------------')
    print('doing %02d vs. %02d!'%(i+1, i))
    print('-------------------------------------------------')

    lc_table1 = ascii.read(lc_list[i-1])
    cts1 = (lc_table1.field(1) * dt)[idx]

    lc_table2 = ascii.read(lc_list[i])
    cts2 = (lc_table2.field(1) * dt)[idx]

    ccf = CCF(dt, cts2, cts1, np.sqrt(cts2), np.sqrt(cts1), save=True, save_dir=pvp_dir, suffix='')
    ccf.mc_simulation(nmc)
    lagpm = ccf.fit_tau_ccf()
    lags.append(lagpm[0])
    lags_err[0].append(lagpm[1])
    lags_err[1].append(lagpm[2])
    lagpms['lag-%02d%02d'%(i+1, i)] = lagpm

    sys.stdout.flush()
    time.sleep(5)

with open(save_dir + '/Erg_Lags.json', 'w') as f_obj:
    json.dump(lagpms, f_obj, indent=4)
savetxt(file=save_dir + '/Erg_Lags.txt', data=[ergs] + ergs_err + [lags] + lags_err, trans=True)

print('\n-------------------------------------------------')
print('Erg vs lags:')
with open(save_dir + '/Erg_Lags.json') as f_obj:
    for line in f_obj:
        print(line.rstrip())
print('\n-------------------------------------------------')

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111)
ax.errorbar(ergs, lags, xerr=ergs_err, yerr=lags_err, ms=2, fmt='ob', elinewidth=1, capsize=1, alpha=1.0)
ax.axhline(0, c='grey', ls='--', lw=0.5)
ax.set_xscale('log')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Lag (s)')
ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in')
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
plt.savefig(save_dir + '/Erg_Lags.pdf', bbox_inches='tight', pad_inches=0.1, dpi=100)
print('\n+++++Erg_Lags.pdf+++++')
plt.close(fig)

print('\n-------------------------------------------------')
print('Calculated Relay lags!')
print('-------------------------------------------------')
