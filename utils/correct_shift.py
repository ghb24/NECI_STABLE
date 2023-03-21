#!/usr/bin/env python2
"""
    Corrects the population control bias of ths shift in FCIQMC
"""
import math
import numpy as np
import pandas as pd

__author__ = "Khaldoon Ghanem"
__email__ = "k.ghanem@fkf.mpg.de"
__status__ = "Dev"

file_name = "./FCIMCStats" #Path to FCIMCStats file
tau = 1e-3 # Time step used in the simulation
steps =  10 # Number of iterations between samples
HF = 0 # Energy of reference determinant
start_iter = 100000 # Number of iterations to be discarded from the beginning (burn-in period)
start_idx = int(start_iter/steps)
end_iter = None # Last Iteration
BLOCK_SIZE = 1024 # block size used to get uncorrelated estimate of the errorbars
min_n = 32 # minimum correction order (as a multiple of steps paramter)
max_n = 1024 # maximum correction order (as a multiple of steps paramter)
log_scale_n = True # whether to vary n on a linear or logarithmic scale
step_n = 1 # steps between correction orders
BLOCK_SIZE_PER_N = 128 # Block size per correction order n. Estimates with larger n have higher correlation, so we use larger blocks

def get_stats(data, block_size = BLOCK_SIZE):
    ndata = len(data)
    nblocks = int(ndata/block_size)
    shift = ndata%block_size
    bdata = np.zeros(nblocks)
    for i in range(nblocks):
        bdata[i] = np.sum(data[shift+i*block_size: shift+(i+1)*block_size])/block_size
    return np.mean(bdata), np.std(bdata)/math.sqrt(nblocks-1)

def get_ratio_stats(num, den, block_size = BLOCK_SIZE):
    ndata = len(num)
    nblocks = int(ndata/block_size)
    shift = ndata%block_size
    bnum = np.zeros(nblocks)
    bden = np.zeros(nblocks)
    for i in range(nblocks):
        bnum[i] = np.sum(num[shift+i*block_size: shift+(i+1)*block_size])/block_size
        bden[i] = np.sum(den[shift+i*block_size: shift+(i+1)*block_size])/block_size

    num_avg = np.mean(bnum)
    num_var = np.var(bnum)/(nblocks-1)
    den_avg = np.mean(bden)
    den_var = np.var(bden)/(nblocks-1)
    cov = np.cov(bnum, bden)[0][1]/(nblocks)

    ratio = num_avg/den_avg
    ratio_rel_err = math.sqrt(num_var/num_avg**2 + den_var/den_avg**2 - 2*cov/(num_avg*den_avg))
    return ratio, ratio_rel_err, nblocks

def blocking(data):
    ndata = len(data)
    max_log_block_size = int(math.log2(ndata))-1
    for block_size in np.logspace(0, max_log_block_size, max_log_block_size+1, base=2.0):
        mean, err = get_stats(shift, block_size)
        nblocks = int(ndata/block_size)
        print("{}\t{}\t{}\t{:e}".format(block_size, nblocks, mean, err))


if(end_iter):
    end_idx = int(end_iter/steps)
else:
    end_idx = -1
#data = np.loadtxt(file_name, usecols = (1, 4))[start_idx:end_idx,:]
data = pd.read_csv(file_name, delim_whitespace=True, comment='#', usecols=[1,4]).values[start_idx:end_idx,:]
nsamples = data.shape[0]
shift = data[:, 0]
nwalkers = data[:, 1]

print("================")
print("Basic Estimates:")
print("================")
print("Number of Blocks:  {}".format(nsamples/BLOCK_SIZE))
shift_avg, shift_avg_err = get_stats(shift)
#walkers_avg = np.mean(nwalkers)
nwalkers_avg, nwalkers_avg_err = get_stats(nwalkers)
#walkers_var = np.var(nwalkers)
nwalkers_var, nwalkers_var_err = get_stats((nwalkers-nwalkers_avg)**2)
#shift_var = np.var(shift)
shift_var, shift_var_err = get_stats((shift-shift_avg)**2)
#cov = np.cov([shift[:-1], nwalkers_data[1:]])[0,1]
cov, cov_err = get_stats((shift[:-1]-shift_avg)*(nwalkers[1:]-nwalkers_avg))

# This estimates the uniform projected energy
ratio = cov/nwalkers_avg
upe = shift_avg+ratio
ratio_err = ratio*math.sqrt((cov_err/cov)**2+(nwalkers_avg_err/nwalkers_avg)**2)
upe_err = math.sqrt(shift_avg_err**2+ratio_err**2)

print("Num. Walkers Avg:  {:.5f}+/-{:.3e}".format(nwalkers_avg, nwalkers_avg_err))
print("Num. Walkers Var:  {:.5f}+/-{:.3e}".format(nwalkers_var, nwalkers_var_err))
print("Shift-Walkers Cov: {:.5f}+/-{:.3e}".format(cov, cov_err))
print("Shift Var:         {:.5f}+/-{:.3e}".format(shift_var, shift_var_err))
print("")
print("Shift Avg:         {:.5f}+/-{:.3e}".format(HF+shift_avg,shift_avg_err))
print("UPE Approx.:       {:.5f}+/-{:.3e}".format(HF+upe,upe_err))

print("==========================")
print("Corrected Shift Estimates:")
print("==========================")
shifts_cum = np.cumsum(data[:,0])
shift_const = shift_avg

if(log_scale_n):
    min_log2_n = int(math.log(min_n, 2))
    max_log2_n = int(math.log(max_n, 2))
    n_vals = [2**log2_n for log2_n in range(min_log2_n, max_log2_n+step_n, step_n)]
else:
    n_vals = list(range(min_n, max_n+step_n, step_n))

for n in n_vals:
    shifts_sum = shifts_cum[n:-1]-shifts_cum[:-n-1]
    weights = np.exp(-tau*steps*(shifts_sum-n*shift_const))

    shifts_sum1= shifts_cum[n+1:]-shifts_cum[:-n-1]
    weights1 = np.exp(-tau*steps*(shifts_sum1-(n+1)*shift_const))

    num = weights1*nwalkers[n+1:]
    den = weights*nwalkers[n:-1]
    if(2*BLOCK_SIZE_PER_N*n>nsamples):
        print("Not enough samples for blocking using n={}".format(n*steps))
        break
    ratio, ratio_rel_err, nblocks = get_ratio_stats(num, den, BLOCK_SIZE_PER_N*n)
    correction = -math.log(ratio)/(steps*tau)
    err = ratio_rel_err/(steps*tau)
    energy = shift_const+correction

    print("Energy (tau = {:.1e}, n={:5d}, blocks={:5d}): {:.5f}+/-{:.3e}".format(n*tau*steps, n*steps, nblocks, HF+energy, err))

