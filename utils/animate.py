import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import math
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from bisect import bisect_left, bisect_right
import sys
import glob
import argparse
import traceback
import warnings
import sys

###############################################################################
#Parse arguments
###############################################################################

parser = argparse.ArgumentParser(description= "Animate NECI wavefunctions \
                                 alongside energy, reference population and \
                                 total number of walkers.")

parser.add_argument('-s','--save', help="Save animation into a gif file.",
                    type=str, metavar='FILE_NAME', dest="file_name")
parser.add_argument('--itr-first', help="First plotted iteration.",
                    type=int)
parser.add_argument('--itr-last', help="Last plotted iteration.",
                    type=int)
parser.add_argument('--itr-step', help="Steps between plotted iterations.",
                    type=int)
parser.add_argument('--interval', help="Display duration of an iteration \
                    (in miliseconds).", type=int,default=100)
parser.add_argument('--ex-low', help="Lowest plotted excitation level.",
                    type=int,default=1)
parser.add_argument('--ex-high', help="Highest plotted excitation level.",
                    type=int)
parser.add_argument('--ex-min-dets', help="Minimum number of determinants\
                    in an excitation level in order to be plotted.",
                    type=int,default=8)
parser.add_argument('--ex-scale', help="Scaling factor for excitation levels.\
                    Lower values make the ampitudes smaller but the \
                    seperation between levels clearer.",
                    type=float,default=0.5)
parser.add_argument('--wf-min', help="The lower limit of the y-axis of the \
                    wavefunction plot.",type=float)
parser.add_argument('--wf-max', help="The upper limit of the y-axis of the \
                    wavefunction plot.",type=float)
args=parser.parse_args()
###############################################################################
#Setup parameters
###############################################################################
#File names:
#-----------
symdets = "SymDETS"
fcimcstats = "FCIMCStats"
energies = "ENERGIES"
spawns_prefix = "SpawnHist-"

#Size of figure
#---------------
fig_size = (16,10)

#Layout of different subplots inside the figure:
#-----------------------------------------------
G = gridspec.GridSpec(4, 2)
excits_block = G[:3,1]
fullWF_block = G[3,1]
E_block = G[:2,0]
N0_block = G[2,0]
N_block = G[3,0]

#Line colors:
#------------
#Exact wavefunction
exact_color = 'blueviolet'
#Average wavefunction
avg_color = 'seagreen'
#Instantenous wavefunction
instant_color = 'gold'
#Reference population
N0_color = 'b'
#Total number of walkers
N_color = 'g'
#Projected energy
PE_color = 'r'
#Exact energy
EE_color = 'k'
###############################################################################
#Read FCIMCStats file
###############################################################################
itrs = []
N0 = []
N = []
PE = []
with open(fcimcstats,"r") as f:
    lines=f.readlines()
    for l in lines:
        if(l.strip().startswith('#')):
            continue
        columns = l.split()
        itrs.append(int(columns[0]))
        N0.append(float(columns[11]))
        N.append(float(columns[4]))
        PE.append(float(columns[22]))

###############################################################################
#Read exact energy
###############################################################################
with open(energies,"r") as f:
    line=f.readline()
    exact_E = float(line)

###############################################################################
#Read exact wave function and excitation levels
###############################################################################
with open(symdets,"r") as f:
    lines=f.readlines()
    #Sequential identifies of determinants
    ids = []
    #Populations of exact wave function
    exact_wf=[]
    #Excitation levels of the determinants
    excits=[]
    for l in lines:
        if(l.strip().startswith('#')):
            continue
        columns = l.split()
        ids.append(int(columns[0]))
        excits.append(int(columns[1]))
        exact_wf.append(float(columns[3]))

#For full wavefunction limits, ignore HF because its population is very large
if(args.wf_min):
    min_y = args.wf_min
else:
    min_y = min(exact_wf[1:])*1.1
if(args.wf_max):
    max_y = args.wf_max
else:
    max_y = max(exact_wf[1:])*1.1

min_ex = max(1, args.ex_low)
if(args.ex_high):
    max_ex = min(args.ex_high, max(excits))
else:
    max_ex = max(excits)
if(min_ex>max_ex):
    sys.exit('Lowest excitation level should not exceed the highest one!')
###############################################################################
#Extract the index of first determinant of each excitation level
###############################################################################
ex_idx=np.zeros(max_ex+2, dtype=int)
idx = 0
for ex in range(max_ex+2):
    while(idx<len(excits) and excits[idx]<ex):
        idx = idx + 1
    ex_idx[ex] = idx
###############################################################################
#Routines to convert wavefunction to polar coordinates
###############################################################################
def x2angle(ids, ex):
    begin = ex_idx[ex]
    end = ex_idx[ex+1]+1
    res = [(2*math.pi*(x-begin)/(end-begin)) for x in ids[begin:end]]
    #repeat first point to close the loop
    res.append(res[0])
    return res

def y2radius(wf, ex):
    begin = ex_idx[ex]
    end = ex_idx[ex+1]+1
    #renormalize each excitation level individually
    #and multipy by ex_factor to make a clearer visual seperation of levels
    max_y = max(map(abs,exact_wf[begin:end]))
    res =  [(ex+(y*args.ex_scale/max_y)) for y in wf[begin:end]]
    #repeat first point to close the loop
    res.append(res[0])
    return res

###############################################################################
#Initialize figure
###############################################################################
fig = plt.figure(figsize=fig_size)

#Intitialize the full plot of wavefunctions
#-----------------------------------------
fullWF_ax = fig.add_subplot(fullWF_block)
exact_line, = fullWF_ax.plot(ids, exact_wf, label="Exact WF", 
                             color=exact_color)
instant_line, = fullWF_ax.plot(ids, [0]*len(ids), label="Instant. WF", 
                               color=instant_color)
avg_line, = fullWF_ax.plot(ids, [0]*len(ids), label="Avg. WF", 
                           color=avg_color)
fullWF_ax.set_ylim(min_y,max_y)
fullWF_ax.set_xlabel("Determinant")

#Intitialize the excitations plot of wavefunctions
#--------------------------------------------------
excits_ax = fig.add_subplot(excits_block, polar=True)
exact_plines = dict()
instant_plines = dict()
avg_plines = dict()
for ex in range(min_ex, max_ex+1):
    angles = x2angle(ids,ex)
    #Skip levels with very few determinants
    if(len(angles)<args.ex_min_dets):
        continue
    exact_plines[ex], = excits_ax.plot(angles, y2radius(exact_wf, ex), 
                                      color=exact_color)
    instant_plines[ex], = excits_ax.plot(angles, [ex]*len(angles), 
                                        color=instant_color)
    avg_plines[ex], = excits_ax.plot(angles, [ex]*len(angles), 
                                    color=avg_color)

excits_ax.plot([], [], label="Exact WF", color=exact_color)
excits_ax.plot([], [], label="Instant. WF", color=instant_color)
excits_ax.plot([], [], label="Avg. WF", color=avg_color)
excits_ax.set_title("Iteration: ")
excits_ax.legend(loc='upper right', bbox_to_anchor=(1.1,1.1))
excits_ax.axes.get_xaxis().set_visible(False)
excits_ax.set_rticks(range(min_ex,max_ex+1))

#Initialize energy plot
#----------------------
E_ax = fig.add_subplot(E_block)
EE_line, = E_ax.plot(itrs, [exact_E]*len(itrs), label="Exact Energy", 
                     color=EE_color)
PE_line, = E_ax.plot(itrs, [0]*len(itrs), label="Projected Energy", 
                     color=PE_color)
E_ax.set_ylim(min(PE),max(PE))
plt.setp(E_ax.get_xticklabels(), visible=False)

#Initialize reference population plot
#------------------------------------
N0_ax = fig.add_subplot(N0_block, sharex=E_ax)
N0_line, = N0_ax.plot(itrs, [0]*len(itrs), color=N0_color)
N0_ax.set_ylim(min(N0),max(N0)*1.1)
plt.setp(N0_ax.get_xticklabels(), visible=False)

#Initialize total walkers plot
#-----------------------------
N_ax = fig.add_subplot(N_block, sharex=N0_ax)
N_line, = N_ax.plot(itrs, [0]*len(itrs), color=N_color)
N_ax.set_ylim(min(N),max(N)*1.1)
N_ax.set_xlabel("Iteration")

E_ax.legend((EE_line, PE_line, N0_line, N_line), 
            ('Exact Energy', 'Projected Energy', 'Reference Population',
              'Total Walkers'), loc='upper right')

#Let python decide on the optimal spacing between plots
plt.tight_layout(pad=1.3)
###############################################################################
#Animation routines
###############################################################################
#Read the average and instantanous wave functions:
#-------------------------------------------------
def get_wfs(itr):
    with open(spawns_prefix+str(itr),"r") as f:
        lines=f.readlines()
        instant_wf=[]
        avg_wf=[]
        for l in lines:
            if(l.strip().startswith('#')):
                continue
            columns = l.split()
            avg_wf.append(float(columns[1]))
            instant_wf.append(float(columns[3]))
        #Ensure same sign convention by comparing the sign of HF
        if(instant_wf[0]*exact_wf[0]<0):
            instant_wf = [-x for x in instant_wf]
            avg_wf = [-x for x in avg_wf]
        return instant_wf, avg_wf

#Dummy initialization routine:
#----------------------------
def init():
    return

#Update figure:
#--------------
def update(itr):
    print "Iteration: "+str(itr)
    excits_ax.set_title("Iteration: "+str(itr))

    #Get wavefuntions of this iteration
    instant_wf, avg_wf = get_wfs(itr)

    #Update full wavefunction plot
    instant_line.set_ydata(instant_wf)
    avg_line.set_ydata(avg_wf)

    #Update excitations plot
    for ex in range(min_ex, max_ex+1):
        if(ex in instant_plines):
            instant_plines[ex].set_ydata(y2radius(instant_wf, ex))
            avg_plines[ex].set_ydata(y2radius(avg_wf, ex))

    #Find the index of an iteration
    itr_idx = bisect_left(itrs, itr)

    #Update projected energy
    PE_line.set_xdata(itrs[:itr_idx])
    PE_line.set_ydata(PE[:itr_idx])

    #Update reference population
    N0_line.set_xdata(itrs[:itr_idx])
    N0_line.set_ydata(N0[:itr_idx])

    #Update total walkers
    N_line.set_xdata(itrs[:itr_idx])
    N_line.set_ydata(N[:itr_idx])

###############################################################################
#Generate animation
###############################################################################
#Figure out the avalaible iterations automatically
spawnshist_files  = glob.glob(spawns_prefix+"*")
plot_itrs = [int(s.replace(spawns_prefix, "")) for s in  spawnshist_files]
if(len(plot_itrs)==0):
    sys.exit('No SpawnHist files availible!')
plot_itrs.sort()
#Now modify them according to input
if(args.itr_first):
    itr_idx = bisect_left(plot_itrs, args.itr_first)
    plot_itrs = plot_itrs[itr_idx:]

if(args.itr_last):
    itr_idx = bisect_right(plot_itrs, args.itr_last)
    plot_itrs = plot_itrs[:itr_idx]

if(len(plot_itrs)==0):
    sys.exit('No SpawnHist files matching specified iterations are availible!')

if(args.itr_step):
    first_itr = plot_itrs[0]
    plot_itrs = [itr for itr in plot_itrs if (itr-first_itr)%args.itr_step==0]

if(len(plot_itrs)==0):
    sys.exit('No SpawnHist files matching specified iterations are availible!')

#Create animation object
anim = FuncAnimation(fig, update, plot_itrs, init, interval=args.interval)

#Without the following, an annoying warning is printed. After tracking this 
#strange warning, I found it happens when using a polar plot with GridSpec. 
#It seems like an internal issue of matplotlib! 
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    if(args.file_name):
        anim.save(args.file_name, dpi=80, writer='imagemagick')
    else:
        plt.show()
