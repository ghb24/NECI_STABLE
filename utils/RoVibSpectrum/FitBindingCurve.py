#!/usr/bin/env python

#To calculate the rotational/vibrational from a binding curve of a diatomic the following inputs are needed:
#RESULTS file:
#Containing the distance between nuclei (in bohr) in the first column, and the total energy relative to the
#completely dissociated limit in the second.
#inputVib file:
#Example of which is given in this folder
#This script is then run (it calls CalcVibSpectrum.x).
#The method followed can be found in Bytautas et al. J.Chem.Phys 127 (2007)

import sys,os
from pylab import *
from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
from numpy import array,mean,std

params = {'legend.fontsize': 11}
rcParams.update(params)

dirname = "PlotsofFittings"
if not os.path.isdir("./" + dirname + "/"):
    os.mkdir("./" + dirname + "/")

print "\n** FITTING THE BINDING CURVE TO A GAUSSIAN EXPANSION POTENTIAL **"
print "These have the form g(R) = sum_k a_k exp (-alpha beta^k R^2)"
print "Data is read in from the RESULTS file with the R in bohr, and g(R) = E_tot(dimer) - 2E_tot(atom) in hartrees."

def residuals(a, y, x, alpha, beta):
    err = y-peval(x,a,alpha,beta)
    return err

def residualsexp(expon, y, x, coeffs):
    err = y-peval(x,coeffs,expon[0],expon[1])
    return err

# p is the fitting function parameters from 0 -> K
def peval(x, a, alpha, beta):
    fit=0
    for i in range(0, len(a)):
#        print i
        fit=fit+a[i]*exp(-alpha*(beta**i)*(x**2))   #gaussian fit
    return fit

filename='RESULTS'
data = scipy.io.array_import.read_array(filename)

y = data[:,2]   #these are the energies
x = data[:,0]   #these are the bond lengths
Etots = data[:,1]
n=len(Etots)

#print x, y

#All parameters are given some initial guesses. 
#E_0=0.0   #this one wants to be the energy at dissociation
Minx=2.154287627
Maxx=15.117807909
k=4
a0_0=-0.1788169627
a1_0=-1.5169391210
a2_0=5.0621415634
a3_0=-2.9439346580
a4_0=16.719611821
#alpha=0.82
#beta=1.55
#a0_0=-0.1
#a1_0=-1.5
#a2_0=5.0
#a3_0=-2.9
#a4_0=16.7
alpha=0.82/3.5710643135020281
beta=1.55

#An array with the names of the parameters is created for printing the results and all initial guesses are also stored in an array. 
aname = (['a_0','a_1','a_2','a_3','a_4'])
a = array([a0_0 , a1_0, a2_0, a3_0, a4_0])

a_0 = array([a0_0 , a1_0, a2_0, a3_0, a4_0])
alpha_0=alpha
beta_0=beta

print "First optimising the a_k values..."
print "Initial parameters"
for i in range(len(aname)):
    print "%s = %.4f " % (aname[i], a[i])

(fit, ierr) = leastsq(residuals, a, args=(y, x, alpha, beta), maxfev=2000)

print "Success (1-4): "
print ierr

print "Final parameters"
for i in range(len(aname)):
    print "%s = %.4f " % (aname[i], fit[i])

#print "Value of function at R=1000:"
#print peval(1000,fit,alpha,beta)

print "Now optimising alpha and beta..."

exponents = array([alpha,beta])

print "Initial parameters"
#print "%s = %.4f " % ('E_0', exponents[0])
print "%s = %.4f " % ('alpha', exponents[0])
print "%s = %.4f " % ('beta', exponents[1])

(fit2, ierr) = leastsq(residualsexp, exponents, args=(y, x, fit), maxfev=2000)

print "Success (1-4): "
print ierr

print "Final parameters"
#print "%s = %.4f " % ('E_0', fit2[0])
print "%s = %.4f " % ('alpha', fit2[0])
print "%s = %.4f " % ('beta', fit2[1])

#a[0]=fit2[0]
alpha=fit2[0]
beta=fit2[1]

errors = residuals(fit, y, x, alpha, beta)
initerrors = residuals(a_0, y, x, alpha_0, beta_0)
#av = mean(errors)
#RMS = sqrt(mean((errors-av)**2))

print 'Initial RMS (mEh)',std(initerrors)*1000
print 'Final RMS (mEh)', std(errors)*1000
#print errors,av,RMS,std(errors),errors-av

figure()
clf()
ax=gca()
plot(x,y*1000,linewidth=0,marker='o',label='orig points')
plot(arange(Minx,Maxx,0.0001),peval(arange(Minx,Maxx,0.0001),fit,alpha,beta)*1000,label='fit')
plot(arange(Minx,Maxx,0.0001),peval(arange(Minx,Maxx,0.0001),a_0,alpha_0,beta_0)*1000,label='orig_guess')
ax.legend(loc=1, markerscale=0.1)
savefig('PlotsofFittings/PotentialFit.eps')
#show()
print 'Saved plot of potential to PotentialFit.eps in the PlotsofFittings folder'

Requil=nanargmin(peval(arange(Minx,Maxx,0.000001),fit,alpha,beta))
Requil=Requil*0.000001+Minx
print "Equilibrium bond length = ", Requil


f = open('PotData', 'w')
f.write (" %i      %s " % (k,'kmax\n'))
f.write (" %.10f      %s " % (Requil,'R_equilibrium (bohr)\n'))
f.write (" %.10f      %s " % (alpha,'alpha\n'))
f.write (" %.10f      %s " % (beta,'beta\n'))
f.write (" %.10f      %s " % (Etots[n-1],'Dissociated Energy\n'))
for i in range(len(aname)):
    f.write ("%.10f      %s \n" % (fit[i], aname[i]))
f.close()

print "\nCalling CalcVibSpectrum program..."
print "This calculates the pure vibrational spectrum, which is printed in the PUREVIB file."
print "A plot of the vibrational energy levels and wavefunction solutions can be viewed by loading GnuplotVib.gpi"
print "A ROTVIBSPECTRUM file is also printed which includes the rotational as well as vibrational energy levels."
print "A plot of the full Ro-Vibrational spectrum can be viewed by loading GnuplotRotVib.gpi"
os.system("./CalcVibSpectrum.x") 

#============== Calculating spectral information ===============================================================
#Out of the above program, we get data so that we can fit various things to find the rotational parameters etc.

print "\n** CALCULATING SPECTRAL INFORMATION **\n"
print "Step 1/2: Values for Bv and Dv are first found using the following equation:"
print "Fv(J)/J(J+1) = Bv - Dv[J(J+1)]"
print "(the required values have been printed to ROTDATA in the above program)"
dirname = "BvDvPlots"
if not os.path.isdir("./" + dirname + "/"):
    os.mkdir("./" + dirname + "/")

inputfile = open('inputVib','r')
line = inputfile.readline()
Maxv = int(line.split(' ')[0])
inputfile.close()

filename='ROTDATA'
data = scipy.io.array_import.read_array(filename)

BvDvFile = open('BvDvValues','w')
BvDvFile.write ("%s     %s               %s               %s            %s\n" % ('# v','Bv (E_h)','Dv (E_h)','Bv (cm-1)','Dv (10^-6 cm-1)'))

for v in range(0,Maxv):

    def residualsBvDv(BvDv,Fv, J):
        err = Fv-peval(J,BvDv)
        return err

    # p is the fitting function parameters from 0 -> K
    def peval(J, BvDv):
        fit3=0
#        fit3=(BvDv[0]*(J*(J+1)))-(BvDv[1]*(J*((J+1)**2)))
        fit3=BvDv[0]-(BvDv[1]*(J*(J+1)))
        return fit3

    print "Fitting Fv vs J curve for v value: ",v

    Fv = data[:,(v+1)]      #these are the Fv (y values)
    J = data[:,0]           #these are the J values
#    for jiter in range(len(J)):
#        Fv[jiter]=Fv[jiter]/(jiter*(jiter+1))

    #print x, y

    #All parameters are given some initial guesses. 
    Bv=0.8833/219474.63
    Dv=(3.4E-06)/219474.63

    BvDv=array([Bv,Dv])
    BvDv_0=array([Bv,Dv])

    print "Initial parameters"
    print "%s = %.15f " % ('Bv',Bv )
    print "%s = %.15f " % ('Dv',Dv )

    (fit3, ierr) = leastsq(residualsBvDv, BvDv, args=(Fv, J), maxfev=2000)

    print "Success Bv / Dv: "
    print ierr

    print "Final parameters"
    print "%s = %.15f " % ('Bv', fit3[0])
    print "%s = %.15f " % ('Dv', fit3[1])

    Bv=fit3[0]
    Dv=fit3[1]

    errors = residualsBvDv(fit3, Fv, J)
    initerrors = residualsBvDv(BvDv_0, Fv, J)
    #av = mean(errors)
    #RMS = sqrt(mean((errors-av)**2))

    print 'Initial RMS (mEh)',std(initerrors)*1000
    print 'Final RMS (mEh)', std(errors)*1000
    #print errors,av,RMS,std(errors),errors-av

    plotname="BvDvPlots/F_"
    plotname=plotname+str(v)
    plotname=plotname+".eps"

    figure()
    clf()
    ax=gca()
    plot(J,Fv,linewidth=0,marker='o',label='orig points')
    plot(arange(0,11,0.0001),peval(arange(0,11,0.0001),fit3),label='fit')
    plot(arange(0,11,0.0001),peval(arange(0,11,0.0001),BvDv_0),label='orig_guess')
    ax.legend(loc=1, markerscale=0.1)
    savefig(plotname)
    #show()

    BvDvFile.write (" %i      %.15f      %.15f      %.10f         %.10f\n" % (v,fit3[0],fit3[1],(fit3[0]*219474.63),((fit3[1]*219474.63)*(10**6))))

BvDvFile.close()
print "Plots of these Fv vs J fits for each value of v are saved in the folder named BvDvPlots"    
print "Bv and Dv values for each v are printed in the file BvDvValues\n"


print "Step 2/2: Using G_v to find the Dunham coefficients of the form Y_k0, and omega_e (w_e)."
print "The equation is of the form G_v = Sum_k Y_k0 ( v + 1/2)^k = Y_00 + omega_e(v + 1/2) - ..."

filename='PUREVIB'
data = scipy.io.array_import.read_array(filename)

def residuals(Yk_0,Gv, v):
    err = Gv-peval(v,Yk_0)
    return err

# p is the fitting function parameters from 0 -> K
def peval(v, Yk_0):
    fit=0
    for i in range(0, len(Yk_0)):
        fit=fit+Yk_0[i]*((v+0.5)**i)
    return fit

# This range of v from 0:8 is chosen because the paper finds this to be appropriate.
# May need to be adjusted!!!
Gv = data[0:8,1]            #these are the Gv (y values)
v = data[0:8,0]             #these are the v values

#All parameters are given some initial guesses. 
Y0_0=-9.5136280671711353e-07
Y1_0=0.0041695115285078736
Y2_0=-5.0311054175145435e-05
Y3_0=-4.5198845989625311e-07
Y4_0=0.0000
Y5_0=0.0000
Y6_0=0.0000

Yk0name = (['Y_00','Y_10','Y_20','Y_30','Y_40','Y_50','Y_60'])
Yk_0 = array([Y0_0 , Y1_0, Y2_0, Y3_0, Y4_0, Y5_0, Y6_0])
Yk_0init = array([Y0_0 , Y1_0, Y2_0, Y3_0, Y4_0, Y5_0, Y6_0])

print "Fitting Gv vs v curve for k up to: ",len(Yk_0)-1
print "The range of v's considered is 0 to ",len(v)

print "Initial parameters"
for i in range(len(Yk0name)):
    print "%s = %.10f " % (Yk0name[i], Yk_0[i])

(fit4, ierr) = leastsq(residuals, Yk_0, args=(Gv, v), maxfev=2000)

print "Success Yk_0: "
print ierr

for i in range(len(Yk0name)):
    Yk_0[i]=fit4[i]

print "Final parameters"
for i in range(len(Yk0name)):
    print "%s = %.10f " % (Yk0name[i], Yk_0[i])

errors = residuals(fit4, Gv, v)
initerrors = residuals(Yk_0init, Gv, v)
#av = mean(errors)
#RMS = sqrt(mean((errors-av)**2))

print 'Initial RMS (mEh)',std(initerrors)*1000
print 'Final RMS (mEh)', std(errors)*1000
#print errors,av,RMS,std(errors),errors-av

figure()
clf()
ax=gca()
plot(v,Gv,linewidth=0,marker='o',label='orig points')
plot(arange(0,9,0.0001),peval(arange(0,9,0.0001),fit4),label='fit')
plot(arange(0,9,0.0001),peval(arange(0,9,0.0001),Yk_0init),label='orig_guess')
ax.legend(loc=1, markerscale=0.1)
savefig('PlotsofFittings/Gv_vs_v.eps')
#show()

YkmFile = open('YkmDunhamCoeffs','w')
YkmFile.write ("%s     %s             %s                     %s                   %s\n" % ('# m','Y_0m (cm-1)','Y_1m','Y_2m','Y_3m'))

YkmFile.write (" %i      %.15f      %.15f      %.10f         %.10f\n" % (0,(fit4[0]*219474.63),(fit4[1]*219474.63),(fit4[2]*219474.63),(fit4[3]*219474.63)))

print "This fit gives an omega_e value of: ",(fit4[1]*219474.63)
print "This fit gives an omega_e.x_e value of: ",(fit4[2]*(-1.0)*219474.63)
print "The fit is saved to Gv_vs_v.eps in the PlotsofFittings folder, and the coefficients and calculated constants printed in YkmDunhamCoeffs"

print "The final, main spectroscopic constants are printed in the SPEC_CONSTS file."

FinalFile = open('SPEC_CONSTS','a')
FinalFile.write ("%s        %.1f      %s     \n" % ('omega_e (in cm-1)',(fit4[1]*219474.63),'The harmonic frequency.'))
FinalFile.write ("%s     %.2f      %s     \n" % ('omega_e.x_e (in cm-1)',(fit4[2]*(-1.0)*219474.63),'The anharmonicity.'))
FinalFile.close()


#print "\nStep 3/4: Using B_v values calculated in step 1 to find the Dunham coefficients of the form Y_k1, B_e and alpha_e."
#print "The equation is of the form B_v = Sum_k Y_k1 ( v + 1/2)^k = B_e - alpha_e(v + 1/2) + ..."

#filename='BvDvValues'
#data = scipy.io.array_import.read_array(filename)

#def residuals(Yk_1,Bv, v):
#    err = Bv-peval(v,Yk_1)
#    return err

## p is the fitting function parameters from 0 -> K
#def peval(v, Yk_1):
#    fit=0
#    for i in range(0, len(Yk_1)):
#        fit=fit+Yk_1[i]*((v+0.5)**i)
#    return fit

## This range of v from 0:8 is chosen because the paper finds this to be appropriate.
## May need to be adjusted!!!
#Bv = data[0:8,1]            #these are the Bv (y values)
#v = data[0:8,0]             #these are the v values

##All parameters are given some initial guesses. 
#Y0_1=-9.5136280671711353e-07
#Y1_1=0.0041695115285078736
#Y2_1=-5.0311054175145435e-05
#Y3_1=-4.5198845989625311e-07
#Y4_1=0.0000
#Y5_1=0.0000
#Y6_1=0.0000

#Yk1name = (['Y_01','Y_11','Y_21','Y_31','Y_41','Y_51','Y_61'])
#Yk_1 = array([Y0_1 , Y1_1, Y2_1, Y3_1, Y4_1, Y5_1, Y6_1])
#Yk_1init = array([Y0_1 , Y1_1, Y2_1, Y3_1, Y4_1, Y5_1, Y6_1])

#print "Fitting Bv vs v curve for k up to: ",len(Yk_1)-1
#print "The range of v's considered is 0 to ",len(v)

#print "Initial parameters"
#for i in range(len(Yk1name)):
#    print "%s = %.10f " % (Yk1name[i], Yk_1[i])

#(fit5, ierr) = leastsq(residuals, Yk_1, args=(Bv, v), maxfev=2000)
#
#print "Success Yk_1: "
#print ierr

#for i in range(len(Yk1name)):
#    Yk_1[i]=fit5[i]
#
#print "Final parameters"
#for i in range(len(Yk1name)):
#    print "%s = %.10f " % (Yk1name[i], Yk_1[i])

#errors = residuals(fit5, Bv, v)
#initerrors = residuals(Yk_1init, Bv, v)
##av = mean(errors)
##RMS = sqrt(mean((errors-av)**2))

#print 'Initial RMS (mEh)',std(initerrors)*1000
#print 'Final RMS (mEh)', std(errors)*1000
##print errors,av,RMS,std(errors),errors-av

#figure()
#clf()
#ax=gca()
#plot(v,Bv,linewidth=0,marker='o',label='orig points')
#plot(arange(0,9,0.0001),peval(arange(0,9,0.0001),fit5),label='fit')
#plot(arange(0,9,0.0001),peval(arange(0,9,0.0001),Yk_1init),label='orig_guess')
#ax.legend(loc=1, markerscale=0.1)
#savefig('PlotsofFittings/Bv_vs_v.eps')
##show()

#YkmFile.write (" %i      %.15f       %.15f       %.10f          %.10f\n" % (1,(fit5[0]*219474.63),(fit5[1]*219474.63),(fit5[2]*219474.63),(fit5[3]*219474.63)))
#
#print "This fit gives a B_e value of: ",(fit5[0]*219474.63)
#print "and an alpha_e value of: ",(fit5[1]*219474.63*(-1.0))
#print "The fit is saved to Bv_vs_v.eps in the PlotsofFittings folder, and the coefficients and calculated constants printed in YkmDunhamCoeffs"
#
#print "\nStep 4/4: Using D_v values calculated in step 1 to find the Dunham coefficients of the form Y_k2, D_e and beta_e."
#print "The equation is of the form D_v = Sum_k Y_k2 ( v + 1/2)^k = D_e + beta_e(v + 1/2) + ..."
#
#filename='BvDvValues'
#data = scipy.io.array_import.read_array(filename)
#
#def residuals(Yk_2,Dv, v):
#    err = Dv-peval(v,Yk_2)
#    return err
#
## p is the fitting function parameters from 0 -> K
#def peval(v, Yk_2):
#    fit=0
#    for i in range(0, len(Yk_2)):
#        fit=fit+Yk_2[i]*((v+0.5)**i)
#    return fit

## This range of v from 0:8 is chosen because the paper finds this to be appropriate.
## May need to be adjusted!!!
#Dv = data[0:8,2]            #these are the Dv (y values)
#v = data[0:8,0]             #these are the v values
#
##All parameters are given some initial guesses. 
#Y0_2=-9.5136280671711353e-07
#Y1_2=0.0041695115285078736
#Y2_2=-5.0311054175145435e-05
#Y3_2=-4.5198845989625311e-07
#Y4_2=0.0000
#Y5_2=0.0000
#Y6_2=0.0000
#
#Yk2name = (['Y_02','Y_12','Y_22','Y_32','Y_42','Y_52','Y_62'])
#Yk_2 = array([Y0_2 , Y1_2, Y2_2, Y3_2, Y4_2, Y5_2, Y6_2])
#Yk_2init = array([Y0_2 , Y1_2, Y2_2, Y3_2, Y4_2, Y5_2, Y6_2])
#
#print "Fitting Dv vs v curve for k up to: ",len(Yk_2)-1
#print "The range of v's considered is 0 to ",len(v)
#
#print "Initial parameters"
#for i in range(len(Yk2name)):
#    print "%s = %.10f " % (Yk2name[i], Yk_2[i])
#
#(fit6, ierr) = leastsq(residuals, Yk_2, args=(Dv, v), maxfev=2000)
#
#print "Success Yk_2: "
#print ierr

#for i in range(len(Yk2name)):
#    Yk_2[i]=fit6[i]
#
#print "Final parameters"
#for i in range(len(Yk2name)):
#    print "%s = %.10f " % (Yk2name[i], Yk_2[i])
#
#errors = residuals(fit6, Dv, v)
#initerrors = residuals(Yk_2init, Dv, v)
##av = mean(errors)
##RMS = sqrt(mean((errors-av)**2))
#
#print 'Initial RMS (mEh)',std(initerrors)*1000
#print 'Final RMS (mEh)', std(errors)*1000
##print errors,av,RMS,std(errors),errors-av
#
#figure()
#clf()
#ax=gca()
#plot(v,Dv,linewidth=0,marker='o',label='orig points')
#plot(arange(0,9,0.0001),peval(arange(0,9,0.0001),fit6),label='fit')
#plot(arange(0,9,0.0001),peval(arange(0,9,0.0001),Yk_2init),label='orig_guess')
#ax.legend(loc=1, markerscale=0.1)
#savefig('PlotsofFittings/Dv_vs_v.eps')
##show()
#
#YkmFile.write (" %i      %.15f       %.15f        %.10f           %.10f\n" % (2,(fit6[0]*219474.63),(fit6[1]*219474.63),(fit6[2]*219474.63),(fit6[3]*219474.63)))
#
#print "This fit gives a D_e value of: ",(fit6[0]*219474.63)
#print "and a beta_e value of: ",(fit6[1]*219474.63)
#print "The fit is saved to Dv_vs_v.eps in the PlotsofFittings folder, and the coefficients and calculated constants printed in YkmDunhamCoeffs"
#print "Note: The very small values of D_v make the it difficult to get an accurate fit, therefore the D_e values from the minimum of the"
#print "potential curve (printed at the bottom of PUREVIB) are the values usually used."
#
#YkmFile.write (" %s                 %.15f                    %.10f\n" % ('omega_e (E_h : cm-1)',(fit4[1]),(fit4[1]*219474.63)))
#YkmFile.write (" %s                 %.15f                    %.10f\n" % ('B_e (E_h : cm-1)',(fit5[0]),(fit5[0]*219474.63)))
#YkmFile.write (" %s                 %.15f                    %.10f\n" % ('alpha_e (E_h : cm-1)',(fit5[1]*(-1.0)),(fit5[1]*219474.63*(-1.0))))
#YkmFile.write (" %s                 %.15f                    %.10f\n" % ('D_e (E_h : cm-1)',(fit6[0]),(fit6[0]*219474.63)))
#YkmFile.write (" %s                 %.15f                    %.10f\n" % ('beta_e (E_h : cm-1)',(fit6[1]),(fit6[1]*219474.63)))
#
YkmFile.close()    
