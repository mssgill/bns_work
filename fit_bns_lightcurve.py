# MSSG
# 9-25-2017

# Fitting the BNS lightcurve
# By request of R.Kessler and D.Scolnic
# (And with Juan Garcia-Bellido)

import matplotlib.pyplot as plt
import sys
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.constants import h,k,c

# File with the data in it
# (From R.Kessler, removed the one row with a negative u-band flux, and the header rows)
fname="bns_fluxes_ugrizY.dat"
#fname="gw.dat"

# Read them in simply as strs
dat = np.genfromtxt(fname, dtype='str')

# From the header, these are the cols
#     MJD  FLT  FIELD FLUXCAL  FLUXCALERR PHOTFLAG PHOTPROB ZPFLUX PSF SKYSIG SKYSIG_T GAIN XPIX YPIX  NITE EXPNUM CCDNUM OBJID
     
try:
    # Dump the columns into str arrays
    # Initze all vecs
    fnumvec = dat[:,0]
    mjd = dat[:,1]
    filt = dat[:,2]
    field = dat[:,3] 
    fluxcal = dat[:,4] 
    fluxcalerr = dat[:,5] 
    photflag = dat[:,6] 
    photprob = dat[:,7]
    zpflux = dat[:,8]
    psf = dat[:,9] 
    SKYSIG= dat[:,10]
    SKYSIG_T= dat[:,11] 
    GAIN= dat[:,12]
    XPIX= dat[:,13];     YPIX= dat[:,14]    
    NITE= dat[:,15];
    EXPNUM= dat[:,16]
    CCDNUM= dat[:,17];
    OBJID= dat[:,18]
    
    print "*********************** Have read it all in"

except:
    print "############### Not working!"


############### Get the mags as a function of day
    
nmax = len(filt)
print "File has this many entries: ", nmax

### Convert the needed columns to floats
mjdnum = mjd.astype(np.float)
flux = fluxcal.astype(np.float)

######### Now convert fluxes into mags, and put into arrays
umjd=[]
umag=[]
for n in range(1, nmax):
    if filt[n] == 'u':
        print n, filt[n], mjdnum[n]
        umjd.append(mjdnum[n])
        umag.append(27.5 - 2.5*math.log10(flux[n]))

gmjd=[]
gmag=[]
for n in range(1, nmax):
    if filt[n] == 'g':
        print n, filt[n], mjdnum[n]
        gmjd.append(mjdnum[n])
        gmag.append(27.5 - 2.5*math.log10(flux[n]))

rmjd=[]
rmag=[]
for n in range(1, nmax):
    if filt[n] == 'r':
        print n, filt[n], mjdnum[n]
        rmjd.append(mjdnum[n])
        rmag.append(27.5 - 2.5*math.log10(flux[n]))
        
imjd=[]
imag=[]
for n in range(1, nmax):
    if filt[n] == 'i':
        print n, filt[n], mjdnum[n]
        imjd.append(mjdnum[n])
        imag.append(27.5 - 2.5*math.log10(flux[n]))
        
zmjd=[]
zmag=[]
for n in range(1, nmax):
    if filt[n] == 'z':
        print n, filt[n], mjdnum[n]
        zmjd.append(mjdnum[n])
        zmag.append(27.5 - 2.5*math.log10(flux[n]))

ymjd=[]
ymag=[]
for n in range(1, nmax):
    if filt[n] == 'Y':
        print n, filt[n], mjdnum[n]
        ymjd.append(mjdnum[n])
        ymag.append(27.5 - 2.5*math.log10(flux[n]))


######################################## Plot mags as funct of day, and fit to poly
        
fig = plt.figure()
ax = fig.add_subplot(111)

u = ax.scatter(umjd,umag, color = 'b')
g = ax.scatter(gmjd, gmag, color = 'g')
r = ax.scatter(rmjd, rmag, color = 'r')
i = ax.scatter(imjd, imag, color = 'crimson')
z = ax.scatter(zmjd, zmag, color = 'fuchsia')
y = ax.scatter(ymjd, ymag, color = 'gold')


# calculate polynomial
ufit = np.polyfit(umjd, umag, 3)
gfit = np.polyfit(gmjd, gmag, 2)
rfit = np.polyfit(rmjd, rmag, 4)
ifit = np.polyfit(imjd, imag, 10)
zfit = np.polyfit(zmjd, zmag, 2)
yfit = np.polyfit(ymjd, ymag, 2)

uff = np.poly1d(ufit)
gff = np.poly1d(gfit)
rff = np.poly1d(rfit)
iff = np.poly1d(ifit)
zff = np.poly1d(zfit)
yff = np.poly1d(yfit)


# calculate new x's and y's
utimefit = np.linspace(umjd[0], umjd[-1], 50) ; umagfit = uff(utimefit) ; ax.plot(utimefit,umagfit)
gtimefit = np.linspace(gmjd[0], gmjd[-1], 50) ; gmagfit = gff(gtimefit) ; ax.plot(gtimefit,gmagfit)
rtimefit = np.linspace(rmjd[0], rmjd[-1], 50) ; rmagfit = rff(rtimefit) ; ax.plot(rtimefit,rmagfit)
itimefit = np.linspace(imjd[0], imjd[-1], 50) ; imagfit = iff(itimefit) ; ax.plot(itimefit,imagfit)
ztimefit = np.linspace(zmjd[0], zmjd[-1], 50) ; zmagfit = zff(ztimefit) ; ax.plot(ztimefit,zmagfit)
ytimefit = np.linspace(ymjd[0], ymjd[-1], 50) ; ymagfit = yff(ytimefit) ; ax.plot(ytimefit,ymagfit)

plt.legend((u,g,r,i,z,y),('u','g','r','i','z','Y'))
plt.show()



################################# Fit BB to the mags on a given day  (not fully working yet -- 9-25-2017)
# Using code from: http://python4esac.github.io/fitting/example_blackbody.html


def blackbody_lam(lam, T):
    """ Blackbody as a function of wavelength (um) and temperature (K).
    returns units of erg/s/cm^2/cm/Steradian
    """
    from scipy.constants import h,k,c   # h = 6.626e-34 ; c = 3.0e+8 ; k = 1.38e-23

    lam = 1e-6 * lam # convert to metres
    return 2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))

wa = np.linspace(0.1, 2, 100)   # wavelengths in um

# From: http://www.sdss.org/instruments/camera/#Filters
umid=3551; gmid = 4686 ; rmid = 6166 ; imid =7480 ; zmid = 8932  # In Angstroms, clearly
cf = 1e-4
xdat = [ umid*cf , gmid*cf , rmid*cf , imid*cf, zmid*cf]  # Now in um

########## 
ndays=1
for n in range(0, ndays):
    ydat = [ umag[n], gmag[n], rmag[n], imag[n], zmag[n] ]
    plt.scatter(xdat,ydat)
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Intensity (erg/s/cm$^2$/cm/Steradian)')
    plt.show()

# plot the input model and synthetic data

#plt.figure()
#plt.plot(wa, y1, ':', lw=2, label='T1=%.0f' % T1)
#plt.plot(wa, ydata, ls='steps-mid', lw=2, label='Fake data')
#plt.xlabel('Wavelength (microns)')
#plt.ylabel('Intensity (erg/s/cm$^2$/cm/Steradian)')

# Note the initial guess values for T1 (p0 keyword below). They
# are quite different to the known true values, but not *too*
# different. If these are too far away from the solution curve_fit()
# will not be able to find a solution. This is not a Python-specific
# problem, it is true for almost every fitting algorithm for
# non-linear models. The initial guess is important!

def func(wa, T1):
    return blackbody_lam(wa, T1)

popt, pcov = curve_fit(func, wa, ydat, p0=10000) # , sigma=sigma)

# get the best fitting parameter values and their 1 sigma errors
# (assuming the parameters aren't strongly correlated).


bestT1 = popt
sigmaT1 = np.sqrt(np.diag(pcov))

ybest = blackbody_lam(wa, bestT1) 

print 'True model values'
print '  T1 = %.2f' % T1

print 'Parameters of best-fitting model:'
print '  T1 = %.2f +/- %.2f' % (bestT1, sigmaT1)

degrees_of_freedom = len(wa) - 2
resid = (ydat - func(wa, *popt)) / sigmaT1
chisq = np.dot(resid, resid)

print degrees_of_freedom, 'dof'
print 'chi squared %.2f' % chisq
print 'nchi2 %.2f' % (chisq / degrees_of_freedom)

# plot the solution

plt.plot(wa, ybest, label='Best fitting\nmodel')
plt.legend(frameon=False)
plt.savefig('fit_bb.png')
plt.show()

sys.exit()
