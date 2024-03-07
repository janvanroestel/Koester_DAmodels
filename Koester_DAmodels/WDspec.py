import glob
import os
import pkg_resources
import copy
import numpy as np
import scipy.ndimage.filters as filters
import extinction



# load a spectrum
def _loadspec(f,wmax=50000,Nw=50):
    ''' load a individual spectrum file. Extrapolate the IR wavelengths in loglog
    '''
    spec = np.loadtxt(f) # load wl grid

    # extrapolate WD

    if wmax > 30000. and Nw is not None:
        # fit from 25000 and higher
        m = (spec[:,0]>25000)
        z = np.polyfit(np.log10(spec[m,0]),np.log10(spec[m,1]),1.)
        p = np.poly1d(z)

        # sample the extrapolated function
        ewl = np.linspace(30000,wmax,Nw)
        es = 10**p(np.log10(ewl))

        # add the extrapolated part to the spectrum
        spec = np.r_[spec,np.c_[ewl,es]]

    return spec



def _bilinear_interpolation(x, y, X, Y):
    '''Interpolate (x,y) in a rectangular grid defined by X and Y

    input: 
        x : float, x position
        y : float, y position
        X : 1D-array of x-grid values
        Y : 1D-array of y-grid values

    output:
        tuple of corner points, x,y,w

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    # find corners
    x_1 = X[np.searchsorted(X,x,side='right')-1] 
    x_2 = X[np.searchsorted(X,x,side='right')] 
    y_1 = Y[np.searchsorted(Y,y,side='right')-1] 
    y_2 = Y[np.searchsorted(Y,y,side='right')] 

    # calculate weights
    dx1 = (x-x_1)/(x_2-x_1)
    dx2 = 1-dx1
    dy1 = (y-y_1)/(y_2-y_1)
    dy2 = 1-dy1

    w_a = dx2*dy2 # lower left
    w_b = dx1*dy2 # lower right
    w_c = dx2*dy1 # upper left
    w_d = dx1*dy1 # upper right 

    # return the corner values, x,y and weight
    return ([x_1,y_1,w_a],
            [x_2,y_1,w_b],
            [x_1,y_2,w_c],
            [x_2,y_2,w_d])



class WDmodels:
    ''' The class that contains the WD spectra models.
    When used to create an instance, it will load the model spectra. 
    The models will be generated on a given grid, from an min to max WL range
    and with a smoothing in A pre-applied (because it is slow)
    '''
    
    def __init__(self,ddir=None,
            dw=None,smooth=0,wmin=None,wmax=None):
        """ Initiate the class with the data. dw is the wavelength stepsize and
        smooth is the sigma of the Gaussian smoothing kernel
        

        input:
            ddir    : str, the directory which contains the KoesterDA models. You 
                      can download them here
                      http://svo2.cab.inta-csic.es/theory/newov/index.php
                      These should be provided with the install
            dw      : float, wavelength step. Use None to keep the orginal grid
            smooth  : float, Gaussian sigma for convolution kernel. Check your 
                      spectragraph+grating for what this value should be 
                      (be carefull to convert FWHM to sigma)
            wmin    : float, the mimimum wavelength for the model
            wmax    : float, the maximum wavelength for the model

        """

        # get the datafiles
        if ddir is None:
            ddir = pkg_resources.resource_filename('Koester_DAmodels', 'models_koester2_1571804487/')
        print(ddir)
        files = glob.glob(ddir+'*.dat.txt')
        files.sort()

        # get T and logg values from the filenames!!!
        T = np.array([float(os.path.basename(f)[2:7]) for f in files])
        logg = np.array([float(os.path.basename(f)[8:11]) for f in files])*0.01

        print('Loading %d model data...' %len(files))
        modelspectra = [_loadspec(f) for f in files]
        print('Spectra loaded')
        
        # calculate statistics
        N = np.array([np.size(d[:,0]) for d in modelspectra])
        
        # get wavelength range
        w = copy.deepcopy(modelspectra[0][:,0])

        # set or get wavelenght range (can make things more efficient)
        if wmin is None:
            #print('changing wmin')
            wmin = np.min(w)
        if wmax is None:
            #print('changing wmax')
            wmax = np.max(w)
        if dw is not None:
            # make a regular grid
            w = np.arange(wmin,wmax,dw)
        else:
            w = w[(w>wmin)&(w<wmax)]

        # smooth models, this is needed if you are fitting an observed spectrum 
        # at a certain resolution
        if smooth>0:
            wf = np.arange(np.min(w),np.max(w),smooth/5.)
            k = smooth*5 # smoothing kernel for regular grid
            s = lambda y: filters.gaussian_filter(y,k)
            l = lambda d: np.interp(w,wf,s(np.interp(wf,d[:,0],d[:,1])))
        else:
            l = lambda d: np.interp(w,d[:,0],d[:,1])

        # store in dict for easy lookup
        print('Smoothing and interpolating...')
        modeldict = dict()
        for _T,_logg,m in zip(T,logg,modelspectra):
            modeldict[(_T,_logg)] = l(m)
        print('Done')

        # save the results to the Instance
        self.w = w
        self.spectra = modeldict
        self.T_grid = np.sort(np.unique(T))
        self.logg_grid = np.sort(np.unique(logg))



    # for a given T and logg, return a spectrum
    def get_spectrum(self,T,logg,RV=0,kernel=0,R=None,dist=None,EBV=None,Fnu=False):
        """ generate a spectrum by interpolating the model
        input:
            T    : float, temperature of the WD in K
            logg : float, surface gravity in log[cgs]
            RV   : float, RV shift in km/s 
            R    : float, radius in solar units, rescales the spectrum
            dist : float, distance in pc, rescales the spectrum
            EBV  : float, extinction in B-V 
            Fnu  : bool, convert to Fnu

        output:
            s    : 1D-array, the flux in 4pi*Eddington flux in erg/cm2/s/A
                   use self.w to get the wavelength grid

        """

        # set constants
        pc = 3.086e+16 # parsec in meters     
        Rsun = 695700000.0 # solar radius in meters

        # check bounds
        if T>np.max(self.T_grid) or T<np.min(self.T_grid):
            raise ValueError('%dK is out of bounds' %T)
        if logg>np.max(self.logg_grid) or logg<np.min(self.logg_grid):
            raise ValueError('%g is out of bounds' %logg)

        # do bilinear interpolation in T logg grid
        corners = _bilinear_interpolation(T,logg,self.T_grid,self.logg_grid)
        s = np.array([c[2]*self.spectra[(c[0],c[1])] for c in corners])
        s = np.sum(s,axis=0)

        # set wavelength
        wl = self.w

        # apply velocity shift
        if RV != 0:
            wn = wl * (1.+RV/(3.*10**5))
            s = np.interp(wl,wn,s) # resample to wl grid
        
        # convert from Fl to Fnu
        if Fnu:
            s *= 3.33564095E+04*wl**2

        # radius scaling
        if R is not None:
            s *= 4*np.pi*(R*Rsun)**2

        # distance scaling
        if dist is not None:
            s /= 4*np.pi*(dist*pc)**2

        # extinction
        if EBV is not None:
            A_V = EBV*3.1
            ext = extinction.fitzpatrick99(wl, A_V, 3.1) # wavelength in AA
            s *= 10**(-0.4 * ext)

        return np.c_[wl,s]

