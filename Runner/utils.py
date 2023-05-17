import numpy as np
from scipy import interpolate
import healpy as hp
import argparse
from collections import OrderedDict as odict
from tqdm import tqdm
import joblib
from astropy.io import fits



class RadiusMagRelation(object):

    def __init__(self, mag, r):
    
        self.interpolator = interpolate.interp1d(mag, r, fill_value="extrapolate")

    def __call__(self, x):

        return self.interpolator(x)



class MaskableObjects(object):

    def __init__(self, filename, cushion, mag = None, m1 = None, m2 =  None, r1=None, r2=None, maxrad=None, minrad=None):
        
        #Possible names of RA and DEC in the catalogs.
        #We will search for these and standardize the column naming
        self.possible_ra_names  = ['_RAJ2000', 'RA',  'RAJ2000',  'ra']
        self.possible_dec_names = ['_DEJ2000', 'DEC', 'DECJ2000', 'dec']

        
        self.filename = filename
        self.cushion = cushion
        self.mag = mag
        self.m1 = m1
        self.m2 = m2
        self.r1 = r1
        self.r2 = r2
        self.maxrad = maxrad
        self.minrad = minrad


        self.cat = self.setup_catalog()


    def setup_catalog(self):

        # Extract data
        hdu = fits.open(self.filename)
        cat = hdu[1].data

        # Homogenize RA, DEC column names
        for r, d in zip(self.possible_ra_names, self.possible_dec_names):

            if r in cat.dtype.names:
                cat.columns[r].name = 'raj2000'
                cat.columns[d].name = 'decj2000'

        if self.mag is not None:
            cat.columns[self.mag].name = 'mag'

        return cat



class MakeMask(object):

    def __init__(self, NSIDE, MaskableObjects, RadiusMagRelation):

        self.NSIDE = NSIDE
        self.obj   = MaskableObjects
        self.maskrad = RadiusMagRelation


    def process(self, verbose = False):

        bad_pixel = np.zeros(hp.nside2npix(self.NSIDE), dtype = bool)
        
        if self.obj.mag is None: 
            maskrad = self.obj.cat['radius']
            
        else:
            maskrad = self.maskrad(self.obj.cat['mag'])
            maskrad = np.clip(maskrad, self.obj.minrad, self.obj.maxrad)

        # Convert RA, DEC info to 3D position vector    
        vec = hp.ang2vec(self.obj.cat['raj2000'], self.obj.cat['decj2000'],lonlat=True)

        for i in np.arange(maskrad.size): #for each object in coarse footprint

            #check all pixels inside the avoidance radii + some cushion in high resolution
            pixint=hp.query_disc(self.NSIDE, vec[i,:], (maskrad[i] + self.obj.cushion) * np.pi/180., inclusive = False)
            bad_pixel[pixint] = True

        return bad_pixel


class SplitJoinParallel(object):
    '''
    Takes a given Runner.
    Splits the object catalog into many parts
    and then runs it separately on maps
    Joins all the maps together in the end with OR operator
    '''
    def __init__(self, Runner, njobs = -1):
        
        self.Runner = Runner
        self.obj    = self.Runner.obj
        self.njobs = njobs if njobs != -1 else joblib.externals.loky.cpu_count()
        
        self.Runner_list = self.split_run(self.Runner)
        
        
    def split_run(self, Runner):
        
        catalog   = self.obj.cat
        Nsplits   = self.njobs
        Ntotal    = len(catalog)
        Npersplit = int(np.ceil(Ntotal/Nsplits))
        
        Runner_list = []
        for i in range(Nsplits):
            
            start = i*Npersplit
            end   = (i + 1)*Npersplit
            
            New_cat = type(self.obj)(self.obj.filename, self.obj.cushion, self.obj.mag, self.obj.m1, 
                                             self.obj.m2, self.obj.r1, self.obj.r2, self.obj.maxrad, self.obj.minrad)
            New_cat.cat = New_cat.cat[start:end]
            
            #Create a new Runner for just a subset of catalog. Has same model, map size etc.
            #Force verbose to be off as we don't want outputs for each subrun of parallel process.
            New_Runner = type(Runner)(self.Runner.NSIDE, New_cat, self.Runner.maskrad)
            
            Runner_list.append(New_Runner)
        
        return Runner_list
        

    def single_run(self, Runner):        
        
        return Runner.process()
    
    
    def process(self):
        
        jobs = [joblib.delayed(self.single_run)(Runner) for Runner in self.Runner_list]
        with joblib.parallel_backend("loky"):
            outputs = joblib.Parallel(n_jobs = self.njobs, verbose=10)(jobs)

        #Start by saying all pixels are NOT bad  pixels.
        #Then go through each map and add in pixel masking
        map_final = False
        for o in outputs:
            map_final = map_final | o
        
        return map_final