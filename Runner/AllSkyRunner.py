import numpy as np
from scipy import interpolate
import healpy as hp
import argparse
from collections import OrderedDict as odict
from tqdm import tqdm
from .utils import RadiusMagRelation, MaskableObjects, MakeMask, SplitJoinParallel
import os


class AllSkyRunner(object):

    def __init__(self, NSIDE = 4096):

        self.NSIDE = NSIDE
        self.bits  = self.setup_data()


    def setup_data(self):

        file_path = os.path.dirname(__file__) + '/../radii_mag_relation/Y6A2/'

        gaia_bright_x      = np.loadtxt(path + 'Gaia_bright_magnitudes_mid_bin.txt')
        gaia_bright_y      = np.loadtxt(path + 'Gaia_bright_conservative_radii_cuts.txt')
        gaia_bright_interp = RadiusMagRelation(gaia_bright_x, gaia_bright_y)

        gaia_mid_x         = np.loadtxt(path + 'Gaia_mid_magnitudes_mid_bin.txt')
        gaia_mid_y         = np.loadtxt(path + 'Gaia_mid_conservative_radii_cuts.txt')
        gaia_mid_interp    = RadiusMagRelation(gaia_mid_x, gaia_mid_y)

        yale_x             = np.loadtxt(path + 'yale_magnitudes_mid_bin.txt')
        yale_y             = np.loadtxt(path + 'yale_conservative_radii_cuts.txt')
        yale_interp        = RadiusMagRelation(yale_x, yale_y)

        mass_bright_x      = np.loadtxt(path + '2mass_bright_magnitudes_mid_bin.txt')
        mass_bright_y      = np.loadtxt(path + '2mass_bright_conservative_radii_cuts.txt')
        mass_bright_interp = RadiusMagRelation(mass_bright_x, mass_bright_y)

        mass_faint_x       = np.loadtxt(path + '2mass_faint_magnitudes_mid_bin.txt')
        mass_faint_y       = np.loadtxt(path + '2mass_faint_conservative_radii_cuts.txt')
        mass_faint_interp  = RadiusMagRelation(mass_faint_x, mass_faint_y)

        leda_x             = np.loadtxt(path + 'leda_magnitudes_mid_bin.txt')
        leda_y             = np.loadtxt(path + 'leda_conservative_radii_cuts.txt')
        leda_interp        = RadiusMagRelation(leda_x, leda_y)
        
        FOREGROUND_BITS = odict([
		(256, dict(name = 'Famous stars',
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/famous_stars/list_famous_stars_in_DES.fits',
                   interp   = None,
                   mag = None, m1 = None, m2 = None, r1 = None, r2 = None, maxrad = None, minrad = None, cushion = 30./3600.)),
                   
		# (128, dict(name = "Near the LMC",
        #            filename = None,
        #            mag = None, m1 = None, m2 = None, r1 = None, r2 = None, maxrad = None, minrad = None, cushion = None)),

		(64,  dict(name = 'Globular clusters',
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/harris_globular_cluster_cat/harris_globclust_cat_in_DES.fit',
                   interp   = None,
                   mag = None, m1 = None, m2 = None, r1 = None, r2 = None, maxrad = None, minrad = None, cushion = 5./3600.)),

		(32,  dict(name = "Large nearby galaxy (HyperLEDA catalog)",
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/hyperleda_db/hyperleda_in_des.fits',
                   interp   = leda_interp,
                   mag = 'BT', m1 = min(leda_x), m2 = max(leda_x), r1 = 0.1, r2 = 0.03, maxrad = max(leda_y), minrad = min(leda_y), cushion = .5/3600.)),

		(16,  dict(name = "2MASS fainter star region (8<J<12)",
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/2mass_stars/twomass_in_des_faint.fits',
                   interp   = mass_faint_interp,
                   mag = 'J', m1 = min(mass_faint_x), m2 = max(mass_faint_x), r1 = 0.03, r2 = 0.008, 
                   maxrad = max(mass_faint_y), minrad = min(mass_faint_y), cushion = 5./3600.)),

		(8,   dict(name = "Gaia mid stars (7.<G<11.5)",
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_mid.fits',
                   interp   = gaia_mid_interp,
                   mag = 'PHOT_G_MEAN_MAG', m1 = min(gaia_mid_x), m2 = max(gaia_mid_x), r1 = 0.03, r2 = 0.008, 
                   maxrad = max(gaia_mid_y), minrad = min(gaia_mid_y), cushion = 5./3600.)),

		(4,   dict(name = "2MASS bright star region (J<8.)",
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/2mass_stars/twomass_in_des_bright.fits',
                   interp   = mass_bright_interp,
                   mag = 'J', m1 = min(mass_bright_x), m2 = max(mass_bright_x), r1 = 0.075, r2 = 0.03, 
                   maxrad = max(mass_bright_y), minrad = min(mass_bright_y), cushion = 5./3600.)),

		(2,   dict(name = "Yale bright stars",
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/yale_stars/yale5th_in_DES.fits',
                   interp   = yale_interp,
                   mag = 'Vmag', m1 = min(yale_x), m2 = max(yale_x), r1 = 0.4, r2 = 0.15, 
                   maxrad = max(yale_y), minrad = min(yale_y), cushion = .5/3600.)),

		(1,   dict(name = "Gaia bright stars (G<7.)",
                   filename = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_bright.fits',
                   interp   = gaia_bright_interp,
                   mag = 'PHOT_G_MEAN_MAG', m1 = min(gaia_bright_x), m2 = max(gaia_bright_x), r1 = 0.03, r2 = 0.008, 
                   maxrad = max(gaia_bright_y), minrad = min(gaia_bright_y), cushion = .5/3600.))
	    ])

        return FOREGROUND_BITS


    def process(self):

        map = np.zeros(hp.nside2npix(self.NSIDE), dtype = np.uint16)

        for bit, data in self.bits:

            name    = data.pop('name')
            maskrad = data.pop('interp')
            Objects = MaskableObjects(**data)
            Runner  = MakeMask(NSIDE = 4096, MaskableObjects = Objects, RadiusMagRelation = maskrad)
            Runner  = SplitJoinParallel(Runner, njobs = 10)
            
            bad_pix = Runner.process()
            bad_pix = np.where(bax_pix, bit, 0)
            map     = map | bad_pix

        
        return map


if __name__ == '__main__':


    X = AllSkyRunner()
    m = AllSkyRunner.process()

    hp.write_map('/scratch/midway2/dhayaa/Foreground_Mask_Allsky.fits', m, overwrite = True)

        
    
        