from collections import OrderedDict as odict
#from astropy import coordinates as coo
from scipy import interpolate
import numpy as np

def init():
    global FOREGROUND_BITS
    global BAD_BITS
    global gaia_bright_interp,gaia_mid_interp,gaia_faint_interp,yale_interp,mass_bright_interp,mass_faint_interp, leda_interp
    
    gaia_bright_x = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/gaia_stars/bright/Gaia_bright_magnitudes_mid_bin.txt')
    gaia_bright_y = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/gaia_stars/bright/Gaia_bright_conservative_radii_cuts.txt')
    gaia_bright_interp = interpolate.interp1d(gaia_bright_x,gaia_bright_y,fill_value="extrapolate")

    gaia_mid_x = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/gaia_stars/mid/Gaia_mid_magnitudes_mid_bin.txt')
    gaia_mid_y = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/gaia_stars/mid/Gaia_mid_conservative_radii_cuts.txt')
    gaia_mid_interp = interpolate.interp1d(gaia_mid_x,gaia_mid_y,fill_value="extrapolate")

    yale_x = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/yale_stars/yale_magnitudes_mid_bin.txt')
    yale_y = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/yale_stars/yale_conservative_radii_cuts.txt')
    yale_interp = interpolate.interp1d(yale_x,yale_y,fill_value="extrapolate")

    mass_bright_x = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/2mass_stars/bright/2mass_bright_magnitudes_mid_bin.txt')
    mass_bright_y = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/2mass_stars/bright/2mass_bright_conservative_radii_cuts.txt')
    mass_bright_interp = interpolate.interp1d(mass_bright_x,mass_bright_y,fill_value="extrapolate")

    mass_faint_x = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/2mass_stars/faint/2mass_faint_magnitudes_mid_bin.txt')
    mass_faint_y = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/2mass_stars/faint/2mass_faint_conservative_radii_cuts.txt')
    mass_faint_interp = interpolate.interp1d(mass_faint_x,mass_faint_y,fill_value="extrapolate")

    leda_x = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/hyperleda_db/leda_magnitudes_mid_bin.txt')
    leda_y = np.loadtxt('/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/hyperleda_db/leda_conservative_radii_cuts.txt')
    leda_interp = interpolate.interp1d(leda_x,leda_y,fill_value="extrapolate")
    
    FOREGROUND_BITS = odict([
		(512, dict(name='Total area')),
        (256, dict(name='Famous stars',filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/famous_stars/list_famous_stars_in_DES.fits',mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=30./3600.)),
#		(64, dict(name='Famous stars',filename='famous_stars_footprint.fits',mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=30./3600.)),
		(128,  dict(name="Near the LMC",filename=None,mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=None)),
		(64, dict(name='Globular clusters',filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/harris_globular_cluster_cat/harris_globclust_cat_in_DES.fit',mag=None,m1=None,m2=None,r1=None,r2=None,maxrad=None,minrad=None,cushion=5./3600.)),
		(32,   dict(name="Large nearby galaxy (HyperLEDA catalog)",filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/hyperleda_db/hyperleda_in_des.fits',mag='BT',m1=min(leda_x),m2=max(leda_x),r1=0.1,r2=0.03,maxrad=max(leda_y),minrad=min(leda_y),cushion=.5/3600.)),
		(16,   dict(name="2MASS fainter star region (8<J<12)",filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/2mass_stars/twomass_in_des_faint.fits',mag='J',m1=min(mass_faint_x),m2=max(mass_faint_x),r1=0.03,r2=0.008,maxrad=max(mass_faint_y),minrad=min(mass_faint_y),cushion=5./3600.)),	
		(8,   dict(name="Gaia mid stars (7.<G<11.5)",filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_mid.fits',mag='PHOT_G_MEAN_MAG',m1=min(gaia_mid_x),m2=max(gaia_mid_x),r1=0.03,r2=0.008,maxrad=max(gaia_mid_y),minrad=min(gaia_mid_y),cushion=5./3600.)),
		(4,   dict(name="2MASS bright star region (J<8.)",filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/2mass_stars/twomass_in_des_bright.fits',mag='J',m1=min(mass_bright_x),m2=max(mass_bright_x),r1=0.075,r2=0.03,maxrad=max(mass_bright_y),minrad=min(mass_bright_y),cushion=5./3600.)),
		(2,  dict(name="Yale bright stars",filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/yale_stars/yale5th_in_DES.fits',mag='Vmag',m1=min(yale_x),m2=max(yale_x),r1=0.4,r2=0.15,maxrad=max(yale_y),minrad=min(yale_y),cushion=.5/3600.)),
		(1,   dict(name="Gaia bright stars (G<7.)",filename='/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_bright.fits',mag='PHOT_G_MEAN_MAG',m1=min(gaia_bright_x),m2=max(gaia_bright_x),r1=0.03,r2=0.008,maxrad=max(gaia_bright_y),minrad=min(gaia_bright_y),cushion=.5/3600.))
	])

	
