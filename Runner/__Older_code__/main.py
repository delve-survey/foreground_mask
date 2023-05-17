
import healpy as hp
import numpy as np
from astropy.io import fits
import sys
import os
import matplotlib.pyplot as plt
import smatch
# from smatch import Catalog
import esutil
from optparse import OptionParser
import fitsio
import time

def calc_radii(goodobjfile, badobjfile, extfile, maxrad, nside, magname, catname, ra_name, dec_name, nbin_sample, path_aux, automated_radius, radii_to_save_opt, radii_to_save_con, magnitudes_array):
  """ Load external catalogue """
  hdu = fits.open(extfile)
  extcat = hdu[1].data
  # Unify RA,DEC naming
  extcat.columns[ra_name].name = 'RAJ2000'
  extcat.columns[dec_name].name = 'DEJ2000'
  """ We need to go back to RA in [0, 360] """
  extcat['RAJ2000'] = np.where(extcat['RAJ2000'] < 0, extcat['RAJ2000'] + 360, extcat['RAJ2000'])
  sort = np.argsort(extcat[magname])
  extcat = extcat[sort]
  """ Definig nbins for each sample, equal number of objects per bin """
  nbin = nbin_sample 
  nperbin=np.ceil(extcat.size/float(nbin)).astype(np.int32)
  """ Loading Good & Bad object files """
  hdu = fits.open(goodobjfile)
  goodobjs=hdu[1].data
  hdu = fits.open(badobjfile)
  badobjs=hdu[1].data
 
  cmap = plt.get_cmap('rainbow')
  # wrapping RA,DEJ2000
  goodrawrap=goodobjs['RA'].copy()
  hi,=np.where(goodrawrap > 180.0)
  if (hi.size > 0) : goodrawrap[hi] = goodrawrap[hi] - 360.0
  badrawrap = badobjs['RA'].copy()
  hi,=np.where(badrawrap > 180.0)
  if (hi.size > 0) : badrawrap[hi] = badrawrap[hi] - 360.0
  extcatrawrap=extcat['RAJ2000'].copy()
  hi,=np.where(extcatrawrap > 180.0)
  if (hi.size > 0) : extcatrawrap[hi] = extcatrawrap[hi] - 360.0

  for i in np.arange(nbin):
    print("working on bin: %d" % (i))
    use=np.arange(nperbin,dtype=np.int32)+i*nperbin
    good,=np.where(use < extcat.size)
    use=use[good]
    cat = smatch.Catalog(extcat['RAJ2000'][use],extcat['DEJ2000'][use],maxrad,nside=nside)
    print(cat)
    cat.match(goodobjs['RA'],goodobjs['DEC'],maxmatch=0)
    good_matches = cat.matches
    print(good_matches)
    dra=(goodrawrap[good_matches['i2']] - extcatrawrap[use[good_matches['i1']]])*np.cos(extcat['DEJ2000'][use[good_matches['i1']]]*np.pi/180.)
    ddec=(goodobjs['DEC'][good_matches['i2']] - extcat['DEJ2000'][use[good_matches['i1']]]) 
    try:
    	h=esutil.stat.histogram(np.arccos(good_matches['cosdist'])*180.0/np.pi,binsize=maxrad/50,min=0,max=maxrad,more=True)
    except:
        continue

    hcent = h['center']

    hdens = h['hist']/(np.pi*(h['high']**2-h['low']**2.))

    hcent=hcent[0:hcent.size-1]

    hdens=hdens[0:hdens.size-1]

    hdens=hdens/hdens[hdens.size-1]


    cat.match(badobjs['RA'],badobjs['DEC'],maxmatch=0)

    bad_matches = cat.matches


    bdra=(badrawrap[bad_matches['i2']] - extcatrawrap[use[bad_matches['i1']]])*np.cos(extcat['DEJ2000'][use[bad_matches['i1']]]*np.pi/180.)

    bddec=(badobjs['DEC'][bad_matches['i2']] - extcat['DEJ2000'][use[bad_matches['i1']]]) 

    bh=esutil.stat.histogram(np.arccos(bad_matches['cosdist'])*180.0/np.pi,binsize=maxrad/50,min=0,max=maxrad,more=True)

    bhcent = bh['center']

    bhdens = bh['hist']/(np.pi*(bh['high']**2-bh['low']**2.))

    bhcent=bhcent[0:bhcent.size-1]

    bhdens=bhdens[0:bhdens.size-1]

    bhdens=bhdens/bhdens[bhdens.size-1]

    """ 
    Defining the radius to cut: radius at which the 95% of bad object density is collected OR
    when the relative density of bad objects is under 2
    
    """
    if automated_radius:
        max_arg = np.argmax(bhdens)
#        ind = np.argmin(np.abs(np.cumsum(bhcent*bhdens/np.sum(bhcent*bhdens)) -0.85*np.ones_like(bhcent)))
        ind = np.argmin(np.abs(bhdens[max_arg:] - 2*np.ones_like(bhdens[max_arg:])))
        radius_con = bhcent[max_arg+ind]
        radius_opt = bhcent[max_arg]
	
        radii_to_save_con = np.append(radii_to_save_con, radius_con)
        radii_to_save_opt = np.append(radii_to_save_opt, radius_opt)
     

    """
    Plots
    """

    fig=plt.figure(1,figsize=(17,5))

    fig.clf()

  

    ax=fig.add_subplot(131)

    ax.set_xlabel('delta RA (deg)')

    ax.set_ylabel('delta DEC (deg)')

    ax.set_aspect('equal')

    ax.set_xlim(-maxrad,maxrad)

    ax.set_ylim(-maxrad,maxrad)

    ax.set_title('%.2f < %s < %.2f, Gal Density' % (extcat[magname][use].min(),magname,extcat[magname][use].max()))


    ax.hexbin(dra,ddec,gridsize=200,cmap=cmap,bins='log')

    for div in np.arange(1,6):
        
      circle=plt.Circle((0,0),np.round(maxrad -div*maxrad/5,2),color='k',fill=False)

      ax.add_patch(circle)
    
    
    if automated_radius:
        
        circle=plt.Circle((0,0),radius_con, color='r',fill=False)

        ax.add_patch(circle)
        
        circle=plt.Circle((0,0),radius_opt, color='r',fill=False)

        ax.add_patch(circle)
        
#    else: 
#    
#      """ Optimistic and conservative estimations """
#
#      circle=plt.Circle((0,0),radius_con[i],color='r',fill=False)
#
#      ax.add_patch(circle)
#    
#      circle=plt.Circle((0,0),radius_opt[i],color='m',fill=False)
#
#      ax.add_patch(circle)

  

    ax=fig.add_subplot(132)

    ax.set_xlabel('delta RA (deg)')

    ax.set_ylabel('delta DEC (deg)')

    ax.set_aspect('equal')

    ax.set_xlim(-maxrad,maxrad)

    ax.set_ylim(-maxrad,maxrad)

    ax.set_title('%.2f < %s < %.2f, Bad Objects' % (extcat[magname][use].min(),magname,extcat[magname][use].max()))


    ax.hexbin(bdra,bddec,gridsize=200,cmap=cmap,bins='log')

    for div in np.arange(1,6):
        
      circle=plt.Circle((0,0),np.round(maxrad -div*maxrad/5,2),color='k',fill=False)

      ax.add_patch(circle)
      
    if automated_radius:
        
       circle=plt.Circle((0,0),radius_con, color='r',fill=False)

       ax.add_patch(circle)
       
       circle=plt.Circle((0,0),radius_opt, color='r',fill=False)

       ax.add_patch(circle)
    
#    else: 
#      """ Optimistic and conservative estimations """
#
#      circle=plt.Circle((0,0),radius_con[i],color='r',fill=False)
#
#      ax.add_patch(circle)
#    
#      circle=plt.Circle((0,0),radius_opt[i],color='m',fill=False)
#
#      ax.add_patch(circle)


    # and the radial plots...

    ax=fig.add_subplot(133)

    ax.set_xlabel('Distance (deg)')

    ax.set_ylabel('Relative Density')

    ax.set_xlim(0,maxrad)

    #ax.set_ylim(0,2)
    ax.set_yscale('log')

    ax.set_title('%.2f < %s < %.2f' % (extcat[magname][use].min(),magname,extcat[magname][use].max()))

    ax.grid()

    ax.plot(hcent,hdens,'b-',label='Galaxy Density')

    ax.plot(bhcent,bhdens,'y--',label='Bad Object Density')

    ax.plot(np.array([0,maxrad]),np.array([1,1]),'k:')
    
    if automated_radius:
        
      ax.axvline(radius_con, linestyle='--',c='r')
      ax.axvline(radius_opt, linestyle='--',c='r')
      
#    else:
        
#     ax.axvline(radius_con[i],linestyle='--',c='r')
#     ax.axvline(radius_opt[i],linestyle='--',c='m')


    ax.legend()
    if not os.path.exists(path_aux):
        # Si no existe, crea el directorio
        os.makedirs(path_aux)
    else:
        print(f"{path_aux} already exists")
    fig.savefig(os.path.join(path_aux,'%s_%.2f_%.2f.png' % (catname, extcat[magname][use].min(),extcat[magname][use].max())))
    #plt.show()
    
    magnitudes_array = np.append(magnitudes_array, np.mean(np.array([extcat[magname][use].min(),extcat[magname][use].max()])))

    if i == nbin - 1:
      np.savetxt(os.path.join(path_aux, '%s_conservative_radii_cuts.txt') % (catname), radii_to_save_con)
      np.savetxt(os.path.join(path_aux, '%s_optimistic_radii_cuts.txt') % (catname), radii_to_save_opt)
      np.savetxt(os.path.join(path_aux, '%s_magnitudes_mid_bin.txt') % (catname), magnitudes_array)
      
def main():

  '''

  Run code with options

  '''
  # Read external options
  usage = "%prog [options]"
  parser = OptionParser(usage=usage)
  parser.add_option("--gaia_bright",action="store_true",dest="toggle_gaia_bright",help="Compute radii for GAIA: bright sample stars",default=False)
  parser.add_option("--gaia_mid",action="store_true",dest="toggle_gaia_mid",help="Compute radii for GAIA: mid sample stars",default=False)
  parser.add_option("--yale",action="store_true",dest="toggle_yale",help="Compute radii for Yale stars",default=False)
  parser.add_option("--2mass_bright",action="store_true",dest="toggle_2mb",help="Compute radii for bright 2MASS stars",default=False)
  parser.add_option("--2mass_faint",action="store_true",dest="toggle_2mf",help="Compute radii for faint 2MASS stars",default=False)
  parser.add_option("--leda",action="store_true",dest="toggle_leda",help="Compute radii for HyperLEDA galaxies",default=False)
  """ Good an Bad objects files """
  goodobjfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/good-objects/y6a2_good_objects.fits'
  badobjfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/bad-objects/y6a2_bad_objects.fits'
  """ QUESTION: Compute with automated radius? """
  # automated_radius = True / False
  automated_radius=True
  """ Initialize the array where saving radii cuts """
  radii_to_save_opt = np.array([])
  radii_to_save_con = np.array([])
  magnitudes_array = np.array([])

  (options, args) = parser.parse_args() 
  # Nside used previous to the high resolution update (moving to Healsparse)
  # nside = 512 #for matches
  nside = 1024
  path_to_save = f'/scratch/davsan06/Doctorado/foreground_mask_Y6A2/results/nside_{nside}'
  print(f'Working with nside={nside}')
  if options.toggle_gaia_bright:
    # Point to GAIA bright sample file  
    extfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_bright.fits'
    magname = 'PHOT_G_MEAN_MAG'
    catname = 'Gaia_bright'
    maxrad = 0.25
    ra_name = 'RA'
    dec_name = 'DEC'
    nbin_sample = 20
#    """ Optimistic and conservative radius """
#    radius_opt = [0.05, 0.06, 0.05, 0.05, 0.07, 0.05, 0.05, 0.025, 0.03, 0.03, 0.03, 0.03, 0.025, 0.025, 0.03, 0.03, 0.03, 0.025, 0.025]
#    radius_con = [0.15, 0.10, 0.10, 0.08, 0.09, 0.09, 0.08, 0.06, 0.06, 0.06, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03, 0.025, 0.025]
    path_aux = os.path.join(path_to_save, "gaia_stars", "bright")
  elif options.toggle_gaia_mid:
    # Point to GAIA mid bright sample
    extfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_mid.fits'
    magname = 'PHOT_G_MEAN_MAG'
    catname = 'Gaia_mid'
    maxrad = 0.04
    ra_name = 'RA'
    dec_name = 'DEC'
    nbin_sample = 20
#    """ Optimistic and conservative radius """
#    radius_opt = [0.021, 0.025, 0.018, 0.018, 0.017, 0.015, 0.015, 0.014, 0.013, 0.012, 0.011, 0.011, 0.011, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010]
#    radius_con = [0.030, 0.030, 0.025, 0.020, 0.020, 0.017, 0.016, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.012]
    path_aux = os.path.join(path_to_save, "gaia_stars", "mid")
  elif options.toggle_yale:
    # Point to Yale stars sample
    extfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/yale_stars/yale5th_in_DES.fits'
    magname = 'Vmag'
    catname = 'yale'
    maxrad = 0.3
    ra_name = 'RAJ2000'
    dec_name = 'DEJ2000'
    nbin_sample = 20
#    """ Optimistic and conservative radius """
#    radius_opt = [0.090, 0.080, 0.070, 0.050, 0.050, 0.050, 0.050, 0.045, 0.045, 0.045, 0.040, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035, 0.030]
#    radius_con = [0.150, 0.150, 0.150, 0.100, 0.100, 0.080, 0.080, 0.070, 0.070, 0.070, 0.065, 0.060, 0.050, 0.050, 0.050, 0.050, 0.050, 0.045, 0.050, 0.040, 0.040]
    path_aux = os.path.join(path_to_save, "yale_stars")
  elif options.toggle_leda:
    # Point to HyperLeda sample
    extfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/hyperleda_db/hyperleda_in_des.fits'
    magname = 'BT'
    catname = 'leda'
    maxrad = 0.1
    ra_name = 'RAJ2000'
    dec_name = 'DEJ2000'
    nbin_sample = 12
#    """ Optimistic and conservative radius """
#    radius_opt = []
#    radius_con = []
    path_aux = os.path.join(path_to_save, "hyperleda_db")
  elif options.toggle_2mb:
    # Point to 2MASS bright stars
    extfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/2mass_stars/twomass_in_des_bright.fits'
    magname = 'J'
    catname = '2mass_bright'
    maxrad = 0.11
    ra_name = 'ra'
    dec_name = 'dec'
    nbin_sample = 20
#    """ Optimistic and conservative radius """
#    radius_opt = [0.030, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.020, 0.020, 0.018, 0.018, 0.017, 0.017, 0.016, 0.016, 0.015, 0.015, 0.015, 0.014, 0.014]
#    radius_con = [0.055, 0.055, 0.050, 0.045, 0.045, 0.040, 0.040, 0.040, 0.030, 0.025, 0.025, 0.025, 0.025, 0.025, 0.024, 0.024, 0.024, 0.024, 0.023, 0.022]
    path_aux = os.path.join(path_to_save, "2mass_stars", "bright")
  elif options.toggle_2mf:
    # Point to 2MASS faint sample
    extfile = '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/2mass_stars/twomass_in_des_faint.fits'
    magname = 'J'
    catname = '2mass_faint'
    maxrad = 0.03
    ra_name = 'ra'
    dec_name = 'dec'
    nbin_sample = 10
#    """ Optimistic and conservative radius """
#    radius_opt = [0.013, 0.008, 0.006, 0.005, 0.005, 0.004, 0.004, 0.003, 0.003, 0.003]
#    radius_con = [0.015, 0.010, 0.008, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005]
    path_aux = os.path.join(path_to_save, "2mass_stars", "faint")
  else:
    print('No catalog specified. Exiting...')
    sys.exit()
  start_time = time.time()
  # Compute the radii
  calc_radii(goodobjfile, badobjfile, extfile, maxrad, nside, magname, catname, ra_name, dec_name, nbin_sample, path_aux, automated_radius, radii_to_save_opt, radii_to_save_con, magnitudes_array)
  end_time = time.time()
  elapsed_time = end_time - start_time
  print(f'Elapsed time -> {elapsed_time} secs.')
if __name__ == '__main__':

  sys.exit(main())

