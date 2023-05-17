import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


mag_j,radius_b = np.loadtxt('T2_radius_bright',unpack=True)
print mag_j,radius_b

z = np.polyfit(mag_j,radius_b,2)

pt2b = np.poly1d(z)

plt.plot(mag_j,radius_b,'.',label='T2 config',c='orange')

plt.plot(mag_j,pt2b(mag_j),'--',label='T2 config fit',c='orange')

#plt.show()
print 'T2 bright polinomial ',z
mag_j,radius_b = np.loadtxt('T2_radius_faint',unpack=True)

print mag_j,radius_b

z = np.polyfit(mag_j,radius_b,2)
pt2f = np.poly1d(z)

print z
print 'T2 faint polinomial ',z
plt.plot(mag_j,radius_b,'.',c='orange')
plt.plot(mag_j,pt2f(mag_j),'--',c='orange')

mag_j,radius_b = np.loadtxt('T3_radius_bright',unpack=True)
print mag_j,radius_b
z = np.polyfit(mag_j,radius_b,2)
pt3b = np.poly1d(z)
plt.plot(mag_j,radius_b,'.',label='T3 config', c='purple')
plt.plot(mag_j,pt3b(mag_j),'--',label='T3 config fit', c='purple')
#plt.show()
print z
print 'T3 bright polinomial ',z
mag_j,radius_b = np.loadtxt('T3_radius_faint',unpack=True)
print mag_j,radius_b
z = np.polyfit(mag_j,radius_b,2)
pt3f = np.poly1d(z)
print z
print 'T3 faint polinomial ',z
plt.plot(mag_j,radius_b,'.',c='purple')
plt.plot(mag_j,pt3f(mag_j),'--',c='purple')

plt.xlabel('J 2mass magnitude')
plt.ylabel('Radius [degrees]')
plt.legend()
plt.savefig('2mass_foreground_T2T3_config.png')


'''
import astropy.io.fits as pf

bright_2mass = pf.open('2mass_stars_jlt8_y3foot.fits')[1].data
jmag = bright_2mass['jmag']

radius_t2 = []
radius_t3 = []
area_t2 = []
area_t3 = []

for j in jmag:
   if j<5.85:
	radius_t2.append(0.03)
	radius_t3.append(0.03)
	area_t2.append(0.03*0.03*np.pi)
	area_t3.append(0.03*0.03*np.pi)
   else:
	radius_t2.append(pt2b(j))
	radius_t3.append(pt3b(j))	
	area_t2.append(pt2b(j)*pt2b(j)*np.pi)
	area_t3.append(pt3b(j)*pt3b(j)*np.pi)


cols = []
cols.append(pf.Column(name='radius_t2', format='D', array=np.array(radius_t2)))
cols.append(pf.Column(name='radius_t3', format='D', array=np.array(radius_t3)))
cols.append(pf.Column(name='area_t2', format='D', array=np.array(area_t2)))
cols.append(pf.Column(name='area_t3', format='D', array=np.array(area_t3)))

orig_cols = bright_2mass.columns
new_cols = pf.ColDefs(cols)
hdu = pf.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto('2mass_stars_jlt8_y3foot_radius.fits')

faint_2mass = pf.open('2mass_stars_jgt8_y3foot.fits')[1].data
jmag = faint_2mass['jmag']

area_t2 = []
area_t3 = []

radius_t2 = []
radius_t3 = []
for j in jmag:
   if j<9.49:
	radius_t2.append(0.01)
	
	radius_t3.append(0.013)
	area_t2.append(0.01*0.01*np.pi)
	area_t3.append(0.013*0.013*np.pi)

   else:
	radius_t2.append(pt2f(j))
	radius_t3.append(pt3f(j))
	area_t2.append(pt2f(j)*pt2f(j)*np.pi)
        area_t3.append(pt3f(j)*pt3f(j)*np.pi)


cols = []
cols.append(pf.Column(name='radius_t2', format='D', array=np.array(radius_t2)))
cols.append(pf.Column(name='radius_t3', format='D', array=np.array(radius_t3)))
cols.append(pf.Column(name='area_t2', format='D', array=np.array(area_t2)))
cols.append(pf.Column(name='area_t3', format='D', array=np.array(area_t3)))
orig_cols = faint_2mass.columns
new_cols = pf.ColDefs(cols)
hdu = pf.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto('2mass_stars_jgt8_y3foot_radius.fits')
'''