# 17/7/24
# compute power spectrum of true EoR from lightcone with t2c for scoring
# 
# conda activate /data-archive/sdc/venvs/t2c_21cmf/ to load the new t2c version

 
import numpy as np
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from scipy import interpolate
import os, sys
import tools21cm as t2c
from astropy.io import fits 

print (t2c.__path__)



#t2c.set_hubble_h(1)  #this is to eliminate the Mpc vs Mpc/h difference
t2c.set_hubble_h(0.6766)  
t2c.set_omega_lambda(0.69)
t2c.set_omega_matter(0.31)
t2c.set_omega_baryon(0.048)
t2c.set_sigma_8(0.82)
t2c.set_ns(0.97)

 
print('Adopted cosmology:')
print('Om',t2c.const.Omega0)
print('Ov',t2c.const.lam)
print('H0',t2c.const.H0)



#read lightcone saved as fits

filename = '/data-archive/sdc/SDC3/foregrounds/NewProducts/EoR_H21cm_truth.fits'
output_dir='/data-archive/sdc/SDC3/foregrounds/NewProducts/powerspectra/'
    
if not os.path.exists(output_dir):
    os.system('mkdir '+output_dir)

print(filename)
print(output_dir)

h = t2c.const.H0/100.

hdul = fits.open(filename)
hdul.info()
hdr=hdul[0].header
Tb_fits=hdul[0].data


#determine FoV
npix=hdul[0].header['NAXIS1']
pixreso=hdul[0].header['CDELT1']
pixreso=np.absolute(pixreso)
pixreso_arcmin=pixreso*60.
unit=hdul[0].header['CUNIT1']
halfnpix=npix/2 #use even so that this is accurate

print(npix,pixreso,unit)
if (unit == 'deg'):
    fov_deg=npix*pixreso

print('Fov [deg]=',fov_deg)    

#select the centre of field
fov_small=4
fov_name=str(fov_small)
fov_small=float(fov_small)
npix_small=npix*fov_small/fov_deg


#determine frequency range

dnu=hdul[0].header['CDELT3']
check=hdul[0].header['CRPIX3']
nfreqs=numin=hdul[0].header['NAXIS3']

if (check ==1):
    numin=hdul[0].header['CRVAL3']
    numax=numin+nfreqs*dnu

nus=np.arange(numin, numax, dnu)
nus_mhz=nus*1.e-6 # frequencies in MHz

print(Tb_fits.shape)

Tb=Tb_fits.transpose() #this is the format of the lighcone object in tools21cm, with the frequency axis as the third field

Tb_all=Tb
nus_mhz_all=nus_mhz

print(Tb.shape)

# load frequency slices
nulimits=np.loadtxt('ps_nurange.txt')
min_freqs=nulimits[:,0]
max_freqs=nulimits[:,1]


for i in range(0,len(min_freqs)):

    Tb=Tb_all
    nus_mhz=nus_mhz_all
    
    min_freq = min_freqs[i]
    max_freq = max_freqs[i]


    print('Slice frequencies (MHz):', min_freq, max_freq)

    min_freq_cdist = t2c.z_to_cdist(t2c.nu_to_z(min_freq))
    max_freq_cdist = t2c.z_to_cdist(t2c.nu_to_z(max_freq))
    print('Slice comoving distances', min_freq_cdist, max_freq_cdist)

    
    min_slice = np.abs(nus_mhz-min_freq).argmin()
    max_slice = np.abs(nus_mhz-max_freq).argmin()
    #avoid double counting the max freq pixel
    if (max_slice != 900):
        max_slice=max_slice-1 
    
    print('Slice indices:', min_slice, max_slice)
    npixels_par=max_slice-min_slice+1

    Tb=Tb[:,:,min_slice:max_slice]
    
    nus_mhz=nus_mhz[min_slice:max_slice]

    #I also need to change frequencies so that they are decreasing
    Tb_i=np.flip(Tb,axis=2)
    nus_mhz=np.flip(nus_mhz)


    Tb=Tb[int(npix/2-npix_small/2):int(npix/2+npix_small/2),int(npix/2-npix_small/2):int(npix/2+npix_small/2),:]
   

    print('transforming the lightcone into comoving units')
    data_phy, z_phy, cell_phy = t2c.observational_lightcone_to_physical(Tb, nus_mhz, pixreso_arcmin)

    print('done')
    print(data_phy.shape)



    
    kbins_fix_par=np.loadtxt('kpar_t2c_v2_augmented.txt')
    kbins_fix_per=np.loadtxt('kper_t2c_v2_augmented.txt')

    
    dkper=(kbins_fix_per[1]-kbins_fix_per[0])/2.
    dkpar=(kbins_fix_par[1]-kbins_fix_par[0])/2.


    kbins_fix_par=np.array(kbins_fix_par)-dkpar #t2c specifies the beginning of the bins, not centres
    kbins_fix_per=np.array(kbins_fix_per)-dkper

    box_perp=cell_phy*data_phy.shape[0]
    box_par=abs(min_freq_cdist-max_freq_cdist)

    print('**********')
    print('Mpc size in kpar direction',box_par)
    print('Mpc resolution in kpar direction',box_par/npixels_par)
    print('Min/max kpar',1./box_par,npixels_par/box_par)
    print('Mpc size in kper direction (4 deg)',box_perp/2.)
    print('Mpc resolution in kper direction (4deg)',cell_phy)
    print('Min/max kper',2./box_perp,1./cell_phy)
    print('**********')


    print('computing cylindrical PS')

    Pk, kper_bins, kpar_bins = t2c.power_spectrum_2d(data_phy,kbins=[kbins_fix_per,kbins_fix_par],box_dims=[box_perp,box_perp,box_par],nu_axis=2,)#window='blackmanharris')
    print('done')


    fig = plt.figure(figsize=(8,6))
    plt.pcolormesh(kper_bins,kpar_bins,np.log10(Pk).T,cmap="plasma",vmin=-4,vmax=-0.5)
    plt.title("Cylindrical PS P(k)")
    plt.colorbar()
    plt.xlabel('kper', fontdict=None, labelpad=None)
    plt.ylabel('kpar', fontdict=None, labelpad=None)
    plt.savefig(output_dir+"Pk_"+str(min_freq)+"_"+str(max_freq)+".png")
    plt.show()

    file = output_dir+"Pk_EoR_"+str(min_freq)+"_"+str(max_freq)+".txt"
    np.savetxt(file, np.transpose(Pk), fmt='%e') # transpose becausse the SDC3 convention is transpose of the t2c convention


    file = output_dir+"kpar_EoR_"+str(min_freq)+"_"+str(max_freq)+".txt"
    np.savetxt(file, kpar_bins, fmt='%e')

    file = output_dir+"kper_EoR_"+str(min_freq)+"_"+str(max_freq)+".txt"
    np.savetxt(file, kper_bins, fmt='%e')





