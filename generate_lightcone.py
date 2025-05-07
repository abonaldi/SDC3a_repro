# script to generate the EoR cube for SDC3a.
# differences wrt the released map are due to new tools21cm version.
# uses custom tools21cm version with the angular_coordinates.py modified as in tools21cm_modified

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import py21cmfast as p21c
from py21cmfast import plotting
from astropy.cosmology import Planck15
from scipy import interpolate
import os, sys
import tools21cm as t2c
from astropy.io import fits 

#venv
#/data-archive/sdc/venvs/t2c_21cmf

#-------------------------------------------------------------------------------
# for interacting with the cache 
from py21cmfast import cache_tools


# set default output directory 
if not os.path.exists('_cache'):
    os.mkdir('_cache')
p21c.config['direc'] = '_cache'

# clear the cache so that we get the same result each time
cache_tools.clear_cache(direc='_cache')
tag='EoR_v9'

output_dir = './Data/Skymodels/'+tag+'/'
filename=output_dir+'lightcone_example.save'

output_dir2=output_dir
x
#-------------------------------------------------------------------------------
## we change the default level of the logger 
## so that we can see what's happening with caching 
import logging
logger = logging.getLogger("21cmFAST")
logger.setLevel(logging.INFO)

print(f"21cmFAST version is {p21c.__version__}")

#-------------------------------------------------------------------------------

# do I need to run the 21cmfast cube?
#run_all=True
run_all=False


random_seed = 9719354
ra_deg = 0.     # in degrees
dec_deg = -30.  # in degrees




def z_from_freq(nu):
    # redshift corresponding given frequency
    nu0 = 1420.40575e6
    return (nu0 /nu) - 1.

#-------------------------------------------------------------------------------
print('setting cosmology')
# 21cmFAST default cosmology is Planck 2018 
# from https://arxiv.org/pdf/1807.06209.pdf
# Table 2, last column [TT+TE+EE+lowE+lensing+BAO]

Om0 = (0.02242 + 0.11933) / 0.6766 ** 2.  #0.309641
Ob0 = 0.02242 / 0.6766 ** 2.              #0.0489747
H0 = 67.66


t2c.set_hubble_h(H0/100.)  
t2c.set_omega_matter(Om0)
t2c.set_omega_baryon(Ob0)


print('Adopted cosmology:')
print('Om',t2c.const.Omega0)
print('Ov',t2c.const.lam)
print('H0',t2c.const.H0)
print(vars(p21c.CosmoParams))

#-------------------------------------------------------------------------------

print('setting frequency properties')
# wanted frequency range and resolution 
min_freq_MHz = 106  # in MHz
max_freq_MHz = 196 # in MHz
dfreq_MHz = 0.1   # in MHz

min_freq=min_freq_MHz*1.e6
max_freq=max_freq_MHz*1.e6
dfreq=dfreq_MHz*1.e6


#freqs = np.arange(min_freq, max_freq+dfreq, dfreq)
freqs = np.arange(min_freq, max_freq, dfreq)
redshift = z_from_freq(freqs).tolist()
zmin = np.amin(redshift)
zmax = np.amax(redshift)
#print('redshifts before',redshift)


# wanted field of view at highest redshift
fov_deg = 8. # in degrees

fov_mpc=t2c.z_to_cdist(zmax)* np.deg2rad(fov_deg)
BOX_LEN = fov_mpc  # length of the box in Mpc (simulation size) : this is the comoving size # here this is 1.5 Gpc
fov_mpc_2 = t2c.deg_to_cdist(fov_deg, zmax)
print('zmax=',zmax)
print('fov_deg=',fov_deg)
print('cdist=',t2c.z_to_cdist(zmax))
print('fov_mpc check',fov_mpc,fov_mpc_2)
print('lumdist',t2c.luminosity_distance(zmax))
angular_size_deg = t2c.angular_size_comoving(BOX_LEN, zmax)

HII_DIM = 512      # number of cells per side for the low res box (output cube)
output_dtheta = (fov_deg /(HII_DIM +1.)) * 60.  # [arcmin]


# 


DIM = 3 * HII_DIM         # number of cells for the high res box (sampling initial conditions) : DIM should be at least 3 times HII_DIM

user_params = {"HII_DIM": HII_DIM, 
               "BOX_LEN": BOX_LEN, 
               "DIM": DIM,
               "N_THREADS": 44,  # keep it below 48, beyond which there is no significant gain in speed because of the scaling limitation 
               "USE_FFTW_WISDOM": True,             # True to make FFT faster
               "PERTURB_ON_HIGH_RES": True,         # True to perform the Zeldovich or 2LPT perturbation on the low or high res grid 
               "USE_INTERPOLATION_TABLES": True,    # True or code is too slow
}

if (run_all ==True):
    print('initial conditions')

    # fixing the initial condition

    initial_conditions = p21c.initial_conditions(user_params = user_params,
                                                 random_seed = random_seed, 
                                                 direc = output_dir,
    )

    print('done')
    # set necessary flags and astro parameter for the model 
    print('flag')
    flag_options = {'INHOMO_RECO': True,
                    'USE_MASS_DEPENDENT_ZETA': True, 
                    'USE_TS_FLUCT': True,
                    'USE_MINI_HALOS': False,
                    'PHOTON_CONS': True,
    }

    print('done')
    # 9 important astro params to be explored - this should be kept blind to participants
    # realization 3 with a little perturbarion
    print('astro_params')
    astro_params = {'F_STAR10': -1.509,        # [-3, 0]
                    'ALPHA_STAR': 0.496,       # [-0.5, 1.0]
                    'F_ESC10': -1.046,         # [-3, 0]
                    'ALPHA_ESC': 0.043,       # [-1.0, 0.5]
                    'M_TURN': 8.277,           # [8, 10]
                    't_STAR': 0.152,          # [0.01, 1]     these 6 are default params to reproduce the EoR history of the map model 
                        }

    lightcone_quantities = ('brightness_temp','xH_box')

    #-------------------------------------------------------------------------------
    # create a lightcone
    print('call run_lightcone')

    lightcone = p21c.run_lightcone(
        redshift = zmin,
        init_box = initial_conditions,
        flag_options = flag_options,
        astro_params = astro_params,
        lightcone_quantities = lightcone_quantities,
        global_quantities = lightcone_quantities,
        random_seed = random_seed,
        direc = output_dir,
    )

    print('done')
    #-------------------------------------------------------------------------------
    lightcone.save(filename)

lightcone = p21c.LightCone.read(filename)

lc0 = getattr(lightcone, 'brightness_temp')
print(lc0.shape) 

lc = lightcone.brightness_temp


zs0 = lightcone.lightcone_redshifts
xh=lightcone.global_xH
z_xh=lightcone.node_redshifts

print('Writing reio frac')
res = "\n".join("{} {}".format(x, y) for x, y in zip(z_xh, xh))

f = open(output_dir2+"xH_z.txt", "x")
f.write(res)
f.close()
print('Done')


# cut a lightcone in a redshift (=freq) range I want 
z_start_index = min(range(len(zs0)), key=lambda i: abs(zs0[i] - zmin))
z_end_index = min(range(len(zs0)), key=lambda i: abs(zs0[i] - zmax))
print(z_start_index, z_end_index) 

zs = zs0[z_start_index:z_end_index]
lc = lc0[:,:,z_start_index:z_end_index]
print('shape of lightcone before conversion',lc.shape) 
print('min,max z',zs.min(), zs.max())


# converting physical to observational coordinates - given cosmology is different here
angular_size_deg = t2c.angular_size_comoving(BOX_LEN, zs)
print('Minimum angular size: {:.2f} degrees'.format(angular_size_deg.min()))
print('Maximum angular size: {:.2f} degrees'.format(angular_size_deg.max()))

physical_freq = t2c.z_to_nu(zs) # redshift to frequencies in MHz
print('Minimum frequency gap in the physical light-cone data: {:.2f} MHz'.format(np.abs(np.gradient(physical_freq)).min()))
print('Maximum frequency gap in the physical light-cone data: {:.2f} MHz'.format(np.abs(np.gradient(physical_freq)).max()))


zmin = zs.min()
zmax = zs.max()

#-------------------------------------------------------------------------------
# convert the lightcone dimension 
obs_lc, obs_freq = t2c.physical_lightcone_to_observational(lc, 
                                                          zs, 
                                                          dfreq_MHz, 
                                                          output_dtheta, 
                                                          input_box_size_mpc=BOX_LEN,
                                                          mode='crop')

print(len(obs_freq)) 
print(obs_freq.min(), obs_freq.max()) 
print(obs_lc.shape) 

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

dec_deg = -30.  # in degrees
xyrefpix = HII_DIM/2.
xypix_deg = fov_deg/HII_DIM

# save to a fits file
lc_out = np.float32(obs_lc.transpose()[::-1]) #order in increasing frequency
lc_out /= 1000.  # mK to K

print(lc_out.shape)
#cut to the desired number of pixels
HII_DIM_obs=obs_lc.shape[0]

print(int(HII_DIM_obs/2-HII_DIM/2),int(HII_DIM_obs/2+HII_DIM/2))

if (HII_DIM_obs>HII_DIM):
     lc_out=lc_out[:,int(HII_DIM_obs/2-HII_DIM/2):int(HII_DIM_obs/2+HII_DIM/2),int(HII_DIM_obs/2-HII_DIM/2):int(HII_DIM_obs/2+HII_DIM/2)]


print('final shape',lc_out.shape)

sdim = '%i' % int(HII_DIM)
sfov = '%i' % int(fov_deg)

outputname = output_dir2+'obs_lightcone_'+sfov+'_'+sdim+'.fits'

# change data type
lc_out = np.float32(lc_out)

# prepare FITS header
hdu = fits.PrimaryHDU(lc_out)
hdul = fits.HDUList([hdu])

hdul[0].header.set('CTYPE1', 'RA---SIN')
hdul[0].header.set('CTYPE2', 'DEC--SIN')
hdul[0].header.set('CTYPE3', 'FREQ    ')
hdul[0].header.set('CRVAL1', ra_deg)
hdul[0].header.set('CRVAL2', dec_deg)
hdul[0].header.set('CRVAL3', min_freq)
hdul[0].header.set('CRPIX1', xyrefpix)
hdul[0].header.set('CRPIX2', xyrefpix)
hdul[0].header.set('CRPIX3', 1)
hdul[0].header.set('CDELT1', -fov_deg/HII_DIM)
hdul[0].header.set('CDELT2', fov_deg/HII_DIM)
hdul[0].header.set('CDELT3', dfreq)
hdul[0].header.set('CUNIT1', 'deg     ')
hdul[0].header.set('CUNIT2', 'deg     ')
hdul[0].header.set('CUNIT3', 'Hz      ')
hdul[0].header.set('BUNIT',  'K      ')


hdul.writeto(outputname, overwrite=True)

# BITPIX : number of bits per data pixel 
#   8 : character or unsigned binary integer 
#  16 : 16 bit twos complement binary integer
#  32 : 32 bit twos complement binary integer
# -32 : IEEE single precision floating point
# -64 : IEEE double precision floating point 

