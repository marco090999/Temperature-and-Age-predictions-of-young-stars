# This macro computes AV (here indicated as A0) by applying the Vergely+22 maps 
# starting from the galactic coordinates l, b and the distance dist_par of a list of stars.
# It also returns the extinctions in other bands.

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

# the following function derives the extinction in different photometric bands by applying the reddening law 
# suitable for Gaia 
# https://www.cosmos.esa.int/web/gaia/edr3-extinction-law
def extlaw(data):
    #GK = data['Gmag']-data['Ks']
    BpRp = data['phot_bp_mean_mag']-data['phot_rp_mean_mag']
	#  ../TABLES/Fitz19_EDR3_extinctionlawcoefficients/Fitz19_EDR3_MainSequence.csv
    #cG_GK = [0.993548688906439,-0.110149052160837,0.0264447715065468,-0.00571010222810317,-0.0374363031107716,0.00151447309438712,-2.52364537395156e-05,0.00623824755961677,-0.000123598316318183,-0.000158499801004388]
	# add also J and K magnitudes
    cRp_BpRp =[0.66320787941067,-0.0179847164933981,0.000493769449961458,-0.00267994405695751,-0.00651422146709376,3.30179903473159e-05,1.57894227641527e-06,-7.9800898337247e-05,0.000255679812110045,1.10476584967393e-05]
    cBp_BpRp = [1.15363197483424,-0.0814012991657388,-0.036013023976704,0.0192143585568966,-0.022397548243016,0.000840562680547171,-1.31018008013549e-05,0.00660124080271006,-0.000882247501989453,-0.000111215755291684]
    cG_BpRp = [0.995969721536602,-0.159726460302015,0.0122380738156057,0.00090726555099859,-0.0377160263914123,0.00151347495244888,-2.52364537395142E-05,0.0114522658102451,-0.000936914989014318,-0.000260296774134201]
    cJ_BpRp = [0.340345410744913,-0.00150278241491212,-0.000664573369992718,0.000417668208427758,-0.000156212937136542,1.82216155475468e-07,2.33928990110555e-09,3.41610470017311e-06,-2.13961419094946e-07,4.78655838716295e-08]
    cH_BpRp = [0.25505435103948,0.000238094101040351,-0.0009480815739429,0.000286492979762418,-6.84496601067957e-05,-2.40939063887753e-09,7.2663178194737e-10,1.87920964708149e-07,1.01009401360882e-06,1.11538757526755e-08]
    cK_BpRp = [0.19404852029171,-0.000259572240803642,0.000495770920767284,-0.00026750735610992,-2.92327041145282e-05,8.01889842444919e-09,1.22296149058604e-10,6.4675963409363e-08,1.50135637674828e-07,2.70438735386413e-09]
    #Av = 3.09*data['E(B-V)'] 
    Av=data['A0']
    #The AJ, AH and AK values are for a G2V star, using Cardelli et al. (1989) + O'Donnell (1994) extinction curve with Rv=3.1.
    # using the UBVRIJHK (Maiz-Apellaniz 2006 bessel 1990)
    #AJ=0.28688*Av
    #AH=0.18103*Av
    #AKs=0.11265*Av
    #The AJ, AH and AK values are for a G2V star, using Cardelli et al. (1989) + O'Donnell (1994) extinction curve with Rv=3.1.
    # using Pan-Starrs1
    AgP1=1.17994*Av	
    ArP1=0.86190*Av	
    AiP1=0.67648*Av	
    AzP1=0.51296*Av	
    AyP1=0.42905*Av	
    AwP1=0.83639*Av    
    AG = Av*(cG_BpRp[0]+cG_BpRp[1]*BpRp+cG_BpRp[2]*BpRp**2+cG_BpRp[3]*BpRp**3+cG_BpRp[4]*Av+cG_BpRp[5]*Av**2+cG_BpRp[6]*Av**3+
            cG_BpRp[7]*Av*BpRp+cG_BpRp[8]*Av*BpRp**2+cG_BpRp[9]*BpRp*Av**2)
    ARp = Av*(cRp_BpRp[0]+cRp_BpRp[1]*BpRp+cRp_BpRp[2]*BpRp**2+cRp_BpRp[3]*BpRp**3+cRp_BpRp[4]*Av+cRp_BpRp[5]*Av**2+
              cRp_BpRp[6]*Av**3+cRp_BpRp[7]*Av*BpRp+cRp_BpRp[8]*Av*BpRp**2+cRp_BpRp[9]*BpRp*Av**2)
    ABp = Av*(cBp_BpRp[0]+cBp_BpRp[1]*BpRp+cBp_BpRp[2]*BpRp**2+cBp_BpRp[3]*BpRp**3+cBp_BpRp[4]*Av+cBp_BpRp[5]*Av**2+
              cBp_BpRp[6]*Av**3+cBp_BpRp[7]*Av*BpRp+cBp_BpRp[8]*Av*BpRp**2+cBp_BpRp[9]*BpRp*Av**2)
    AJ = Av*(cJ_BpRp[0]+cJ_BpRp[1]*BpRp+cJ_BpRp[2]*BpRp**2+cJ_BpRp[3]*BpRp**3+cJ_BpRp[4]*Av+cJ_BpRp[5]*Av**2+cJ_BpRp[6]*Av**3+
            cJ_BpRp[7]*Av*BpRp+cJ_BpRp[8]*Av*BpRp**2+cJ_BpRp[9]*BpRp*Av**2)
    AH = Av*(cH_BpRp[0]+cH_BpRp[1]*BpRp+cH_BpRp[2]*BpRp**2+cH_BpRp[3]*BpRp**3+cH_BpRp[4]*Av+cH_BpRp[5]*Av**2+cH_BpRp[6]*Av**3+
            cH_BpRp[7]*Av*BpRp+cH_BpRp[8]*Av*BpRp**2+cH_BpRp[9]*BpRp*Av**2)
    AKs = Av*(cK_BpRp[0]+cK_BpRp[1]*BpRp+cK_BpRp[2]*BpRp**2+cK_BpRp[3]*BpRp**3+cK_BpRp[4]*Av+cK_BpRp[5]*Av**2+cK_BpRp[6]*Av**3+
            cK_BpRp[7]*Av*BpRp+cK_BpRp[8]*Av*BpRp**2+cK_BpRp[9]*BpRp*Av**2)
 
    
    data = data.assign(AG=AG, ARp=ARp, ABp=ABp, AJ=AJ, AH=AH, AKs=AKs, AgP1=AgP1, ArP1=ArP1, AiP1=AiP1, AzP1=AzP1, AyP1=AyP1)
    #Gmag0 = data['Gmag']-AG
    #Bp_Rp0 = data['BPmag']-data['RPmag']-(ABp-ARp)
    #Gmag_Kmag0 = data['Gmag']-data['Ks']-(AG-AKs)
    #Gmag_Rp0 = data['Gmag']-data['RPmag']-(AG-ARp)

    #data.add_columns([Gmag0, Bp_Rp0, Gmag_Kmag0, Gmag_Rp0], names = ['Gmag0', 'Bp_Rp0', 'G_K0', 'G_Rp0'])
    return data


outputfile='extinction_3dmap.csv'


# open the FITS file containing the Vergely et al. 2022 map
dustmap1 = '/Users/prisinzano/GAIADR3/TABLES/explore_cube_density_values_010pc_v2.fits'
dustmap2 = '/Users/prisinzano/GAIADR3/TABLES/explore_cube_density_values_025pc_v2.fits'
dustmap3 = '/Users/prisinzano/GAIADR3/TABLES/explore_cube_density_values_050pc_v2.fits'

hdul1 = fits.open(dustmap1)
hdul2 = fits.open(dustmap2)
hdul3 = fits.open(dustmap3)

mappa1 = np.transpose(hdul1[0].data, (2, 1, 0))
mappa2 = np.transpose(hdul2[0].data, (2, 1, 0))
mappa3 = np.transpose(hdul3[0].data, (2, 1, 0))

hdul1.close()
hdul2.close()
hdul3.close()

dist_max1 = 1500  # maximum distance for dustmap1
dist_max2 = 3000  # maximum distance for dustmap2

# Function to obtain the map based on the distance
def get_dustmap_and_hdul(distance):
    if distance <= dist_max1:
        return mappa1, hdul1
    elif distance <= dist_max2:
        return mappa2, hdul2
    else:
        return mappa3, hdul3

 

# input file should contain l,b, and distance of the stars , indicated in the code by 'dist_par'   . 
filedata='input_file.csv'
data = pd.read_csv(filedata, sep=',', na_values='')

# check if your database includes stars closer than 20 pc or stars  with a negative parallax
# I used the Gaia distance indicated as 'r_med_photogeo'  derived by Bailer Jones et al. 
#data.loc[data['r_med_photogeo'] < 20, 'r_med_photogeo'] = 20
#data.loc[data['parallax'] <= 0, 'parallax'] = 0.001 # serve per togliere i valori di parallax negative, assegnando un valore vicino a 0
#data['dist_par'] = data['r_med_photogeo'].fillna(1000. / data['parallax'])


extinctions=[]
for l, b, d in zip(data['l'], data['b'], data['dist_par']):
    # from galactic coordinates to cartesian coordinates
    coord = SkyCoord(l=l*u.degree, b=b*u.degree, distance=d*u.parsec, frame='galactic')
    x_star, y_star, z_star = coord.cartesian.xyz.value
#    if i == 10:
#        break
#    print(f"Elemento {i+1}: l={l}, b={b}, d={d}")
#    print(f"Elemento {i+1}: x={x_star}, y={y_star}, z={z_star}")
    #
    # In the code, we are sampling the line of sight from 0 to $d_{\star}$ using
    # s_range, which is a list of equally spaced points at steps of 10 pc along the line of sight. 
    # These points are used to evaluate the numerical integral with the trapezoidal rule, 
    # where ds represents the width of each trapezoid.
    # In particular, the interpolation function interp_func is used to 
    # interpolate the dust density values along the line of sight. 
    # This function takes as input the coordinates $(x,y,z)$ along the line of sight, 
    # which are computed as $(s*x_{\star}/d_{\star}, s*y_{\star}/d_{\star}, s*z_{\star}/d_{\star})$, 
    # where $x_{\star}$, $y_{\star}$ and $z_{\star}$ are the coordinates of the star with respect 
    # to the dust map coordinate system.
    # The result of the numerical integral is then saved in extinction_integral.
    # Compute the extinction integral
     # Select the appropriate extinction map according to the distance
    current_dustmap, current_hdul = get_dustmap_and_hdul(d)
    nx, ny, nz = current_dustmap.shape
    #step = hdul1[0].header['STEP']
    step = current_hdul[0].header['STEP']
    x_range = np.linspace(-(nx//2)*step, (nx//2)*step, nx)
    y_range = np.linspace(-(ny//2)*step, (ny//2)*step, ny)
    z_range = np.linspace(-(nz//2)*step, (nz//2)*step, nz)
    # Define an interpolator for the dust map
    # RegularGridInterpolator takes two arguments as input: 
    # the tuple of arrays defining the grid (one for each dimension) 
    # and the multidimensional array of data corresponding to the grid. 
    # Once the interpolator object is created, it can be used to compute 
    # the interpolated values at new grid points. In practice, to obtain 
    # the interpolated value at a point (x,y,z) you use the syntax 
    # interp_func((x, y, z)), where interp_func is the 
    # interpolator object created with RegularGridInterpolator.
    # The "nearest" method, instead, would be more suitable for representing 
    # discontinuities and sharp variations in density, since it assigns to each 
    # point the density of the nearest neighbor. However, it may fail to capture 
    # smooth density variations.
    interp_func = RegularGridInterpolator((x_range, y_range, z_range), current_dustmap, method='nearest', bounds_error=False, fill_value=None)
    s_range = np.linspace(0, d, int(d/10))
    ds = s_range[1] - s_range[0]
    integrand = interp_func(np.transpose([s_range*x_star/d, s_range*y_star/d, s_range*z_star/d]))    
    extinction_integral = np.sum(integrand*ds)
    extinctions.append(extinction_integral)
    #print("The extinction integral from the Sun to the star is {:.4f} mag".format(extinction_integral))

data["A0"] = extinctions
data["E(B-V)"] = data["A0"]/3.1
print(data.columns)

# END OF THE SECTION THAT COMPUTES A0, which as explained in Vergely+22 
# is defined as AV

# START OF THE SECTION THAT COMPUTES the extinction in the other bands. 
data = extlaw(data)
# Correct Gaia magnitudes for absorption in order to recalculate 
# the extinctions in the other bands using the intrinsic colors 
ABp_tmp=data['ABp']
ARp_tmp=data['ARp']
data['phot_bp_mean_mag']=data['phot_bp_mean_mag']-ABp_tmp
data['phot_rp_mean_mag']=data['phot_rp_mean_mag']-ARp_tmp
# Run the extlaw function again so that it can now apply the formula 
# as a function of the intrinsic colors
data = extlaw(data) 
# Restore the observed colors back to the original values
data['phot_bp_mean_mag']=data['phot_bp_mean_mag']+ABp_tmp
data['phot_rp_mean_mag']=data['phot_rp_mean_mag']+ARp_tmp
ABp_tmp2=data['ABp']
ARp_tmp2=data['ARp']



residuals=(ABp_tmp-ARp_tmp)-(ABp_tmp2-ARp_tmp2)
media_res = np.mean(residuals)
rms_res = np.std(residuals)
median_res = np.nanmedian(residuals)
# Compute the absolute value of the differences between the residuals and the median
differenze = np.abs(residuals - median_res)
# Compute the median of the differences
mad_res = np.nanmedian(differenze)
print(media_res,' mean')
print(rms_res,' rms')
print(median_res,' median')
print(mad_res, 'mad_res')

# Since the differences are small, only one iteration is sufficient

print(data.columns)  



data.to_csv(outputfile, sep=';', index=False)
