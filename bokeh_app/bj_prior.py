from numpy import *
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates
from zero_point import zpt
from astroquery.gaia import Gaia
import scipy.interpolate, scipy.integrate    
import pyvo, pickle
import astropy.units as u

hpxn, rlnen, alpha, beta = loadtxt('nh3d/prior_summary.csv',skiprows=1,delimiter=',',usecols=(0,5,6,7),unpack=True)
zpt.load_tables()

data = loadtxt('nh3d/prior_summary.csv',skiprows=1,unpack=True,delimiter=',')
# The parameters are reported per healpix pixel on the sky. This includes also Galactic coordinates for each pixel
# We read those to find row number we need to read other parameters from

hpxns,hpxra,hpxde = loadtxt('nh3d/HEALpixel_level5_radec_longlat_coordinates.csv',delimiter=',',skiprows=1,unpack=True,usecols=(0,3,4))
hpxcoord = astropy.coordinates.SkyCoord(hpxra,hpxde,frame='icrs',unit=u.deg)
# hpxhs = hpxn[coord.match_to_catalog_sky(hpxcoord)[0]].astype(int)

def get_bj_distance(pos):
    """docstring for get_bj_distance"""
    try:        
        res = pyvo.dal.scs.search('https://dc.zah.uni-heidelberg.de/gedr3dist/q/cone/scs.xml?',pos,2*u.arcsec)
        return res['r_med_geo'].data[0], res['r_lo_geo'].data[0], res['r_hi_geo'].data[0], res['source_id'][0]
    except:
        return nan,nan,nan, 'not found'

def get_BJ_prior(pos,maxr=5):
    """docstring for get_BJ_prior"""
    dmean, dmin, dmax = array(get_bj_distance(pos)[:3])/1000
    if isnan(dmean):
        return False
    Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
    data = Gaia.query_object_async(coordinate=pos,width=maxr*u.arcsec,height=maxr*u.arcsec) 
    # from https://gitlab.com/icc-ub/public/gaiadr3_zeropoint/-/blob/master/tutorial/ZeroPoint_examples.ipynb
    gmag = data['phot_g_mean_mag'].data.data[0]
    nueffused = data['nu_eff_used_in_astrometry'].data.data[0]
    psc = data['pseudocolour'].data.data[0]
    ecl_lat = data['ecl_lat'].data.data[0]
    soltype = data['astrometric_params_solved'].data.data[0]
    sid = (data['source_id'].data.data[0]).astype(int64)
    sid_str = data['source_id'].data.data[0]
    #to get alpha/beta/L
    phpx = int(floor(sid/562949953421312))
    hpxi = where(hpxn==phpx)[0][0]
    a, b, L = alpha[hpxi],beta[hpxi],rlnen[hpxi]
    # print(a,b,L)
    zpvals = zpt.get_zpt(gmag, nueffused, psc, ecl_lat, soltype)
    # print(zpvals)
    px = data['parallax'].data.data[0]
    pxe = data['parallax_error'].data.data[0]
    # print(dmax)
    def edr3_geoprior(r, alpha=a,beta=b, L=L):
        """parameters set for our healpix. Using info https://www2.mpia-hd.mpg.de/homes/calj/gedr3_distances/main.html"""
        from scipy.special import gamma
        return alpha/gamma((beta+1)/alpha)/(L**(beta+1))*(r*1000)**beta*exp(-((r*1000)/L)**alpha)
    def edr3_pxprob(r,px=px,pxe=pxe,wzp=zpvals):
        """Parallax and parallax error are for Gaia EDR source 5975119332093959552 which matches optical counterpart. r is in kpc"""
        return 1./(sqrt(2*pi)*pxe)*exp(-0.5/pxe**2*(px-wzp-1./r)**2)
    problem_points = [1e-6,50]+[dmin,dmax,dmean] +[dmean - (dmean-dmin)*xx for xx in linspace(0.05,5,20)] + [dmean + (dmax-dmean)*xx for xx in linspace(0.05,5,20)]
    # dmin,dmax = 0.1*dmin/1000, 5*dmax/1000
    # print(dmin,dmax)
    intgeo = scipy.integrate.quad(lambda r: edr3_geoprior(r),0,50,points=problem_points)[0]
    intpx = scipy.integrate.quad(edr3_pxprob,0,50,points=problem_points)[0]
    # print(sid,intgeo,intpx,px,pxe,zpvals,alpha[hpxi],beta[hpxi], rlnen[hpxi],hpxi)
    # here one can select which priors to use
    # Lindegren priors
    total_prob = lambda r,px=px,pxe=pxe,wzp=zpvals: (edr3_geoprior(r)/intgeo)*(edr3_pxprob(r,px=px,pxe=pxe,wzp=wzp)/intpx)
    inttot = scipy.integrate.quad(total_prob,0,50, points=problem_points)[0]
    total_pdf = lambda r,px=px,pxe=pxe,wzp=zpvals: total_prob(r,px=px,pxe=pxe,wzp=zpvals)/inttot    
    # rr = linspace(dmin,dmax,5000)
    rr = concatenate([linspace(problem_points[i],problem_points[i+1],10) for i in range(len(problem_points)-1)])
    rr.sort()
    cumprob = array([scipy.integrate.quad(total_pdf,0,x,points=problem_points)[0] for x in rr])
    ci =scipy.interpolate.interp1d(cumprob[:cumprob.searchsorted(1-1e-4)],rr[:cumprob.searchsorted(1-1e-4)]*1000,bounds_error=False,fill_value=(0,rr[cumprob.searchsorted(1-1e-4)]*1000))
    # mydprior = lambda x: ci(x) # to convert distance to kpc (required by nsx)
    offset = median(ci([0.158655,0.5,0.841345])-1000*array([dmin,dmean,dmax]))
    mydprior = lambda x: ci(x)-offset
    ii = scipy.interpolate.interp1d(linspace(0,1,1000),mydprior(linspace(0,1,1000)))
    return ii, px, pxe, mean(mydprior([0.158655,0.5,0.841345])-1000*array([dmin,dmean,dmax])), sid_str
    

    