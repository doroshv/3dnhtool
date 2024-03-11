import numpy as np
import healpy as hp
import astropy.coordinates as c
import astropy.units as u
import bj_prior
import re
import glob

fn='nh3d/ebv_fin.npz'
maps_mean = np.load(fn)['maps'].T
dbins = np.load(fn)['radius']

# h1_hi4pi = hp.read_map('nh3d/h1_HI4PI_256.hpx')
# h1_lab = hp.read_map('nh3d/h1_LAB_256.hpx')
# h1_dl = hp.read_map('nh3d/h1_DL_256.hpx')
# for x in [h1_hi4pi,h1_dl,h1_lab]:
#     x[x==hp.UNSEEN] = np.nan
    

Rv, calib_avks_mean, calib_avks_std, calib_nhag89_mean,calib_nhag89_std,calib_nhw00_mean,calib_nhw00_std = [np.load(fn)['%s'%x] for x in ['Rv', 'calib_avks_mean', 'calib_avks_std', 'calib_nhag89_mean','calib_nhag89_std','calib_nhw00_mean','calib_nhw00_std']]


import astropy.coordinates as coord

def parse_pos(posstring):
    """function to parse user input to get coordinates and possibly distance given"""
    distance_part = re.findall('@(.*)pc',posstring)
    dmin, dmax = np.nan, np.nan
    if len(distance_part)>0:
        coordinate_part = posstring[:posstring.find('@')].strip()
        try:
            dmin,dmax = re.findall('@(.*)pc',posstring)[0].split('..')
        except:
            dmax = re.findall('@(.*)pc',posstring)[0]
            dmin = dmax.replace('k','')
        if dmax[-1]=='k':
            multfac =1000
            dmin,dmax = float(dmin)*multfac, float(dmax.replace('k',''))*multfac
        else:
            multfac = 1.0
            dmin,dmax = float(dmin)*multfac, float(dmax)*multfac
    else:
        coordinate_part = posstring
    #need to take special care for case when galactic coordinates are given
    if coordinate_part.lower().find('gal')>=0:
        try:
            l, b = [float(x) for x in coordinate_part.replace('gal','').strip().split()]
            pos = c.SkyCoord(l,b,frame='galactic',unit=u.deg).transform_to(c.ICRS) # everything in ICRS
        except:
            return False
    else:
        try:
            pos = c.SkyCoord(coordinate_part, unit=(u.deg, u.deg))
        except:
            try:
                pos = c.SkyCoord(coordinate_part, unit=(u.hourangle, u.deg))
            except:
                try:
                    pos = c.get_icrs_coordinates(coordinate_part)
                except:
                    return False
    return pos, dmin, dmax            

def query_map(pos, nsim=1000):
    """Query given position and return line of sight E(B-V), AV, AK, and NH (wilms abundances) in given direction including uncertainties"""
    gpos = pos.transform_to(c.Galactic)
    # print(gpos)
    gpix = np.array(hp.ang2pix(256,gpos.l.deg,gpos.b.deg,lonlat=True))
    #this function is for single position, so there can be only one unique gpix
    assert np.unique(gpix).shape[0]==1, "function is designed to query single sky position (even with multiple distances)"
    try:
        dist = np.array(gpos.distance.pc)
    except:
        dist = np.array([dbins[-1]])
    #account for scutter associated with distance uncertainty
    # ebv   = maps_mean[gpix,np.minimum(dbins.searchsorted(dist),len(dbins))]
    ebv   = maps_mean[gpix[0],:]
    av = Rv*ebv 
    add_err = calib_avks_mean[0]
    rel_err = calib_avks_mean[1]
    av2ak = calib_avks_mean[3]
    ave = np.sqrt(add_err**2+(av*rel_err)**2)
    # avsim = np.random.normal(loc=av,scale=ave)
    ebve = ave/Rv
    ebve_opt = ebv*0.1355 # estimated scatter between E23 and GNILC
    ebve = (ebve*(1-ebv/ebv[-1])+ebve_opt*ebv/ebv[-1])
    ak = av*av2ak
    ake = ave*av2ak
        
    # Use W00 abundances by default
    nh_conv_fac = calib_nhw00_mean[0]
    nh_add_err = calib_nhw00_mean[1]
    nh_rel_err = calib_nhw00_mean[2]
    nh = ebv*nh_conv_fac
    nhe = np.sqrt(nh_add_err**2+(nh*nh_rel_err)**2)

    # now we can calculate relevant LOS quantities (i.e. for all radial bins)
    los = (np.array((ebv-ebve,ebv,ebv+ebve)),
            np.array((av-ave,av,av+ave)),
            np.array((ak-ake,ak,ak+ake)),
            np.array((nh-nhe,nh, nh+nhe)))
    #set minimal absorbtion to zero
    for x in los:
        for y in x:
            y[y<=0]=1e-10 # should not be zero otherwise plot screw-up
    idx = np.minimum(dbins.searchsorted(dist),len(dbins)-1)
    #and those integrated up to certain distance (distance can be randomly sampled according to BJ priors or otherwise hence here we need quantiles)
    ebvsim = np.random.normal(loc=ebv[idx],scale=ebve[idx])
    avsim = np.random.normal(loc=av[idx],scale=ave[idx])
    aksim = np.random.normal(loc=ak[idx],scale=ake[idx])
    nhsim = np.random.normal(loc=nh[idx],scale=nhe[idx])
    cum = ( np.quantile(ebvsim,[0.15865,0.5,1-0.15865]),
            np.quantile(avsim, [0.15865,0.5,1-0.15865]), 
            np.quantile(aksim, [0.15865,0.5,1-0.15865]),
            np.quantile(nhsim, [0.15865,0.5,1-0.15865]))
    # NH is in units 1e21
    return los, cum#, [x[gpix] for x in [h1_hi4pi,h1_lab,h1_dl]]


def get_bj_posterior(name, nsim=1000):
    """docstring for get_pos_from_name"""
    parsed_pos = parse_pos(name)
    sid = 'not found'
    if parsed_pos:
        pos, dmin, dmax = parsed_pos
    else:
        return False
    if np.isnan(dmin*dmax):        
        try:
            bjp = bj_prior.get_BJ_prior(pos)
            dsim = bjp[0](np.random.uniform(0,1,size=nsim))
            dlow, dmean, dhi, sid = *np.quantile(dsim,[0.15864,0.5,1-0.15865]), bjp[-1]
        except:
            try:
                dmean, dlow, dhi, sid = bj_prior.get_bj_distance(pos)
                dsim = np.random.uniform(dlow,dhi,size=nsim)
            except:
                dsim = [dbins[-1]]*nsim
                dlow = dbins[-1]
                dhi = dbins[-1]
                dmean = dbins[-1]
                sid = 'not found'
    else:
        dmean, dlow, dhi = 0.5*(dmin+dmax), dmin, dmax
        dsim = np.random.uniform(dlow,dhi,size=nsim)
    posnew = c.SkyCoord(pos,distance=dsim*u.pc)
    los, cum = query_map(posnew, nsim=nsim)
    # return posnew, dlow, dmean,dhi, dsim, sid, pos
    return cum, los, dlow, dmean,dhi, dsim, sid, pos#, h1