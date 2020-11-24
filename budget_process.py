### Exectutalbe gernralised budget script


import numpy as np
import pickle
import datetime as dt
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from sys import path
import copy
import os
path.insert(0, '/Users/H/WAVES/geo_data_group/')
from mpl_toolkits.basemap import Basemap
# import data_year as dy
import grid_set as gs
import budget_inputs as bi
import budget as b


## outputs
spath = '/Users/H/PREMELT/Budget/Outputs/out_test/'

save_dir = os.path.dirname(spath)
# check it exists
if not os.path.exists(save_dir):
    # make if it doesn't
    os.makedirs(save_dir)
    print('Creating directory: ',save_dir)
else: print('Existing directory: ',save_dir)

### copy the option from this file to spath
write_l = False
with open("./budget_process.py") as f:
    with open(spath+"config.txt", 'w+') as f1:
        for l in f:
            if len(l)>10:
                if l.split()[0]=="#####":
                    write_l = True
            if write_l and l[0:2]!="# ":
                f1.write(l)
            if "END USER" in l and write_l:
                break

##### USER OPTIONS #####
## dates

ts = dt.datetime(2014,3,11)
te = dt.datetime(2015,4,16)
tw = relativedelta(days = 1)

## message option
load_verbos = False
hole_verbos = False

## diagnostic plots - for the first time point
print_inputs = False
print_budget = False
print_budget_square = False
print_dstats = False
print_history = False
print_history_square = False

## deformation options
calc_dstats = False

## masking options
mask_poleV = False
mask_poleT = False
pole_hole_lat = 88.0

## Thickness/Con limit masking
Tlim = 0.2   ### meters below which we mask
Clim = 0.15  ### conc.  below which we mask

### Time smoothing velocity option
V_time_smooth = False

### First we need grids for the data we're using
## the grids sit with a projection North pole in this case


m = Basemap(projection='stere', 
            lon_0=-45.0, lat_0=83, lat_ts=0, 
            height = 3335000*1.8, width = 3335000*1.8)

### VELOCITY
GV = gs.grid_set(m)
# GV.load_grid('grids/PIOMAS_gs.npz')
# GV.load_mask('grids/PIOMAS_gs_mask.npz')
GV.load_grid('grids/Pathfinder_gs.npz')
GV.load_mask('grids/Pathfinder_gs_mask.npz')
# GV.load_grid('grids/Kimura_gs.npz')
# GV.load_mask('grids/Kimura_gs_mask.npz')
# GV.load_grid('grids/osisaf_gs.npz')
# GV.load_mask('grids/osisaf_gs_mask2.npz')
GV.reproject(m)
if mask_poleV:
    GV.mask[GV.lats>pole_hole_lat] = np.nan

### THICKNESSS 
GT = gs.grid_set(m)
GT.load_grid('grids/AWI_gs.npz')
GT.load_mask('grids/AWI_gs_mask.npz')
GT.reproject(m)
# # if mask_poleT:
#     GT.mask[GT.lats>pole_hole_lat] = np.nan


#### CONCENTRATION
GC = gs.grid_set(m)
# GC.load_grid('grids/AWI_gs.npz')
# GC.load_mask('grids/AWI_gs_mask.npz')
GC.load_grid('grids/NSIDC_gs.npz')
GC.load_mask('grids/NSIDC_gs_mask.npz')
GC.reproject(m)

##### Set up budget input objects
### all the code is in the budget_inputs.py
### first we query the dates of interest
# path = '/Volumes/BU_extra/BUDGET/PIOMAS_data/'
path = '/Volumes/BU_extra/BUDGET/Pathfinderv4/'
# path = '/Volumes/BU_extra/BUDGET/Data/Drift/Kimura/All_binary_data/'
# path = '/Volumes/BU_extra/BUDGET/Data/Drift/OSISAF/osisaf.met.no/archive/ice/drift_lr/merged/'
# P = bi.PIOMAS(path)
P = bi.Pathfinder(path)
# P = bi.Kimura(path)
# P = bi.OSISAF(path)

P.get_dates(ts,te)

### THICKNESS
path = '/Volumes/BU_extra/CryoSat/AWI_thickness/monthly/'
# path = '/Volumes/BU_extra/BUDGET/Data/Thickness/CPOM/OfficialGridded/'
A = bi.AWI_monthly(path)
# A = bi.AWI_SMOS_daily(path)
# A = bi.CPOM_hi(path,GV)

A.get_dates(ts,te)


### CONCENTRATION
nt_path = '/Volumes/BU_extra/NSIDC/daily/'
C = bi.NSIDC_nt(nt_path)

C.get_dates(ts,te)

### use data to set the dates to process
dlist = b.get_loop_list(A.dates,tw)
ts = dlist[0]
ndays = np.shape(dlist)[0]

### compare to additional data
# dlist = b.comp_dates(dlist,P.dates)


# input select
InputV = P # velocity forcing
InputT = A # thickness
InputC = C # concentration

# option for whether we're using effective thickness or not
make_eff_thk = False

# grid select
Gvel = GV
Gthk = GT
Gcon = GC

rgd_thk = True
rgd_con = True

# smoothing options
smthvel = False
smththk = False
smthcon = False

#### old uniform
# velsmrad = 6
# Vsmth = lambda x: b.smooth2a(x,velsmrad)
# thksmrad = 4
# Tsmth = lambda x: b.smooth2a(x,thksmrad)
# consmrad = 1
# Csmth = lambda x: b.smooth2a(x,consmrad)

#### new vary smoothing
velsmthdist = 200e3
thksmthdist = 100e3
consmthdist = 60e3

VSmthObj = gs.geo_vary_smooth(Gvel,velsmthdist,verbos=True)
TSmthObj = gs.geo_vary_smooth(Gvel,thksmthdist,verbos=True)
CSmthObj = gs.geo_vary_smooth(Gvel,consmthdist,verbos=True)

Vsmth = VSmthObj.smooth
Tsmth = TSmthObj.smooth
Csmth = CSmthObj.smooth

#### history fields
budhist = b.budhist(Gvel,dlist,tw,drift_stats = calc_dstats,
                   )
#                     dumptw = relativedelta(days=20),budday0=ts)

### optional history regrid
histrgd = gs.Gs2Gs(GV,GC,vectors=True)

#### END USER OPTIONS


##### below here should be left alone (should)

# regrid objects
if rgd_thk:
    print("Building thickness regridder")
    print(InputT.name+" to "+InputV.name)
    RGDthk = gs.Gs2Gs(Gthk,Gvel)
if rgd_con:
    print("Building concentration regridder")
    print(InputC.name+" to "+InputV.name)
    RGDcon = gs.Gs2Gs(Gcon,Gvel)



Uslices = np.empty([3,Gvel.m,Gvel.n])
Vslices = np.empty([3,Gvel.m,Gvel.n])
Tslices = np.empty([3,Gvel.m,Gvel.n])
Cslices = np.empty([3,Gvel.m,Gvel.n])
Vol_slices = np.empty([3,Gvel.m,Gvel.n])

# divergence advection are runnign means so make array
budadv = np.empty([3,Gvel.m,Gvel.n])
buddiv = np.empty([3,Gvel.m,Gvel.n])

### initalising the slices

Utemp,Vtemp = b.get_vels_array(InputV,ts,tw,verbos=load_verbos)
Utemp = Utemp*Gvel.mask
Vtemp = Vtemp*Gvel.mask
if smthvel:
    Uslices[0] = Vsmth(Utemp)
    Vslices[0] = Vsmth(Vtemp)
else:
    Uslices[0] = Utemp
    Vslices[0] = Vtemp
Ttemp= b.get_hi_array(InputT,ts,tw,verbos=load_verbos)
if print_inputs: Tplot = copy.copy(Ttemp)
Ttemp[Ttemp <Tlim] = np.nan
Ttemp = Ttemp*Gthk.mask
if rgd_thk:
    # call the regrid object
    Ttemp = RGDthk.rg_array(Ttemp)
if smththk:
    Tslices[0] = Tsmth(Ttemp)
else:
    Tslices[0] = Ttemp
Ctemp= b.get_aice_array(InputC,ts,tw,verbos=load_verbos)
if print_inputs: Cplot = copy.copy(Ctemp)
Ctemp[Ctemp <Clim] = np.nan
Ctemp = Ctemp*Gcon.mask
if rgd_con:
    # call the regrid object
    Ctemp = RGDcon.rg_array(Ctemp)
if smthcon:
    Cslices[0] = Csmth(Ctemp)
else:
    Cslices[0] = Ctemp

if print_inputs: b.plot_inputs(Gvel,Gthk,Gcon,InputV,InputT,InputC,
                               Utemp,Vtemp,Tplot,Cplot,spath,ts)
    
if make_eff_thk:
    Vol_slices[0] = Tslices[0]*Cslices[0]
else:
    Vol_slices[0] = Tslices[0]

Uslices[1],Vslices[1] = Uslices[0],Vslices[0]
Tslices[1] = Tslices[0]
Cslices[1] = Cslices[0]
Vol_slices[1] = Vol_slices[0]

#### adv & div at single points
# advection and divergence are calculated slice by slice
# so do it as we load
# ADVECTION
advtemp = b.advection(Uslices[0],Vslices[0],Vol_slices[0],Gvel)
budadv[0] = advtemp

# DIVERGENCE
divtemp = b.divergence(Uslices[0],Vslices[0],Vol_slices[0],Gvel)
buddiv[0] = divtemp

budadv[1] = budadv[0]
buddiv[1] = buddiv[0]

### Now loop along the list of dates
for n,d in enumerate(dlist[1:]):
    ### check list for 
    budhist.check_dump(n,d,hole_verbos)
    
    print("new slice from "+d.strftime('%Y%m%d'))
    # update final slice
    
    Utemp,Vtemp = b.get_vels_array(InputV,d,tw,verbos=load_verbos)
    Utemp = Utemp*Gvel.mask
    Vtemp = Vtemp*Gvel.mask
    if smthvel:
        Uslices[2] = Vsmth(Utemp)
        Vslices[2] = Vsmth(Vtemp)
    else:
        Uslices[2] = Utemp
        Vslices[2] = Vtemp
    Ttemp= b.get_hi_array(InputT,d,tw,verbos=load_verbos)
    Ttemp[Ttemp <Tlim] = np.nan
    Ttemp = Ttemp*Gthk.mask
    if rgd_thk:
        # call the regrid object
        Ttemp = RGDthk.rg_array(Ttemp)
    if smththk:
        Tslices[2] = Tsmth(Ttemp)
    else:
        Tslices[2] = Ttemp
    Ctemp= b.get_aice_array(InputC,d,tw,verbos=load_verbos)
    Ctemp[Ctemp <Clim] = np.nan
    Ctemp = Ctemp*Gcon.mask
    if rgd_con:
        # call the regrid object
        Ctemp = RGDcon.rg_array(Ctemp)
    if smthcon:
        Cslices[2] = Csmth(Ctemp)
    else:
        Cslices[2] = Ctemp

    #### actual budget
    if make_eff_thk:
        Vol_slices[2] = Tslices[2]*Cslices[2]
    else:
        Vol_slices[2] = Tslices[2]
    
    ### time smooth the velocity before doing the budget
    if V_time_smooth:    
        UTsmth = np.nanmean(Uslices,axis=0)
        VTsmth = np.nanmean(Vslices,axis=0)
        VolTsmth = np.nanmean(Vol_slices,axis=0)
        
        budadv3 =  b.advection(UTsmth,VTsmth,VolTsmth,Gvel)
        buddiv3 = b.divergence(UTsmth,VTsmth,VolTsmth,Gvel)
        
    else:
        # advection and divergence are calculated slice by slice
        # so do it as we load
        # ADVECTION
        advtemp = b.advection(Uslices[2],Vslices[2],Vol_slices[2],Gvel)
        budadv[2] = advtemp

        # DIVERGENCE
        divtemp = b.divergence(Uslices[2],Vslices[2],Vol_slices[2],Gvel)
        buddiv[2] = divtemp

    #     print(np.nanmean(buddif))
        # we also need to smooth the div and adv over 3 time steps
        buddiv3 = np.nanmean(buddiv,axis=0)
        budadv3 = np.nanmean(budadv,axis=0)
    
    # CHANGE IN VOL USES ALL SLICES    
    nsecs = b.get_step_size(d,tw)
    buddif=b.intensification(Vol_slices,nsecs)
    
    if n==0 and print_budget: b.plot_budget(Gvel,buddif,buddiv3,
                                            budadv3,Utemp,Vtemp,spath,ts)
    if n==0 and print_budget_square: b.plot_budget_square(Gvel,buddif,buddiv3,
                                            budadv3,Utemp,Vtemp,spath,ts)
    
    budhist.accumulate_budget(d,buddif,budadv3,buddiv3,
                              Tslices[1],Cslices[1],Uslices[1],Vslices[1],
                              hole_verbos)
    
    if calc_dstats:
        ### get components
        Dudx, Dudy = gs.geo_gradient(Uslices[1]*Gvel.mask,Gvel)
        Dvdx, Dvdy = gs.geo_gradient(Vslices[1]*Gvel.mask,Gvel)
        ddiv = Dudx + Dvdy
        dcrl = Dvdx - Dudy
        dshr = Dvdx + Dudy
        dshr = np.sqrt((Dudx + Dvdy)**2 + dshr**2)
        budhist.accumulate_dstats(ddiv,dcrl,dshr)
        if n==0 and print_dstats: b.plot_dstats(Gvel,Utemp,Vtemp,
                                  InputV,ddiv,dcrl,dshr,spath,ts)
    
    # update slices
    Uslices[0] = Uslices[1]
    Uslices[1] = Uslices[2]
    Vslices[0] = Vslices[1]
    Vslices[1] = Vslices[2]
    
    Tslices[0] = Tslices[1]
    Tslices[1] = Tslices[2]
    Cslices[0] = Cslices[1]
    Cslices[1] = Cslices[2]
    
    Vol_slices[0] = Vol_slices[1]
    Vol_slices[1] = Vol_slices[2]
    
    budadv[0] = budadv[1]
    budadv[1] = budadv[2]
    
    buddiv[0] = buddiv[1]
    buddiv[1] = buddiv[2]
    
    # also need to dump on end of list
    # also need to dump when we reach a hole 
    if budhist.hist_dump_now:
        if print_history:
            b.plot_budget(Gvel,budhist.buddifhist,
                          budhist.buddivhist,budhist.buddivhist,
                          budhist.budxvelhist/budhist.hist_count,
                          budhist.budyvelhist/budhist.hist_count,
                          spath,d)
        if print_history_square:
            b.plot_budget_square(Gvel,budhist.buddifhist,
                          budhist.buddivhist,budhist.buddivhist,
                          budhist.budxvelhist/budhist.hist_count,
                          budhist.budyvelhist/budhist.hist_count,
                          spath,d)
        ### optional second dump
        budhist.dump_nc(spath,d,Gvel,InputV,InputT,InputC,
                       outrgd = [histrgd,'_EASE',GC],accu_reset=False)
        
        budhist.dump_nc(spath,d,Gvel,InputV,InputT,InputC)
#         budhist.dump_pkl(spath,d,Gvel)

    if budhist.hist_dump_now and calc_dstats:
        budhist.dump_dstats_nc(spath,d,Gvel,InputV)
#         budhist.dump_pkl(spath,d,dstats = calc_dstats)
        
