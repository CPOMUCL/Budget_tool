#! /home/hds/.conda/envs/harry/bin/python3
### Exectutalbe gernralised budget script


import numpy as np
import pickle
import datetime as dt
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from sys import path
import copy
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import NearestNDInterpolator
# path.insert(0, '/Users/H/WAVES/geo_data_group/')
# from mpl_toolkits.basemap import Basemap
# import data_year as dy
import grid_set as gs
import budget_inputs as bi
# import budget as b
import budget_accum as b

import sys
pFile = sys.stdout

## outputs
# spath = '/home/hds/DATA2/Budgets/Budget_out_2021_06_23/out_OSISAF_nt_AWI_SMOS/'
# spath = '/home/hds/DATA2/Budgets/Budget_out_2022_07_19/out_Pathfinder_nt_Bristol_full/'
# spath = '/home/hds/DATA2/Budgets/Budget_out_2022_07_19/out_Pathfinder_nt_Bristol_full_smooth/'
# spath = '/home/hds/DATA2/Budgets/Budget_out_2022_07_19/out_Pathfinder_nt_Bristol_full_daily_smooth_test/'
spath = '/home/hds/DATA2/Budgets/Budget_out_2023_06_15/out_OSISAF_nt_Bristol_full_year_min/'


save_dir = os.path.dirname(spath)
# check it exists
if not os.path.exists(save_dir):
    # make if it doesn't
    os.makedirs(save_dir)
    print('Creating directory: ',save_dir)
else: print('Existing directory: ',save_dir)

### copy the option from this file to spath
write_l = False
with open("./budget_process_seasons_uncert.py") as f:
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

ts = dt.datetime(2010,10,14)
te = dt.datetime(2021,5,1)
# ts = dt.datetime(2014,10,1)
# te = dt.datetime(2015,10,1)
tw = relativedelta(days = 1)

## message option
load_verbos = False
hole_verbos = False

#### history field options
hist_month_maps = False
hist_month_gates= False
#### set to false to remove the original data for publishing
save_input = True


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
## improved hole masking...
mask_poleV = True
mask_poleT = False
pole_hole_lat = 87.4

## Mask Arctic Only
arco_lat_lim = 50

## New ice in the budget or not set to True to record growth only from 20cm
cut_new_ice = True

## Thickness/Con limit masking
Tlim = 0.2   ### meters below which we split into data and fill
Clim = 0.15  ### conc.  below which we mask

### Time smoothing velocity option
V_time_smooth = False

### First we need grids for the data we're using
## the grids sit with a projection North pole in this case


# m = Basemap(projection='stere', 
#             lon_0=-45.0, lat_0=83, lat_ts=0, 
#             height = 3335000*1.8, width = 3335000*1.8)
# m = ccrs.RotatedPole(pole_latitude=00.0)
m = ccrs.LambertAzimuthalEqualArea(central_latitude=90)
### VELOCITY
GV = gs.grid_set(m)
# GV.load_grid('grids/PIOMAS_gs.npz')
# GV.load_mask('grids/PIOMAS_gs_mask.npz')
# GV.load_grid('grids/Pathfinder_gs.npz')
# GV.load_mask('grids/Pathfinder_gs_mask.npz')
# GV.load_grid('grids/Kimura_gs.npz')
# GV.load_mask('grids/Kimura_gs_mask.npz')
# GV.load_grid('grids/osisaf_gs.npz')
# GV.load_mask('grids/osisaf_gs_mask2.npz')
GV.load_grid('grids/osisaf_455_gs.npz')
GV.get_grid_mask()
if mask_poleV:
    GV.mask[GV.lats>pole_hole_lat] = np.nan

### THICKNESSS 
GT = gs.grid_set(m)
# GT.load_grid('grids/AWI_gs.npz')
# GT.load_mask('grids/AWI_gs_mask.npz')
# GT.load_grid('grids/Bristol_gs.npz')
# GT.load_mask('grids/Bristol_mask_gs.npz')
GT.load_grid('grids/Bristol_nh_80km_gs.npz')
GT.load_mask('grids/Bristol_nh_80km_gs_mask.npz')
if mask_poleT:
    GT.mask[GT.lats>pole_hole_lat] = np.nan


#### CONCENTRATION
GC = gs.grid_set(m)
# GC.load_grid('grids/AWI_gs.npz')
# GC.load_mask('grids/AWI_gs_mask.npz')
GC.load_grid('grids/NSIDC_gs.npz')
GC.load_mask('grids/NSIDC_gs_mask.npz')

### extra grids for history
Gnsidc = gs.grid_set(m)
Gnsidc.load_grid('grids/NSIDC_gs.npz')
Gnsidc.load_mask('grids/NSIDC_gs_mask.npz')

### Gate grids
save_gate_grids=False
gate_grid_dir = './grids/Gate_grids/'
gate_dir = os.path.dirname(gate_grid_dir)
# check it exists
if not os.path.exists(gate_dir) and save_gate_grids:
    # make if it doesn't
    os.makedirs(gate_dir)

### Fram strait
Fram = gs.grid_set(m)
### Barents  out
Barents1 = gs.grid_set(m)
### Barents In
Barents2 = gs.grid_set(m)
### Kara In
Kara1 = gs.grid_set(m)
### Kara Out
Kara2 = gs.grid_set(m)
### Chukchi  out
Chukchi = gs.grid_set(m)
### Baffin Top
BaffinTop = gs.grid_set(m)
### Baffin Bot
BaffinBot = gs.grid_set(m)
### Central_west
Central_west = gs.grid_set(m)
### Central_east
Central_east = gs.grid_set(m)
### Chukchi_Beaufort
Chukchi_Beaufort = gs.grid_set(m)

Gates_list = [
    [Fram,'Gate_Fram/'],
    [Barents1,'Gate_BarentsBot/'],
    [Barents2,'Gate_BarentsTop/'],
    [Kara1,'Gate_KaraTop/'],
    [Kara2,'Gate_KaraBot/'],
    [Chukchi,'Gate_Chukchi/'],
    [BaffinTop,'Gate_BaffinTop/'],
    [BaffinBot,'Gate_BaffinBot/'],
    [Central_east,'Gate_Central_east/'],
    [Central_west,'Gate_Central_west/'],
    [Chukchi_Beaufort,'Gate_Chukchi_Beaufort/'],
]

for gl in Gates_list:
    if save_gate_grids:
        gl[0].get_grid_info(av_ang=False)
        grid_name = gl[1].split('/')[0]+'.npz'
        gl[0].save_grid(gate_grid_dir+grid_name)
        mask_name = gl[1].split('/')[0]+'_mask.npz'
        gl[0].save_mask(gate_grid_dir+mask_name)
    else:
        grid_name = gl[1].split('/')[0]+'.npz'
        gl[0].load_grid(gate_grid_dir+grid_name)
        mask_name = gl[1].split('/')[0]+'_mask.npz'
        gl[0].load_mask(gate_grid_dir+mask_name)


##### Set up budget input objects
### all the code is in the budget_inputs.py
### first we query the dates of interest
# path = '/home/hds/DATA/Pathfinder/v4/'
# path = '/home/hds/DATA2/PIOMAS/'
path = '/home/hds/DATA/OSISAF/OSI-455/'
# path = '/Volumes/BU_extra/BUDGET/PIOMAS_data/'
# path = '/Volumes/BU_extra/BUDGET/Pathfinderv4/'
# path = '/Volumes/BU_extra/BUDGET/Data/Drift/Kimura/All_binary_data/'
# path = '/Volumes/BU_extra/BUDGET/Data/Drift/OSISAF/osisaf.met.no/archive/ice/drift_lr/merged/'
# P = bi.PIOMAS(path)
# P = bi.Pathfinder(path)
# P = bi.Kimura(path)
P = bi.OSISAF(path,version='455')

P.get_dates(ts,te)

### THICKNESS
# path = '/home/hds/DATA2/PIOMAS/'
# path = '/home/hds/DATA2/CPOM_SIT/'
# path = '/home/hds/DATA/AWI_SIT/SMOS_daily/'
# path = '/home/hds/DATA2/ESA_CCI_SIT/'
path = '/home/hds/DATA2/Bristol_SIT/'
# path = '/Volumes/BU_extra/CryoSat/AWI_thickness/monthly/'
# path = '/Volumes/BU_extra/BUDGET/Data/Thickness/CPOM/OfficialGridded/'
# A = bi.AWI_monthly(path)
# A = bi.AWI_SMOS_daily(path)
# A = bi.CPOM_hi(path,GV)
# A = bi.ESA_CCI(path,'ENVISAT')
# A = bi.CPOM_hi(path,GV)
A = bi.Bristol_thickness_seasonal(path)
# A = bi.Bristol_thickness(path,var = 'Sea_Ice_Thickness_incSMOS')
# A = bi.PIOMAS(path)

A.get_dates(ts,te)


### CONCENTRATION
path = '/home/hds/DATA2/NSIDC/nasa_daily/arco/'
C = bi.NSIDC_nt(path)
# path = '/home/hds/DATA2/NSIDC/bootstrap_daily/'
# C = bi.NSIDC_bt(path)

C.get_dates(ts,te)

### use data to set the dates to process
dlist = b.get_loop_list(P.dates,tw)
ts = dlist[0]
ndays = np.shape(dlist)[0]

# dlist = []
# d0 = dt.datetime(2014,10,1)
# [dlist.append(d0+relativedelta(days=d)) for d in range(365) ]
# ts = dlist[0]

### compare to additional data
# dlist = b.comp_dates(dlist,A.dates)
# dlist = b.comp_dates(dlist,C.dates)
# dlist = b.comp_dates(dlist,P.dates)
dlist = dlist[1:]
ts = dlist[0]

# input select
InputV = P # velocity forcing
InputT = A # thickness
InputC = C # concentration

# option for whether we're using effective thickness or not
make_eff_thk = True

# grid select
Gvel = GV
Gthk = GT
Gcon = GC

rgd_thk = True
rgd_con = True

# smoothing options
# smthvel = False
# smththk = False
# smthcon = False
smthvel = True
smththk = True
smthcon = True

#### old uniform
# velsmrad = 6
# Vsmth = lambda x: b.smooth2a(x,velsmrad)
# thksmrad = 4
# Tsmth = lambda x: b.smooth2a(x,thksmrad)
# consmrad = 1
# Csmth = lambda x: b.smooth2a(x,consmrad)

#### new vary smoothing
velsmthdist = 150e3
thksmthdist = 100e3
consmthdist = 60e3

VSmthObj = gs.geo_vary_smooth(Gvel,velsmthdist,in_mask = True,verbos=True)
TSmthObj = gs.geo_vary_smooth(Gthk,thksmthdist,in_mask = True,verbos=True)
CSmthObj = gs.geo_vary_smooth(Gcon,consmthdist,in_mask = True,verbos=True)


Vsmth = VSmthObj.smooth
Tsmth = TSmthObj.smooth
Csmth = CSmthObj.smooth

#### END USER OPTIONS
#### history fields
if hist_month_maps or hist_month_gates:
    budhist = b.budhist(Gvel,dlist,tw,drift_stats = calc_dstats,tot_budget=True,split_div=False,
#                         dumptw = relativedelta(days=1),budday0=dlist[1],  ### comment in/out for daily dumps
                        transport=True,save_input=save_input)

### we need to find most previous 04-15 or 10-15
# dcheck = dlist[0].replace(month = 4,day=15)
dcheckg = dlist[0].replace(month = 4,day=15)
dcheckm = dlist[0].replace(month = 10,day=15)
if dlist[0]>dcheckg and dlist[0]<dcheckm:
    ## then previous march 15
    budday0 = dcheckg
    print('GrowMelt budget, first dump on '+(budday0+relativedelta(months  = 6)).strftime('%Y-%m-%d'),file=pFile,flush=True)
elif dlist[0]>dcheckm:
    ## then this years Sept(Oct) 15
    budday0 = dcheckm
    print('GrowMelt budget, first dump on '+(budday0+relativedelta(months  = 6)).strftime('%Y-%m-%d'),file=pFile,flush=True)
else:
    #### if not previous years september 15
    budday0 = dcheckm-relativedelta(years=1)
    print('GrowMelt budget, first dump on '+(budday0+relativedelta(months  = 6)).strftime('%Y-%m-%d'),file=pFile,flush=True)
# if dlist[0]>dcheck:
#     ## then previous april 15
#     budday0 = dcheck
#     print('GrowMelt budget, first dump on '+(budday0+relativedelta(months  = 6)).strftime('%Y-%m-%d'),file=pFile,flush=True)
# else:
#     #### if not previous years october 15
#     budday0 = dcheck.replace(month = 10)-relativedelta(years=1)
#     print('GrowMelt budget, first dump on '+(budday0+relativedelta(months  = 6)).strftime('%Y-%m-%d'),file=pFile,flush=True)
hist_grow_melt = b.budhist(Gvel,dlist,tw,
                     dumptw = relativedelta(months  = 12),budday0=budday0,split_div = True,transport=True,tot_budget=True,save_input=save_input)
#### need every 6months growht/melt starts on 09-15
#### regrid onto NSIDC
histrgd = gs.Gs2Gs(Gvel,Gnsidc,vectors=True)
#### now we need all the gates
### make list of hist objects
if hist_month_maps or hist_month_gates:
    hist_list = [hist_grow_melt,budhist]
    ### hist_rgd_list is a list of list, one list for each hist_list
    hist_rdg_list = [[[histrgd,'GrowMelt_EASE_NH/',Gnsidc]],
                     [[False,'',Gvel]]]
else:
    hist_list = [hist_grow_melt]
    ### hist_rgd_list is a list of list, one list for each hist_list
    hist_rdg_list = [[[histrgd,'GrowMelt_EASE_NH/',Gnsidc]]]
    

## set GrowMelt gates
for gl in Gates_list:
    ## make regridder
    print('Building Season regrider: ',gl[1])
    Glrgd = gs.Gs2Gs(Gvel,gl[0],vectors = True)
    ## make hist object
    hist_rdg_list[0].append([Glrgd,'GrowMelt_EASE_NH/'+gl[1],gl[0]])
    
## set monthly gates
if hist_month_gates:
    for n,gl in enumerate(Gates_list):
        ### grid and mask info
        ## make regridder
#         print('Building Monthly regrider: ',gl[1])
#         Glrgd = gs.Gs2Gs(Gvel,gl[0],vectors = True)
#         ## make hist object
#         hist_rdg_list[1].append([Glrgd,gl[1],gl[0]])
        hist_rdg_list[1].append([hist_rdg_list[0][n+1][0],gl[1],gl[0]])


### to check all the the hist dirs to make they exist
for hist_rdg in hist_rdg_list:
    for hrl in hist_rdg:
        if type(hrl[0])==bool:
            pass
        else:
            histpath = os.path.dirname(spath+hrl[1])
            if not os.path.exists(histpath):
                # make if it doesn't
                os.makedirs(histpath)
                print('Creating directory: ',histpath,file=pFile,flush=True)
        #### save Gate grids in the dirs
        if 'Gate' in hrl[1]:
            ### hrl[2] is a grid_set
            hrl[2].save_grid(spath+hrl[1]+'Gate_gs.npz')
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

##### all the extra arrays for SIGMAS
hist = {}
hist_sn = {}
### normal terms
# hist_vars = ['vol','U','V','div','adv','int','res','dyn']
# hist_vars = ['volerr','Uerr','Verr']
hist_vars = []
### missing T terms (ice edge, thin ice)
hist_vars.extend(['div_t','adv_t','int_t','res_t','dyn_t'])
### new/old ice terms
hist_vars.extend(['new_ice','old_ice'])
### uncertainties
hist_vars.extend(['div_sig','adv_sig','int_sig','res_sig','dyn_sig'])
hist_vars.extend(['trans_x','trans_y','trans_x_sig','trans_y_sig'])
hist_vars.extend(['U_sig','V_sig','Vol_sig'])
for var in hist_vars:
    hist[var] = b.accumulator(Gvel.shape)
    hist_sn[var] = b.accumulator(Gvel.shape)

## F1 convariants
cor_dVdt = [] 
cor_grdV1 = [] 
cor_grdV2 = [] 
cor_divU1 = []
cor_divU2 = []
## F2 convariants
cor_sigu1grdVol1 = []
cor_sigu2grdVol2 = []
cor_sigdivu1divu2 = []
## F3 convariants
cor_adv1adv2 = []
cor_voldivU = []
## totals
cor_adv_div = []
cor_int_div = []
cor_int_adv = []

### historys for covariance testing
recent = {}
recent['nr'] = 12
# recent_vars = ['int','adv','div']
# #### some more for extra covariances in histories
# for var in recent_vars:
#     recent[var] = b.recent_array(Gvel.shape,recent['nr'])
#     for dr in range(1,recent['nr']):
#         recent[var+str(dr)] = []
### essentials for accumulating the covariances (int and adv)
cov_vars = ['int_c','adv_c','res_c','u_c','v_c','vol_c']
trans_vars = ['trans_x','trans_y']
for var in cov_vars+trans_vars:
    recent[var] = b.recent_array(Gvel.shape,recent['nr'])
    for dr in range(0,recent['nr']+1):
        recent[var].update(np.zeros(Gvel.shape)*Gvel.mask)
    
print("FIRST slice from "+ts.strftime('%Y%m%d'))

Uslices = np.empty([3,Gvel.m,Gvel.n])
Vslices = np.empty([3,Gvel.m,Gvel.n])
Tslices = np.empty([3,Gvel.m,Gvel.n])
Cslices = np.empty([3,Gvel.m,Gvel.n])
Vol_slices = np.empty([3,Gvel.m,Gvel.n])

# divergence advection are runnign means so make array
budadv = np.empty([3,Gvel.m,Gvel.n])
buddiv = np.empty([3,Gvel.m,Gvel.n])

Volerrslices = np.empty([3,Gvel.m,Gvel.n])
Uerrslices = np.empty([3,Gvel.m,Gvel.n])
Verrslices = np.empty([3,Gvel.m,Gvel.n])

### initalising the slices

Utemp,Vtemp = b.get_vels_array(InputV,ts,tw,verbos=load_verbos)
Verrt       =  b.get_err_array(InputV,ts,tw,verbos=load_verbos)
Utemp = Utemp*Gvel.mask
Vtemp = Vtemp*Gvel.mask
Verrt = Verrt.data*Gvel.mask
if smthvel:
    Uslices[0] = Vsmth(Utemp)
    Vslices[0] = Vsmth(Vtemp)
else:
    Uslices[0] = Utemp
    Vslices[0] = Vtemp
### Verr is non-dim, so we mult out (sorry not generalised)
### cehck if variance or sd, we want sd here, sqrt at a key point.
### for OSISAF the opposite is true, so check
Verrt = Verrt.data*Gvel.mask
Verrt = Verrt/np.hypot(Utemp,Vtemp) ## non-dim verrr
Uerrslices[0] = Utemp*Verrt
Verrslices[0] = Vtemp*Verrt

Ttemp= b.get_hi_array(InputT,ts,tw,verbos=load_verbos)
Terrt = b.get_err_array(InputT,ts,tw,verbos=load_verbos)
#### filling uncertainties
Tmask_p = Ttemp>Tlim
True_p = Tmask_p&np.isfinite(Terrt)
Fill_p = Tmask_p&np.isnan(Terrt)
### make the interpolator
xym = np.vstack( (GT.xpts[True_p], GT.ypts[True_p]) ).T
interp0 = NearestNDInterpolator( xym, Terrt[True_p])
### get values and locations (scale up by 50%)
Err_interp = interp0(GT.xpts[Fill_p], GT.ypts[Fill_p])*1.5
xF,yF = np.where(Fill_p)
Terrt[xF,yF] = Err_interp

if print_inputs: Tplot = copy.copy(Ttemp)
Ttemp[Ttemp <Tlim] = np.nan
Terrt[Ttemp <Tlim] = np.nan
Ttemp = Ttemp*Gthk.mask
if smththk:
    Ttemp = Tsmth(Ttemp)
if rgd_thk:
    # call the regrid object
    Ttemp = RGDthk.rg_array(Ttemp)
    Terrt = RGDthk.rg_array(Terrt,method = 'nearest')
Ttemp = Ttemp*Gvel.mask
Terrt = Terrt*Gvel.mask

Ctemp= b.get_aice_array(InputC,ts,tw,verbos=load_verbos)
if print_inputs: Cplot = copy.copy(Ctemp)
Ctemp[Ctemp <Clim] = np.nan
Ctemp = Ctemp*Gcon.mask
if smththk:
    Ctemp = Csmth(Ctemp)
if rgd_con:
    # call the regrid object
    Ctemp = RGDcon.rg_array(Ctemp)
Cslices[0] = Ctemp

if print_inputs: b.plot_inputs(Gvel,Gthk,Gcon,InputV,InputT,InputC,
                               Utemp,Vtemp,Tplot,Cplot,spath,ts)
    
### All the neccessary masking and extras
IsIce = Ctemp>0.15
Tmask = Ttemp>0.2#*np.isfinite(Tmask)
Vmask = np.isfinite(Vtemp)
Thin_mask = IsIce^Tmask
#### force more data
Ttemp[Thin_mask] = 0.2    



Tslices[0] = Ttemp
if make_eff_thk:
    Vol_slices[0] = Tslices[0]*Cslices[0]
else:
    Vol_slices[0] = Tslices[0]
Volerrslices[0] = Terrt

# #### we need to set all open ocean volume to 0.0,
# #### all land points to np.nan - use a grid mask
# if not cut_new_ice:
#     Vol_slices[0][np.isnan(Vol_slices[0])] = 0.0
#     Vol_slices[0] *=Gvel.mask

Uslices[1],Vslices[1] = Uslices[0],Vslices[0]
Uerrslices[1] = Uerrslices[0]
Verrslices[1] = Verrslices[0]
Tslices[1] = Tslices[0]
Cslices[1] = Cslices[0]
Vol_slices[1] = Vol_slices[0]
Volerrslices[1] = Volerrslices[0]

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
#     budhist.check_dump(n,d,hole_verbos)
    for hist_l in hist_list:
        hist_l.check_dump(n,d,hole_verbos)
    
    print("new slice from "+d.strftime('%Y%m%d'))
    # update final slice
    
    Utemp,Vtemp = b.get_vels_array(InputV,d,tw,verbos=load_verbos)
    Verrt        = b.get_err_array(InputV,d,tw,verbos=load_verbos)
    Verrt = Verrt.data*Gvel.mask
    Utemp = Utemp*Gvel.mask
    Vtemp = Vtemp*Gvel.mask
    if smthvel:
        Utemp = Vsmth(Utemp)
        Vtemp = Vsmth(Vtemp)
    Uslices[2] = Utemp
    Vslices[2] = Vtemp
    ### for OSISAF the opposite is true, so check
    Verrt = Verrt.data*Gvel.mask
    Verrt = Verrt/np.hypot(Utemp,Vtemp) ## non-dim verrr
    Uerrslices[2] = Utemp*Verrt
    Verrslices[2] = Vtemp*Verrt
    
    Ttemp = b.get_hi_array(InputT,d,tw,verbos=load_verbos)
    Terrt= b.get_err_array(InputT,d,tw,verbos=load_verbos)
    #### filling uncertainties
    Tmask_p = Ttemp>Tlim
    True_p = Tmask_p&np.isfinite(Terrt)
    Fill_p = Tmask_p&np.isnan(Terrt)
    ### make the interpolator
    xym = np.vstack( (GT.xpts[True_p], GT.ypts[True_p]) ).T
    interp0 = NearestNDInterpolator( xym, Terrt[True_p])
    ### get values and locations (scale up by 50%)
    Err_interp = interp0(GT.xpts[Fill_p], GT.ypts[Fill_p])*1.5
    xF,yF = np.where(Fill_p)
    Terrt[xF,yF] = Err_interp 
    Ttemp[Ttemp <Tlim] = np.nan
    Terrt[Ttemp <Tlim] = np.nan
    Ttemp = Ttemp*Gthk.mask
    Terrt = Terrt*Gthk.mask
    if smththk:
        Ttemp = Tsmth(Ttemp)
    if rgd_thk:
        # call the regrid object
        Ttemp = RGDthk.rg_array(Ttemp)
        Terrt = RGDthk.rg_array(Terrt,method = 'nearest')
    Ttemp = Ttemp*Gvel.mask
    Terrt = Terrt*Gvel.mask
    
    Ctemp= b.get_aice_array(InputC,d,tw,verbos=load_verbos)
    Ctemp[Ctemp <Clim] = np.nan
    Ctemp = Ctemp*Gcon.mask
    if smththk:
        Ctemp = Csmth(Ctemp)
    if rgd_con:
        # call the regrid object
        Ctemp = RGDcon.rg_array(Ctemp)
    Cslices[2] = Ctemp
    
    ### All the neccessary masking and extras
    IsIce = Ctemp>0.15
    Tmask = Ttemp>0.2#*np.isfinite(Tmask)
    Vmask = np.isfinite(Vtemp)
    Thin_mask = IsIce^Tmask
    Miss_mask = IsIce^Vmask
    #### force more data
    Ttemp[Thin_mask] = 0.2

    Tslices[2] = Ttemp


    #### actual budget
    if make_eff_thk:
        Vol_slices[2] = Tslices[2]*Cslices[2]
    else:
        Vol_slices[2] = Tslices[2]
    Volerrslices[2] = Terrt
        
#     #### we need to set all open ocea volume to 0.0,
#     #### all land points to np.nan - use a grid mask
#     if not cut_new_ice:
#         Vol_slices[2][np.isnan(Vol_slices[2])] = 0.0
#         Vol_slices[2] *=Gvel.mask

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
    
#     ### covariances F1
#     tstep = 2
#     x,y = Tslices[2-tstep],Tslices[2]
#     xymask = np.logical_and(x>0.2,y>0.2)
#     cor_dVdt.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     dstep = 2
#     x,y = Vol_slices[1,1:-1,dstep:], Vol_slices[1,1:-1,:-dstep]
#     xymask = np.logical_and(x>0.2,y>0.2)
#     cor_grdV1.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     x,y = Vol_slices[1,dstep:,1:-1], Vol_slices[1,:-dstep,1:-1]
#     xymask = np.logical_and(x>0.2,y>0.2)
#     cor_grdV2.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     x,y = Vslices[1,1:-1,dstep:], Vslices[1,1:-1,:-dstep]
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_divU1.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     x,y = Uslices[1,dstep:,1:-1], Uslices[1,:-dstep,1:-1]
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_divU2.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     ### covariances F2
#     grdVol1,grdVol2,grdVol1_sig,grdVol2_sig =  b.grd_V_sig(
#                     Vol_slices[1],Volerrslices[1],
#                     #### add for smooth covs
#                     cor_dvdx = 0.998,
#                     G=Gvel)
    
#     divV1,divV2,divV1_sig,divV2_sig =  b.div_U_comp_sig(
#                     Uslices[1],Vslices[1],
#                     Uerrslices[1],Verrslices[1],
#                     #### add for smooth covs
#                     cor_dUdx = 0.998,
#                     G=Gvel)
#     x = Uslices[1,1:-1,1:-1]
#     y = grdVol1
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_sigu1grdVol1.append(np.corrcoef(x[xymask],y[xymask])[0][1])
#     x = Vslices[1,1:-1,1:-1]
#     y = grdVol2
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_sigu2grdVol2.append(np.corrcoef(x[xymask],y[xymask])[0][1])
#     x = divV1[1:-1,1:-1]
#     y = divV2[1:-1,1:-1]
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_sigdivu1divu2.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     ### covariances F3
#     adv1,adv2,adv1_sig,adv2_sig = b.adv_comp_sig(
#                     Uslices[1],Vslices[1],Uerrslices[1],Verrslices[1],
#                     grdVol1,grdVol2,grdVol1_sig,grdVol2_sig,
#                     #### add for smooth covs
#                     G=Gvel)
#     divU, divU_sig = b.div_U_sig(
#                     divV1,divV2,divV1_sig,divV2_sig,
#                     #### add for smooth covs
#                     G=Gvel)

#     ### add: adv1adv2, voldivU
#     x = adv1
#     y = adv2
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_adv1adv2.append(np.corrcoef(x[xymask],y[xymask])[0][1])
#     x = Vol_slices[1,1:-1,1:-1]
#     y = divU
#     xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#     cor_voldivU.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
    inttemp,int_sig = b.intensification_sig(Vol_slices,Volerrslices,
                            #### add for smooth covs
                            ### PATHfinder and OSISAF
                            cor_Volt = 0.998,
                            nsecs = nsecs)
    advtemp,advsig = b.advection_sig(Uslices[1],Vslices[1],Vol_slices[1],
                            Uerrslices[1],Verrslices[1],Volerrslices[1],
                            #### add for smooth covs
                            ### PATHfinder
#                             cor_dvdx = 0.998,
#                             cor_sigUdivV=0.028,
#                             cor_adv1adv2 = -0.17,
                            ### OSISAF
                            cor_dvdx = 0.91,
                            cor_sigUdivV=0.028,
                            cor_adv1adv2 = -0.09,
                            G=Gvel)
    divtemp,divsig = b.divergence_sig(Uslices[1],Vslices[1],Vol_slices[1],
                            Uerrslices[1],Verrslices[1],Volerrslices[1],
                            #### add for smooth covs
                            ### PATHfinder
#                             cor_dUdx = 0.995,
#                             cor_divU12=-0.3,
                            ### OSISAF
                            cor_dUdx = 0.96,
                            cor_divU12=-0.56,
                            G=Gvel)
    int_adv_cov = 0.08
    res_sig2 = int_sig**2 + advsig**2 + divsig**2 \
             - 2*int_adv_cov*int_sig*advsig
    ressig = np.sqrt(res_sig2)
    #### EXTRA HISTORIES JUST IN CASE....
#     hist['volerr'].update(Volerrslices[1])
#     hist['Uerr'].update(Uerrslices[1])
#     hist['Verr'].update(Verrslices[1])
#     ### update recent arrays
#     recent['int'].update(inttemp)
#     recent['adv'].update(advtemp)
#     recent['div'].update(divtemp)
    ### correlate the recent arrays
#     for var in recent_vars:
#         for dn in range(1,np.minimum(n,recent['nr'])):
#             x = recent[var].data[recent['nr']-dn-1]
#             y = recent[var].data[-1]
#             xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
#             recent[var+str(dn)].append(np.corrcoef(x[xymask],y[xymask])[0][1])
    ### SANITY CHECK....
#     hist['div'].update(divtemp,mask = ~Tmask)
#     hist['adv'].update(advtemp,mask = ~Tmask)
#     hist['int'].update(inttemp,mask = ~Tmask)
#     hist['res'].update(inttemp-divtemp-advtemp,mask = ~Tmask)
#     hist['dyn'].update(divtemp+advtemp,mask = ~Tmask)
    
    #### estimated covariances for previous n_days values
#     cor_intint = [0.7,0.3,0.3,0.3,0.3]
#     cor_advadv = [0.6,0.2,0.2,0.2,0.2]
    cor_intint = 0.6
    cor_advadv = 0.7
    cor_resres = 0.8
    cor_uu = 0.6
    cor_volvol = 0.95
    #### sum up all the covariances over the n_days
    int_cov_sum = np.zeros(Gvel.shape)
    adv_cov_sum = np.zeros(Gvel.shape)
    res_cov_sum = np.zeros(Gvel.shape)
    u_cov_sum = np.zeros(Gvel.shape)
    v_cov_sum = np.zeros(Gvel.shape)
    vol_cov_sum = np.zeros(Gvel.shape)
#     for dr,dcov in zip(range(0,recent['nr']),cor_intint):
#         int_cov_sum += recent['int_c'].data[recent['nr']-dr-1]*dcov*int_sig
#     for dr,dcov in zip(range(0,recent['nr']),cor_advadv):
#         adv_cov_sum += recent['adv_c'].data[recent['nr']-dr-1]*dcov*advsig
    for dc in range(recent['nr']):
        int_cov_sum += recent['int_c'].data[recent['nr']-dc-1]*cor_intint**(dc+1)*int_sig
        adv_cov_sum += recent['adv_c'].data[recent['nr']-dc-1]*cor_advadv**(dc+1)*advsig
        res_cov_sum += recent['res_c'].data[recent['nr']-dc-1]*cor_resres**(dc+1)*ressig
        u_cov_sum += recent['u_c'].data[recent['nr']-dc-1]*cor_uu**(dc+1)*Uerrslices[1]
        v_cov_sum += recent['v_c'].data[recent['nr']-dc-1]*cor_uu**(dc+1)*Verrslices[1]
        vol_cov_sum += recent['vol_c'].data[recent['nr']-dc-1]*cor_volvol**(dc+1)*Volerrslices[1]
    
    int_adv_cov = 0.08
    int_sig_tot2= int_sig**2 + 2*int_cov_sum
    adv_sig_tot2= advsig**2 + 2*adv_cov_sum
    ### sqrting because we need int->adv cov
    int_sig_tot = np.sqrt(int_sig_tot2)
    adv_sig_tot = np.sqrt(adv_sig_tot2)
    div_sig_tot2= divsig**2 
    
#     hist['int_sig'].update(int_sig**2 + 2*cor_intint*int_sig*prev_int_sig)
#     hist['adv_sig'].update(advsig**2  + 2*cor_advadv*advsig*prev_adv_sig)
    hist['int_sig'].update(int_sig_tot2)
    hist['adv_sig'].update(adv_sig_tot2)
    hist['div_sig'].update(div_sig_tot2)
    hist_sn['int_sig'].update(int_sig_tot2)
    hist_sn['adv_sig'].update(adv_sig_tot2)
    hist_sn['div_sig'].update(div_sig_tot2)
#     res_sig2 = int_sig_tot2 + adv_sig_tot2 + div_sig_tot2 \
#              - 2*int_adv_cov*int_sig_tot*adv_sig_tot
    hist['res_sig'].update(res_sig2)
    hist_sn['res_sig'].update(res_sig2)
    dyn_sig2 = adv_sig_tot2 + div_sig_tot2 
    hist['dyn_sig'].update(dyn_sig2)
    hist['U_sig'].update(Uerrslices[1]**2 + 2*u_cov_sum)
    hist['V_sig'].update(Verrslices[1]**2 + 2*v_cov_sum)
    hist['Vol_sig'].update(Volerrslices[1]**2 + 2*vol_cov_sum)
    hist_sn['dyn_sig'].update(dyn_sig2)
    hist_sn['U_sig'].update(Uerrslices[1]**2 + 2*u_cov_sum)
    hist_sn['V_sig'].update(Verrslices[1]**2 + 2*v_cov_sum)
    hist_sn['Vol_sig'].update(Volerrslices[1]**2 + 2*vol_cov_sum)
    #### need to add all the covariances for int and adv
    recent['int_c'].update(int_sig)
    recent['adv_c'].update(advsig)
    recent['res_c'].update(ressig)
    recent['u_c'].update(Uerrslices[1])
    recent['v_c'].update(Verrslices[1])
    recent['vol_c'].update(Volerrslices[1])
    
    #### accumulate transport
    ### ['trans_x','trans_y','trans_x_sig','trans_y_sig']
    transx,transx_sig= b.mult_sig(Vol_slices[1],Uslices[1],Volerrslices[1],Uerrslices[1])
    transy,transy_sig= b.mult_sig(Vol_slices[1],Vslices[1],Volerrslices[1],Verrslices[1])
    hist['trans_x'].update(transx)
    hist['trans_y'].update(transy)
    hist_sn['trans_x'].update(transx)
    hist_sn['trans_y'].update(transy)
    recent['trans_x'].update(transx_sig)
    recent['trans_y'].update(transy_sig)
    ### transport covariance
#     cor_trans = [0.55,0.2,0.1,0.1,0.1]
    cor_trans = 0.55
    transx_cov_sum = np.zeros(Gvel.shape)
    transy_cov_sum = np.zeros(Gvel.shape)
#     for dr,dcov in zip(range(0,recent['nr']),cor_trans):
#         transx_cov_sum += recent['trans_x'].data[recent['nr']-dr-1]*dcov*transx_sig
#         transy_cov_sum += recent['trans_y'].data[recent['nr']-dr-1]*dcov*transy_sig
    for dc in range(recent['nr']):
        transx_cov_sum  += recent['trans_x'].data[recent['nr']-dc-1]*cor_trans**(dc+1)*transx_sig
        transy_cov_sum  += recent['trans_y'].data[recent['nr']-dc-1]*cor_trans**(dc+1)*transy_sig
    hist['trans_x_sig'].update(transx_sig**2 + 2*transx_cov_sum)
    hist['trans_y_sig'].update(transy_sig**2 + 2*transy_cov_sum)
    hist_sn['trans_x_sig'].update(transx_sig**2 + 2*transx_cov_sum)
    hist_sn['trans_y_sig'].update(transy_sig**2 + 2*transy_cov_sum)
    
    ### need to check if there is any covariance between:
    ### div and adv
    x = divtemp
    y = advtemp
    xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
    cor_adv_div.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    ### int and dyamics...
    x = inttemp
    y = divtemp
    xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
    cor_int_div.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    x = inttemp
    y = advtemp
    xymask = np.logical_and(np.isfinite(x),np.isfinite(y))
    cor_int_adv.append(np.corrcoef(x[xymask],y[xymask])[0][1])
    
#     budhist.accumulate_budget(d,buddif,budadv3,buddiv3,
#                               Tslices[1],Cslices[1],Uslices[1],Vslices[1],
#                               hole_verbos)
    
    ### ice edge terms
    hist['div_t'].update(buddiv3,mask = ~Thin_mask)
    hist['adv_t'].update(budadv3,mask = ~Thin_mask)
    hist['int_t'].update(buddif,mask = ~Thin_mask)
    hist['res_t'].update(buddif-buddiv3-budadv3,mask = ~Thin_mask)
    hist['dyn_t'].update(buddiv3+budadv3,mask = ~Thin_mask)
    hist_sn['div_t'].update(buddiv3,mask = ~Thin_mask)
    hist_sn['adv_t'].update(budadv3,mask = ~Thin_mask)
    hist_sn['int_t'].update(buddif,mask = ~Thin_mask)
    hist_sn['res_t'].update(buddif-buddiv3-budadv3,mask = ~Thin_mask)
    hist_sn['dyn_t'].update(buddiv3+budadv3,mask = ~Thin_mask)
    ### new/old ice
    hist['new_ice'].update(np.isnan(Vol_slices[0])*Vol_slices[1])
    hist['old_ice'].update(-1.0*np.isnan(Vol_slices[1])*Vol_slices[0])
    hist_sn['new_ice'].update(np.isnan(Vol_slices[0])*Vol_slices[1])
    hist_sn['old_ice'].update(-1.0*np.isnan(Vol_slices[1])*Vol_slices[0])
    for hist_l in hist_list:
        hist_l.accumulate_budget(d,buddif,budadv3,buddiv3,
                              Tslices[1],Cslices[1],Uslices[1],Vslices[1],mask=~Tmask,
                              hole_verbos=hole_verbos)
    
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
    
    Uerrslices[0] = Uerrslices[1]
    Uerrslices[1] = Uerrslices[2]
    Verrslices[0] = Verrslices[1]
    Verrslices[1] = Verrslices[2]
    
    Tslices[0] = Tslices[1]
    Tslices[1] = Tslices[2]
    Cslices[0] = Cslices[1]
    Cslices[1] = Cslices[2]
    
    Vol_slices[0] = Vol_slices[1]
    Vol_slices[1] = Vol_slices[2]
    
    Volerrslices[0] = Volerrslices[1]
    Volerrslices[1] = Volerrslices[2]
    
    budadv[0] = budadv[1]
    budadv[1] = budadv[2]
    
    buddiv[0] = buddiv[1]
    buddiv[1] = buddiv[2]
    
    # also need to dump on end of list
    # also need to dump when we reach a hole 
    for hist_l,hist_rgd in zip(hist_list,hist_rdg_list):
        if hist_l.hist_dump_now:
            if print_history:
                b.plot_budget(Gvel,hist_l.buddifhist,
                              hist_l.buddivhist,hist_l.buddivhist,
                              hist_l.budxvelhist/hist_l.hist_count,
                              hist_l.budyvelhist/hist_l.hist_count,
                              spath,d)
            if print_history_square:
                b.plot_budget_square(Gvel,hist_l.buddifhist,
                              hist_l.buddivhist,hist_l.buddivhist,
                              hist_l.budxvelhist/hist_l.hist_count,
                              hist_l.budyvelhist/hist_l.hist_count,
                              spath,d)
            ### optional second dump
#             hist.dump_nc(spath,d,Gvel,InputV,InputT,InputC,
#                            outrgd = [histrgd,'_EASE',GC],accu_reset=False)
            ### now loop the history regridders
            for nh,hrl in enumerate(hist_rgd):
                #### need to check for GrowMelt we add the extra conc and thick arrays
                if 'GrowMelt' in hrl[1]:
                    if 'Gate' in hrl[1]:
                        print_info = False
                        extra_arrays = [[var,hist_sn[var].total()*nsecs**2] for var in ['trans_x_sig','trans_y_sig']]
                    else:
                        extra_arrays = [
                            ['inst_thickness',Tslices[1]],
                            ['inst_thickness_sig',Volerrslices[1]],
                            ['inst_concentration',Cslices[1]],
                                        ] 
                        extra_arrays.extend([[var,hist_sn[var].total()*nsecs**2] for var in 
                                        ['div_sig','adv_sig','int_sig','res_sig','dyn_sig']]) #hist_vars if 'sig' in var]
                        extra_arrays.extend([[var,hist_sn[var].total()*nsecs] for var in 
                                         ['div_t','adv_t','int_t','res_t','dyn_t']])    # hist_vars if '_t' in var])
                        extra_arrays.extend([[var,hist_sn[var].total()] for var in 
                                            ['new_ice','old_ice']]) # hist_vars if '_ice' in var])
                        extra_arrays.extend([[var,hist_sn[var].mean()/hist_sn[var].count] for var in 
                                            ['U_sig','V_sig','Vol_sig']]) # hist_vars if '_ice' in var])
                        print_info = True
                elif 'Gate' in hrl[1]:
                    print_info = False
                    extra_arrays = [[var,hist[var].total()*nsecs**2] for var in ['trans_x_sig','trans_y_sig']]
                else:
                    print_info = True
                    #### all the new hist accumulators
                    ### extra_arrays = [[label,data]]
                    extra_arrays = [[var,hist[var].total()*nsecs**2] for var in 
                                    ['div_sig','adv_sig','int_sig','res_sig','dyn_sig']] #hist_vars if 'sig' in var]
                    extra_arrays.extend([[var,hist[var].total()*nsecs] for var in 
                                     ['div_t','adv_t','int_t','res_t','dyn_t']])    # hist_vars if '_t' in var])
                    extra_arrays.extend([[var,hist[var].total()] for var in 
                                        ['new_ice','old_ice']]) # hist_vars if '_ice' in var])
                    extra_arrays.extend([[var,hist[var].mean()/hist[var].count] for var in 
                                        ['U_sig','V_sig','Vol_sig']]) # hist_vars if '_ice' in var])
                ### only accu_reset on the last one
                ### add all the new acuumulators as extra arrays
                accu_reset = False
                if nh == len(hist_rdg)-1: accu_reset = True
                print('History time for '+hrl[1])
                hist_l.dump_nc(spath,d,hrl[2],InputV,InputT,InputC,
                             outrgd = hrl,extra_arrays = extra_arrays,
                             accu_reset=accu_reset, print_info = print_info)
                if accu_reset and 'GrowMelt' not in hrl[1]:
                    ### resetting the slices
                    print('Cleaning the uncertainty history',file=pFile,flush=True)
                    for var in hist_vars:
                        hist[var].clean()
                    for var in cov_vars+trans_vars:
                        recent[var].clean()
#                         for d in range(1,recent['nr']):
#                             recent[var+str(d)] = []
                elif accu_reset and 'GrowMelt' in hrl[1]:
                    print('Cleaning the uncertainty seasonal history',file=pFile,flush=True)
                    for var in hist_vars:
                        hist_sn[var].clean()

        if hist_l.hist_dump_now and calc_dstats:
            hist_l.dump_dstats_nc(spath,d,Gvel,InputV)
    #         hist.dump_pkl(spath,d,dstats = calc_dstats)

#### special for here: print to screen all the corellations
# with open("covs_to_copy.txt", "a") as pF:
#     print('Extra covariants to check, copy to a notebook, etc',file=pF,flush=True)
#     print('COVARIANT--COPY',file=pF,flush=True)
#     ## F1 convariants
#     print('cor_dVdt = ',file=pF,flush=True)
#     print(cor_dVdt,file=pF,flush=True)
#     print('cor_grdV1 = ',file=pF,flush=True)
#     print(cor_grdV1,file=pF,flush=True)
#     print('cor_grdV2 = ',file=pF,flush=True)
#     print(cor_grdV2,file=pF,flush=True)
#     print('cor_divU1 = ',file=pF,flush=True)
#     print(cor_divU1,file=pF,flush=True)
#     print('cor_divU2 = ',file=pF,flush=True)
#     print(cor_divU2,file=pF,flush=True)
#     ## F2 convariants
#     print('cor_sigu1grdVol1 = ',file=pF,flush=True)
#     print(cor_sigu1grdVol1,file=pF,flush=True)
#     print('cor_sigu2grdVol2 = ',file=pF,flush=True)
#     print(cor_sigu2grdVol2,file=pF,flush=True)
#     print('cor_sigdivu1divu2 = ',file=pF,flush=True)
#     print(cor_sigdivu1divu2,file=pF,flush=True)
#     ## F3 convariants
#     print('cor_adv1adv2 = ',file=pF,flush=True)
#     print(cor_adv1adv2,file=pF,flush=True)
#     print('cor_voldivU = ',file=pF,flush=True)
#     print(cor_voldivU,file=pF,flush=True)
#     ## totals
#     print('cor_adv_div = ',file=pF,flush=True)
#     print(cor_adv_div,file=pF,flush=True)
#     print('cor_int_div = ',file=pF,flush=True)
#     print(cor_int_div,file=pF,flush=True)
#     print('cor_int_adv = ',file=pF,flush=True)
#     print(cor_int_adv,file=pF,flush=True)