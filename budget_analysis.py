from sys import path
path.insert(0, '/Users/h/Github/geo_data_group/')
import data_year as dy
import grid_set as gs
import numpy as np
import copy
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import datetime as dt
from dateutil.relativedelta import relativedelta
import copy
import gc




budget_load_default = ['intensification','advection','divergence',
                       'residual','dynamics','ice_drift_x','ice_drift_y']
class budget_result(dy.data_year):
    def __init__(self, directory,verbos = False,var_load = budget_load_default,regrid = None,year_limit = None,seasonal = False,yearly = False,filter_count = None,daily=False,ignore_file = []):
        self.path = directory
        files = glob.glob(directory+'*.nc')
        files.sort()
        #### find all the files
        self.n_t = 0
        self.fdates = []
        dates_list = []
        self.new_files = False
        self.all_files = files
        self.load_files = []
        for file in files: 
            if file.split('/')[-1] in ignore_file: continue
            if verbos: print(file) 
            try:
                date_str = file.split('budfields_')[1]
                if len(date_str)>12:
                    if seasonal:
                        date_str = date_str.split('--')[1]
                    else:
                        date_str = date_str.split('--')[0]
                    self.new_files = True
                date_str = date_str.split('.nc')[0]
                date_use = dt.datetime.strptime(date_str,'%Y%m%d')
                if verbos: print(date_str) 
            except ValueError:
                pass
            except IndexError:
                pass
            else:
                #### here is the year limit option
                if (year_limit is None) or (date_use.year>=year_limit[0] and date_use.year<=year_limit[1]):
                    self.n_t += 1
                    self.fdates.append(date_use)
                    if date_use.day > 1 and not daily:
                        date_use = date_use.replace(day = 1)
    #                     if date_use.month <6:
    #                         date_use+=relativedelta(months=1)
                    dates_list.append(date_use)
                    self.load_files.append(file)
#                 else: print(date_use.year,year_limit)
        self.fdates.sort()
        dates_list.sort()
        #### create dates array
        #### if we are regridding, set dimensions from reggridder
        if type(regrid) ==  gs.Gs2Gs:
            m,n = regrid.mesh_new[0].shape
            rgA = regrid.rg_array
            rgV = regrid.rg_vecs
        else: ## use read in
            rgA = lambda x: x
            rgV = lambda x,y: (x,y)
            #### open first file 
            if len(self.fdates) == 0:
                print('No dates found, check directory structure and year_limits')
                print('Directory = '+directory)
                return None
            d = self.fdates[0]
            if self.new_files and seasonal:
                #### go for the second_date
                f_name = glob.glob(self.path + 'budfields_*' + d.strftime('%Y%m%d') + '.nc')[0]
            elif self.new_files:
                #### go for the first date
                f_name = glob.glob(self.path + 'budfields_' + d.strftime('%Y%m%d') + '*.nc')[0]
            else:
                f_name = self.path + 'budfields_' + d.strftime('%Y%m%d') + '.nc'
            #### get dims
            if self.new_files:
                m,n = Dataset(f_name).variables['lon'].shape[:]
            else:
                m,n = Dataset(f_name).variables['lon'].shape[1:]
        #### create attributes
        for vl in var_load:
            if verbos: print(vl)
            setattr(self,vl,np.empty([self.n_t,m,n]))
        #### load the data from
        self.days_per_slice = []
        for n,d in enumerate(self.fdates):
            if self.new_files:
                if seasonal:
                    f_name = glob.glob(self.path + 'budfields_*' + d.strftime('--%Y%m%d') + '.nc')[0]
                else:
                    f_name = glob.glob(self.path + 'budfields_' + d.strftime('%Y%m%d') + '*.nc')[0]
            else:
                f_name = self.path + 'budfields_' + d.strftime('%Y%m%d') + '.nc'
            f_nc = Dataset(f_name)
            self.days_per_slice.append(f_nc.getncattr('Budget time slices'))
            for vl in var_load:
#                 if vl == 'ice_drift_x':
                if 'ice_drift' in vl: 
                    ### careful for vectors
                    if 'ice_drift_y' in var_load and 'ice_drift_x' in var_load:
                        ### we want both so we have to check for regridding
                        ### set x and y
                        if vl == 'ice_drift_x':
                            ### regrid when considering x
        #                     x = getattr(self,vl)
                            x = getattr(self,'ice_drift_x')
                            y = getattr(self,'ice_drift_y')
                            xr,yr = rgV(np.squeeze(f_nc.variables['ice_drift_x']),
                                        np.squeeze(f_nc.variables['ice_drift_y']))
                            x[n] = xr
                            y[n] = yr
                        elif vl == 'ice_drift_y':
                            ### already done previously
                            continue
                    else:
                        ### we only want one so do it normally
                        x = getattr(self,vl)
                        x[n] = rgA(np.squeeze(f_nc.variables[vl]))
                else:
                    x = getattr(self,vl)
                    x[n] = rgA(np.squeeze(f_nc.variables[vl]))
        #### set one data
        if len(var_load)>0:
            data = copy.copy(getattr(self,var_load[0]))
        else: ### empty class
            data = np.zeros([len(self.days_per_slice)])
#         data = getattr(self,var_load[0])
        #### initialise the data_year using data
        if seasonal:
            super().__init__(data,dates_list,periods=2)
        elif daily:
            super().__init__(data,dates_list,periods=366)
        else:
            super().__init__(data,dates_list)
        ### set the variables loaded
        self.vars = var_load
        
    def set_var(self,varname):
        """
        use a variable name to set as 'data'
        this allows all the dy methods to work correctly
        """
#         self.data = copy.copy(getattr(self,varname))
#         x = getattr(self,varname)
        self.data[:] = getattr(self,varname)
#         data = getattr(self,varname)
#         gc.collect()
    
    def __sub__(self,other):
        #### check shape
        if (self.n_t != other.n_t or 
            self.m   != other.m or 
            self.n   != other.n):
            print("These budget ain't the same shape so they can be subtracted")
            print("B1.shape = ",self.data.shape,", B2.shape = ",other.data.shape)
            return None
        #### check varaibles
        check_list = [t in other.vars for t in self.vars ]
        if np.sum(check_list) ==0:
            print("These budget don't have any data in common - check the .vars")
            return None
        if np.sum(check_list) != len(check_list):
            print("These budget don't have all data in common")
            subtract_list = [t for t,tt in zip(self.vars,check_list) if tt]
            print("only subtracting:" +
                ', '.join(subtract_list))
            del_list = [t for t,tt in zip(self.vars,check_list) if not tt]
        else:
            subtract_list = self.vars
            del_list= []
        out_bud = copy.copy(self)
        out_bud.vars = subtract_list
        for v in subtract_list:
#             print(v)
            x  = getattr(out_bud,v)
            y1 = getattr(self   ,v)
            y2 = getattr(other  ,v)
            x = y1 - y2
            setattr(out_bud,v,x)
        for v in del_list:
            delattr(out_bud,v)
        return out_bud
            
### copy class for redundancy           
#class budget_result_monthly(budget_result)       
        
    def mask_by_count(self,count_min):
        """
        delete entries from yrpd if there are too few days in the days_per_slice array
        
        Parameters
        ---------
        count_min: int, the number of days_per_slice we want more than.
            for example for a monthly budget, set this to 25 will make it so we ignore all cases with fewer than 25 days (say weird beginning and end half months)
            
        Alternatives
        ------------
        
        I might write a sacle by count method, to boost up under represented months,
        this could well be a bit dodgy though.
        """
        for nt,dps in enumerate(self.days_per_slice):
            if dps < count_min:
                ### find the yrpd entry
                y,p = np.unravel_index(np.abs(self.yrpd.data - nt).argmin(),
                                       self.yrpd.shape)
                self.yrpd.mask[y,p] = True
            

list_months_titles=['October','November','December','January','February','March','April']
list_ylabels=['Intensification m/month','Advection m/month','Divergence m/month','Residual m/month']
def plot_budget_season(bud_month,mstart,nmonths,year,GS,mask=False,
                       save_dir = False,show_plot = True,plot_type='.png',plot_label = '',plevels = (-0.6,0.8,0.2),
                       depth_nn=False,drift_vectors = False,rv=(15,0.2,0.01,1),rvGplot=None,dbug=False,extend_div=False):
    bs0 = plevels[0]
    if bud_month.new_files:
        pscale = 1.0
    else:
        pscale = 60*60*24
    if extend_div:
        bud_list = ['advection','divergence','positive_divergence','negative_divergence']
        bud_label_list = ['Advection m/month','Divergence m/month',
                          'Divergence (+ve part) m/month','Divergence (-ve part) m/month']
    else:
        bud_list = budget_load_default[:4]
        bud_label_list = list_ylabels
    bs1 = plevels[1]  - plevels[2]
    blevels = np.arange(plevels[0],plevels[1],plevels[2])
    mshift = 12-mstart
    if drift_vectors:
        rm = int(GS.m/rv[0])
        rn = int(GS.n/rv[0])
        ra = np.sqrt(rm+rn)
        ra=ra*rv[1]
        ### temp vector variables
#         vec_plot = np.empty([2,nmonths,GS.m,GS.n])
        vec_plot0= np.empty([2,nmonths,GS.m,GS.n])
        vec_plot = np.empty([2,nmonths,rv[0],rv[0]])
        budstr = 'ice_drift_x'
        bud_month.set_var(budstr)
        for mn in range(nmonths):
            year_u = year
            m_u = np.mod(mn-mshift,12)
#             if m_u > 6: year_u = [y-1 for y in year_u]
            if m_u >= mstart: year_u = [y-1 for y in year_u]
            vec_plot0[0,mn,:,:] = bud_month.clim_map(periods=[m_u],year_set=year_u)
        budstr = 'ice_drift_y'
        bud_month.set_var(budstr)
        for mn in range(nmonths):
            year_u = year
            m_u = np.mod(mn-mshift,12)
            if m_u >= mstart: year_u = [y-1 for y in year_u]
            vec_plot0[1,mn,:,:] = bud_month.clim_map(periods=[m_u],year_set=year_u)
        for mn in range(nmonths):
            vec_plot[0,mn],vec_plot[1,mn]=rvGplot.rg_vecs(vec_plot0[0,mn],vec_plot0[1,mn])
    fig = plt.figure(figsize=(3*nmonths,12))
    for bn, budstr in enumerate(bud_list):
        bud_month.set_var(budstr)
        for mn in range(nmonths):
            pn = mn + bn*nmonths  + 1
            ax = fig.add_subplot(4,nmonths,pn,projection=GS.ccrs)
            ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
#             mt.drawparallels(np.arange(60,90,5),linewidth=0.5, dashes=[1,5])
#             mt.drawmeridians(np.arange(0,360,10),linewidth=0.5, dashes=[1,5])
            if type(depth_nn) is not bool:
                llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
                s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
            
            #### get the array of interest - return logical if it exists
            year_u = year
#             m_u = np.mod(mn-3,12)
#             if m_u > 6: year_u = [y-1 for y in year_u]
            m_u = np.mod(mn-mshift,12)
            if m_u >= mstart: year_u = [y-1 for y in year_u]
            try:
                if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                    budget = bud_month.clim_map(periods=[m_u],year_set=year_u,mask=mask)
                    
                    budget[budget==0] = np.nan
                    s = ax.contourf(GS.xpts, GS.ypts, pscale*budget, vmin=bs0, 
                                    vmax=bs1,levels=blevels,cmap=plt.cm.RdYlBu, extend='both')
                else: 
                    print('Skipping year from: ',year_u)
#                 mt.pcolormesh(GS.xptp, GS.yptp, 86400.*12.*budget, vmin=-6, vmax=6,cmap=plt.cm.RdYlBu, rasterized=True)
            except TypeError:
                pass
            except ValueError:
                pass
            try:
#                 mt.quiver(B.xpts[::rm,::rn],B.ypts[::rm,::rn],ur[::rm,::rn],vr[::rm,::rn],scale 
                if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                    ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                              vec_plot[0,mn],vec_plot[1,mn],
                              scale = ra,width=rv[2],alpha=rv[3])
            except TypeError:
                pass
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.COASTLINE)
#             if (mn==0):
# #                 ax.set_ylabel(list_ylabels[bn], rotation=0, size='large')
#                 ax.set_ylabel('test')
            #### top row print dates
            if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                if (bn==0 and year[0]==year[1]):
                    ax.set_title(bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%B %Y'))
                elif (bn==0):
    #                 print(bud_month.yrpd[year_u[0],m_u],bud_month.yrpd[year_u[1],m_u])
    #                 print(year_u[0],m_u,year_u[1],m_u)
                    ### for ranges sometimes the first entry is missing
                    y_now = year_u[0]
                    ## so print next year if so
                    if bud_month.yrpd.mask[year_u[0],m_u] == True: y_now+=1
                    ax.set_title(bud_month.print_date(bud_month.yrpd[y_now,m_u],'%B %Y - ')+
                                        bud_month.print_date(bud_month.yrpd[year_u[1],m_u],'%Y'))


            if (mn == nmonths-1):
#                 divider = make_axes_locatable(ax)
#                 cax = divider.append_axes('right', size='5%', pad=0.05)
#                 plt.colorbar(s, cax=ax, pad=0.1)
                cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
                cax.ax.yaxis.set_ticks_position('left')
                cax.set_label(bud_label_list[bn])
    
#             else:
#                 divider = make_axes_locatable(ax)
#                 cax = divider.append_axes('right', size='5%', pad=0.05)
#                 cax.axis('off')
            if dbug:
                print('dbug: stopping early')
                plt.show()
                return False
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02)
    if type(save_dir) == str: 
        if (year[0]==year[1]):
            print_date = bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%Y')
        else:
            print_date = (bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%Y-') + 
                          bud_month.print_date(bud_month.yrpd[year_u[1],m_u],'%Y'))
        fig.savefig(save_dir+plot_label+print_date+plot_type,bbox_inches='tight')
#     plt.subplots_adjust(wspace=0.05, hspace=0.05)
    if show_plot: plt.show()

clim_labels = ['Total Intensification m','Total Advection m',
               'Total Divergence m','Total Residual m']
def plot_budget_clim(bud_month,mstart,nmonths,years,GS,mask=False,
                     save_dir = False,show_plot = True,plot_type='.png',plot_label = '',
                     pdlevels = (-0.5,0.6,0.1),plevels = (-1.6,2.0,0.4),
                     depth_nn=False,drift_vectors = False,vector_anomaly=False,
                     rv=(15,0.2,0.01),rvGplot=None,dbug=False,
                     var_list = budget_load_default[:4],label_list=clim_labels,
                     year_list_labels = ['A','B','C','D']):
    if bud_month.new_files:
        pscale = 1.0
    else:
        pscale = 60*60*24
    bs0 = plevels[0]
    bs1 = plevels[1]  - plevels[2]
    blevels = np.arange(plevels[0],plevels[1],plevels[2])
    ds0 = pdlevels[0]
    ds1 = pdlevels[1]  - pdlevels[2]
    dlevels = np.arange(pdlevels[0],pdlevels[1],pdlevels[2])
    nrows = 1+len(years)
    ncols = len(var_list)
    fig = plt.figure(figsize=(3*ncols,3*nrows))
    #### top line is the climatology
    periods = [np.mod(m,12) for m in range(mstart,mstart+nmonths)]
    periods.sort()
    if drift_vectors:
        rm = int(GS.m/rv[0])
        rn = int(GS.n/rv[0])
        ra = np.sqrt(rm+rn)
        ra=ra*rv[1]
        ### temp vector variables
#         vec_plot = np.empty([2,nmonths,GS.m,GS.n])
        vec_plot0= np.empty([2,nrows,GS.m,GS.n])
        vec_plot = np.empty([2,nrows,rv[0],rv[0]])
        budstr = 'ice_drift_x'
        bud_month.set_var(budstr)
        vec_plot0[0,0,:,:] = bud_month.clim_map(periods=periods,mask=mask)
        for d_p,y in enumerate(years):
            if type(y) == list:
                idx = []
                for yu in y:
                    idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if yu+(m>11)<bud_month.nyrs])
            else:
                idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                          if y+(m>11)<bud_month.nyrs]
            idx = [i for i in idx if not np.ma.core.is_masked(i)]
            if vector_anomaly:
                vec_plot0[0,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) - vec_plot0[0,0,:,:]
            else:
                vec_plot0[0,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) 
        budstr = 'ice_drift_y'
        bud_month.set_var(budstr)
        vec_plot0[1,0,:,:] = bud_month.clim_map(periods=periods,mask=mask)
        d_p = 1
        for d_p,y in enumerate(years):
            if type(y) == list:
                idx = []
                for yu in y:
                    idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if yu+(m>11)<bud_month.nyrs])
            else:
                idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if y+(m>11)<bud_month.nyrs]
            idx = [i for i in idx if not np.ma.core.is_masked(i)]
            if vector_anomaly:
                vec_plot0[1,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) - vec_plot0[1,0,:,:]
            else:
                vec_plot0[1,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) 
            d_p +=1
        
        for nr in range(nrows):
            vec_plot[0,nr],vec_plot[1,nr]=rvGplot.rg_vecs(vec_plot0[0,nr],vec_plot0[1,nr])
                
        
    for bn, budstr in enumerate(var_list):
        bud_month.set_var(budstr)
        ### first row is clim map
        ax = fig.add_subplot(nrows,ncols,bn+1,projection=GS.ccrs)
        ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        if type(depth_nn) is not bool:
            llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
            s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05) 
        try:
            #### top row is the mean over the whole bud_month
            #### units here will be m/s*days in month
            budget0 = bud_month.clim_map(periods=periods,method = 'ysum',calc_mask=mask,mask=mask)
            budget0[budget0==0] = np.nan
            s = ax.contourf(GS.xpts, GS.ypts, pscale*budget0, vmin=bs0, vmax=bs1,levels=blevels,cmap=plt.cm.RdYlBu, extend='both')
            if bn == ncols-1:
                cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
                cax.ax.yaxis.set_ticks_position('left')
                cax.set_label('Climatology')
        except TypeError:
            pass
        if drift_vectors:
            try:
                ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                          vec_plot[0,0],vec_plot[1,0],
                          scale = ra,width=rv[2],alpha=rv[3])
            except TypeError:
                pass
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        
        #### titles and colorbars here
        ax.set_title(label_list[bn])
            
        for nr,y in enumerate(years):
            pn = bn + ncols*(nr+1) + 1 
            if type(y) == list:
                idx = []
                for yu in y:
                    idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if yu+(m>11)<bud_month.nyrs])
                ysum = len(y)
            else:
                idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if y+(m>11)<bud_month.nyrs]
                ysum = 1
            idx = [i for i in idx if not np.ma.core.is_masked(i)]
            ax = fig.add_subplot(nrows,ncols,pn,projection=GS.ccrs)
            ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
            if type(depth_nn) is not bool:
                llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
                s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
            
            try:
                if mask:
                    budget = np.nansum([bud_month[i] for i in idx], axis = 0)/ysum - budget0
                else:
                    budget = np.nansum([bud_month.data[i] for i in idx], axis = 0)/ysum - budget0
                budget[budget==0] = np.nan
                s = ax.contourf(GS.xpts, GS.ypts, pscale*budget, vmin=ds0, vmax=ds1,levels=dlevels,cmap=plt.cm.RdYlBu, extend='both')
                if bn == ncols-1:
                    cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
                    cax.ax.yaxis.set_ticks_position('left')
#                     cax.ax.set_yticklabels(['-2','','-1','','0','','1','','2'])
                    d0 = bud_month.dates[idx[0]]
                    d1 = bud_month.dates[idx[-1]]
                    if type(y) == list:
                        clabel = year_list_labels[nr]
                    elif d0.year == d1.year:
                        clabel = d0.strftime('%Y')+' diff.'
                    else:
                        clabel = d0.strftime('%Y/')+d1.strftime('%Y')+' diff.'
                    cax.set_label(clabel)
            except TypeError:
                pass
            if drift_vectors:
                ra_use = copy.copy(ra)
                if vector_anomaly:
                    ra_use *= 0.5
                try:
                    ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                              vec_plot[0,nr+1],vec_plot[1,nr+1],
                              scale = ra_use,width=rv[2],alpha=rv[3])
                except TypeError:
                    pass
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.COASTLINE)
            
            if dbug:
                print('dbug: stopping early')
                plt.show()
                return False
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02)
    if type(save_dir) == str: 
        # list the years
        if type(years[0]) == list:
            year_list = [str(bud_month.dates[bud_month.yrpd[y,mstart]].year) for y in years[0]]
        else:
            year_list = [str(bud_month.dates[bud_month.yrpd[y,mstart]].year) for y in years]
        fig.savefig(save_dir+plot_label+'-'.join(year_list)+'_climatology'+plot_type,bbox_inches='tight')
#     plt.subplots_adjust(wspace=0.05, hspace=0.05)
    if show_plot: plt.show()


budget_state = ['thickness','concentration']
label_state = ['Thickness','Concentration']
def plot_state_season(bud_month,mstart,nmonths,year,GS,save_dir = False,show_plot = True,plot_type='.png',plot_label = '',
            hlevels = np.arange(0,4,0.5),
            clevels = np.array([0,0.5,0.75,0.9,0.95,1.0]),
            vlevels = np.arange(0,0.22,0.02),
                      diffscale=False,
                              depth_nn=False,drift_vectors = False,rv=(15,0.2,0.01,1),rvGplot=None,dbug=False):
    mshift = 12-mstart
    nrows = 2
    if True:
        nrows = 3
        rm = int(GS.m/rv[0])
        rn = int(GS.n/rv[0])
        ra = np.sqrt(rm+rn)
        ra=ra*rv[1]
        ### temp vector variables
#         vec_plot = np.empty([2,nmonths,GS.m,GS.n])
        vec_plot0= np.empty([2,nmonths,GS.m,GS.n])
        vec_plot = np.empty([2,nmonths,rv[0],rv[0]])
        budstr = 'ice_drift_x'
        bud_month.set_var(budstr)
        for mn in range(nmonths):
            year_u = year
            m_u = np.mod(mn-mshift,12)
#             if m_u > 6: year_u = [y-1 for y in year_u]
            if m_u >= mstart: year_u = [y-1 for y in year_u]
            vec_plot0[0,mn,:,:] = bud_month.clim_map(periods=[m_u],year_set=year_u)
        budstr = 'ice_drift_y'
        bud_month.set_var(budstr)
        for mn in range(nmonths):
            year_u = year
            m_u = np.mod(mn-mshift,12)
            if m_u >= mstart: year_u = [y-1 for y in year_u]
            vec_plot0[1,mn,:,:] = bud_month.clim_map(periods=[m_u],year_set=year_u)
        for mn in range(nmonths):
            vec_plot[0,mn],vec_plot[1,mn]=rvGplot.rg_vecs(vec_plot0[0,mn],vec_plot0[1,mn])
    fig = plt.figure(figsize=(3*nmonths,3*nrows))
    for bn, budstr in enumerate(budget_state[:2]):
        bud_month.set_var(budstr)
        if budstr == 'thickness':
            llevels = hlevels
            vm = [llevels[0],llevels[-1]]
            if diffscale:
                cmap=plt.cm.RdYlBu
            else:
                cmap = plt.cm.plasma
        elif budstr == 'concentration':
            llevels = clevels
            vm = [llevels[0],llevels[-1]]
            if diffscale:
                cmap=plt.cm.RdYlBu
            else:
                cmap = plt.cm.YlGnBu 
        else:
            llevels = np.arange(-1,1,2,0.2)
            vm = [llevels[0],llevels[-1]]
            cmap=plt.cm.RdYlBu
            
        for mn in range(nmonths):
            pn = mn + bn*nmonths  + 1
            ax = fig.add_subplot(nrows,nmonths,pn,projection=GS.ccrs)
            ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
            if type(depth_nn) is not bool:
                llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
                s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
            
            year_u = year
            m_u = np.mod(mn-mshift,12)
            if m_u >= mstart: year_u = [y-1 for y in year_u]
            try:
                if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                    budget = bud_month.clim_map(periods=[m_u],year_set=year_u,method = 'ysum')
                    budget[budget==0] = np.nan
#                     llevels = np.arange(-6,6,1)
                    s = ax.contourf(GS.xpts, GS.ypts, budget, vmin=vm[0], vmax=vm[1], levels=llevels,
                                    cmap=cmap, extend='both')
                else: 
                    print('Skipping year from: ',year_u)
            except TypeError:
                pass
            try:
                if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                    ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                              vec_plot[0,mn],vec_plot[1,mn],
                              scale = ra,width=rv[2],alpha=rv[3])
            except TypeError:
                pass
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.COASTLINE)
#             if (mn==0):
# #                 ax.set_ylabel(list_ylabels[bn], rotation=0, size='large')
#                 ax.set_ylabel('test')

            if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                if (bn==0 and year[0]==year[1]):
                    ax.set_title(bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%B %Y'))
                elif (bn==0):
                    ax.set_title(bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%B %Y - ')+
                                        bud_month.print_date(bud_month.yrpd[year_u[1],m_u],'%Y'))


            if (mn == nmonths-1):
#                 divider = make_axes_locatable(ax)
#                 cax = divider.append_axes('right', size='5%', pad=0.05)
#                 plt.colorbar(s, cax=ax, pad=0.1)
                cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
                cax.ax.yaxis.set_ticks_position('left')
#                 cax.ax.set_yticklabels(llevels,rotation=90)
                cax.set_label(label_state[bn])
    
#             else:
#                 divider = make_axes_locatable(ax)
#                 cax = divider.append_axes('right', size='5%', pad=0.05)
#                 cax.axis('off')
            if dbug:
                print('dbug: stopping early')
                plt.show()
                return False
    ### use the vectors to make the bottom plots
    for mn in range(nmonths):
        if diffscale:
            cmap=plt.cm.RdYlBu
        else:
            cmap=plt.cm.viridis
        bn = 2
        pn = mn + bn*nmonths  + 1
        ax = fig.add_subplot(nrows,nmonths,pn,projection=GS.ccrs)
        ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        if type(depth_nn) is not bool:
            llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
            s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)

        year_u = year
        m_u = np.mod(mn-mshift,12)
        if m_u >= mstart: year_u = [y-1 for y in year_u]
        try:
            if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                budget = np.hypot(vec_plot0[0,mn],vec_plot0[1,mn])
                budget[budget==0] = np.nan
                llevels = vlevels
                s = ax.contourf(GS.xpts, GS.ypts, budget, vmin=llevels[0], vmax=llevels[-1], levels=llevels,
                                cmap=cmap, extend='both')
            else: 
                print('Skipping year from: ',year_u)
        except TypeError:
            pass
        try:
            if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
                ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                          vec_plot[0,mn],vec_plot[1,mn],
                          scale = ra,width=rv[2],alpha=rv[3])
        except TypeError:
            pass
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
#             if (mn==0):
# #                 ax.set_ylabel(list_ylabels[bn], rotation=0, size='large')
#                 ax.set_ylabel('test')

        if year_u[0] != year_u[1] or type(bud_month.yrpd[year_u[0],m_u]) == np.int64:
            if (bn==0 and year[0]==year[1]):
                ax.set_title(bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%B %Y'))
            elif (bn==0):
                ax.set_title(bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%B %Y - ')+
                                    bud_month.print_date(bud_month.yrpd[year_u[1],m_u],'%Y'))


        if (mn == nmonths-1):
#                 divider = make_axes_locatable(ax)
#                 cax = divider.append_axes('right', size='5%', pad=0.05)
#                 plt.colorbar(s, cax=ax, pad=0.1)
            cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
            cax.ax.yaxis.set_ticks_position('left')
#             cax.ax.set_yticklabels(llevels,rotation=90)
            cax.set_label('Drift Speed')
    
#             else:
#                 divider = make_axes_locatable(ax)
#                 cax = divider.append_axes('right', size='5%', pad=0.05)
#                 cax.axis('off')
    if type(save_dir) == str: 
        if (year[0]==year[1]):
            print_date = bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%Y')
        else:
            print_date = (bud_month.print_date(bud_month.yrpd[year_u[0],m_u],'%Y-') + 
                          bud_month.print_date(bud_month.yrpd[year_u[1],m_u],'%Y'))
        fig.savefig(save_dir+plot_label+print_date+plot_type,bbox_inches='tight')
#     plt.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02)
    if show_plot: plt.show()



def plot_state_clim(bud_month,mstart,nmonths,years,GS,hscale,vscale,save_dir = False,show_plot = True,plot_type='.png',plot_label = '',
                              depth_nn=False,drift_vectors = False,vector_anomaly=False,rv=(15,0.2,0.01),rvGplot=None,dbug=False,
                     year_list_labels = ['A','B','C','D']):
    nrows = 1+len(years)
    fig = plt.figure(figsize=(9,3*nrows))
    #### top line is the climatology
    periods = [np.mod(m,12) for m in range(mstart,mstart+nmonths)]
    periods.sort()
    if True:
        rm = int(GS.m/rv[0])
        rn = int(GS.n/rv[0])
        ra = np.sqrt(rm+rn)
        ra=ra*rv[1]
        ### temp vector variables
#         vec_plot = np.empty([2,nmonths,GS.m,GS.n])
        vec_plot0= np.empty([2,nrows,GS.m,GS.n])
        vec_plot = np.empty([2,nrows,rv[0],rv[0]])
        budstr = 'ice_drift_x'
        bud_month.set_var(budstr)
        vec_plot0[0,0,:,:] = bud_month.clim_map(periods=periods)
        for d_p,y in enumerate(years):
            if type(y) == list:
                idx = []
                for yu in y:
                    idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if yu+(m>11)<bud_month.nyrs])
            else:
                idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if y+(m>11)<bud_month.nyrs]
            idx = [i for i in idx if not np.ma.core.is_masked(i)]
            if vector_anomaly:
                vec_plot0[0,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) - vec_plot0[0,0,:,:]
            else:
                vec_plot0[0,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) 
        budstr = 'ice_drift_y'
        bud_month.set_var(budstr)
        vec_plot0[1,0,:,:] = bud_month.clim_map(periods=periods)
        d_p = 1
        for d_p,y in enumerate(years):
            if type(y) == list:
                idx = []
                for yu in y:
                    idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if yu+(m>11)<bud_month.nyrs])
            else:
                idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if y+(m>11)<bud_month.nyrs]
            idx = [i for i in idx if not np.ma.core.is_masked(i)]
            if vector_anomaly:
                vec_plot0[1,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) - vec_plot0[1,0,:,:]
            else:
                vec_plot0[1,d_p+1,:,:]  = np.nanmean([bud_month.data[i] for i in idx], axis = 0) 
            d_p +=1
        
        for nr in range(nrows):
            vec_plot[0,nr],vec_plot[1,nr]=rvGplot.rg_vecs(vec_plot0[0,nr],vec_plot0[1,nr])
                
        
    for bn, budstr in enumerate(budget_state[:2]):
        bud_month.set_var(budstr)
        if budstr == 'thickness':
            vm = [0,hscale[0]]
            vmdiff = hscale[1]
            llevels = np.arange(0,hscale[0]+0.5,0.5)
            cmap = plt.cm.plasma
        elif budstr == 'concentration':
            vm = [0,1]
            vmdiff = 0.2
            llevels = np.array([0,0.5,0.75,0.9,0.95,1.0])
            cmap = plt.cm.YlGnBu 
        else:
            vm = [-1,1]
            llevels = np.arange(-1,1.2,0.2)
            cmap=plt.cm.RdYlBu
        ### first row is clim map
        ax = fig.add_subplot(nrows,3,bn+1,projection=GS.ccrs)
        ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        if type(depth_nn) is not bool:
            llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
            s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
        try:
            budget0 = bud_month.clim_map(periods=periods)
            budget0[budget0==0] = np.nan
            s = ax.contourf(GS.xpts, GS.ypts, budget0, vmin=vm[0], vmax=vm[1],levels=llevels,cmap=cmap, extend='max')
            cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
            cax.ax.yaxis.set_ticks_position('left')
            cax.set_label('Climatology')
        except TypeError:
            pass
        if drift_vectors:
            try:
                ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                          vec_plot[0,0],vec_plot[1,0],
                          scale = ra,width=rv[2],alpha=rv[3])
            except TypeError:
                pass
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        
        #### titles and colorbars here
        ax.set_title(label_state[bn])
            
        for nr,y in enumerate(years):
            pn = bn + 3*(nr+1) + 1 
            if type(y) == list:
                idx = []
                for yu in y:
                    idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if yu+(m>11)<bud_month.nyrs])
            else:
                idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                              if y+(m>11)<bud_month.nyrs]
            idx = [i for i in idx if not np.ma.core.is_masked(i)]
            ax = fig.add_subplot(nrows,3,pn,projection=GS.ccrs)
            ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
            if type(depth_nn) is not bool:
                llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
                s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
            
            try:
                budget = np.nanmean([bud_month.data[i] for i in idx], axis = 0) - budget0
                budget[budget==0] = np.nan
                llevels = np.arange(-vmdiff,vmdiff+0.02,0.02)
                s = ax.contourf(GS.xpts, GS.ypts, budget, vmin=-vmdiff, vmax=vmdiff,levels=llevels,cmap=plt.cm.RdYlBu, extend='both')
                cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
                cax.ax.yaxis.set_ticks_position('left')
                d0 = bud_month.dates[idx[0]]
                d1 = bud_month.dates[idx[-1]]
                if d0.year == d1.year:
                    clabel = d0.strftime('%Y')+' diff.'
                else:
                    clabel = d0.strftime('%Y/')+d1.strftime('%Y')+' diff.'
                cax.set_label(clabel)
            except TypeError:
                pass
            if drift_vectors:
                ra_use = copy.copy(ra)
                if vector_anomaly:
                    ra_use *= 0.5
                try:
                    ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                              vec_plot[0,nr+1],vec_plot[1,nr+1],
                              scale = ra_use,width=rv[2],alpha=rv[3])
                except TypeError:
                    pass
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.COASTLINE)
            
    ### first row is clim map
    ### repeat again for velocity
    bn =2
    vmdiff = vscale[1]
    ax = fig.add_subplot(nrows,3,bn+1,projection=GS.ccrs)
    ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
    if type(depth_nn) is not bool:
        llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
        s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
    try:
        budget0 = np.hypot(vec_plot0[0,0],vec_plot0[1,0])
        budget0[budget0==0] = np.nan
        llevels = np.arange(0,vscale[0]+0.02,0.02)
        s = ax.contourf(GS.xpts, GS.ypts, budget0, vmin=0, vmax=vscale[0], levels=llevels,cmap=plt.cm.viridis, extend='both')
        cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
        cax.ax.yaxis.set_ticks_position('left')
        cax.set_label('Climatology')
    except TypeError:
        pass
    if drift_vectors:
        try:
            ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                      vec_plot[0,0],vec_plot[1,0],
                      scale = ra,width=rv[2],alpha=rv[3])
        except TypeError:
            pass
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)

    #### titles and colorbars here
    ax.set_title('Drift speed m/s')

    for nr,y in enumerate(years):
        pn = bn + 3*(nr+1) + 1 
        if type(y) == list:
            idx = []
            for yu in y:
                idx.extend([bud_month.yrpd[yu+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                          if yu+(m>11)<bud_month.nyrs])
        else:
            idx = [bud_month.yrpd[y+(m>11),np.mod(m,12)] for m in range(mstart,mstart+nmonths)
                          if y+(m>11)<bud_month.nyrs]
        idx = [i for i in idx if not np.ma.core.is_masked(i)]
        ax = fig.add_subplot(nrows,3,pn,projection=GS.ccrs)
        ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        if type(depth_nn) is not bool:
            llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
            s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)

        try:
            xvel = np.nanmean([bud_month.ice_drift_x[i] for i in idx], axis = 0)
            yvel = np.nanmean([bud_month.ice_drift_y[i] for i in idx], axis = 0)
            budget = np.hypot(xvel,yvel) - budget0
            budget[budget==0] = np.nan
            llevels = np.arange(-vmdiff,vmdiff+0.02,0.01)
            s = ax.contourf(GS.xpts, GS.ypts, budget, vmin=-vmdiff, vmax=vmdiff,levels=llevels,cmap=plt.cm.RdYlBu, extend='both')
            cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
            cax.ax.yaxis.set_ticks_position('left')
            d0 = bud_month.dates[idx[0]]
            d1 = bud_month.dates[idx[-1]]
            if type(y) == list:
                clabel = year_list_labels[nr]
            elif d0.year == d1.year:
                clabel = d0.strftime('%Y')+' diff.'
            else:
                clabel = d0.strftime('%Y/')+d1.strftime('%Y')+' diff.'
            cax.set_label(clabel)
        except TypeError:
            pass
        if drift_vectors:
            ra_use = copy.copy(ra)
            if vector_anomaly:
                ra_use *= 0.5
            try:
                ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                          vec_plot[0,nr+1],vec_plot[1,nr+1],
                          scale = ra_use,width=rv[2],alpha=rv[3])
            except TypeError:
                pass
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
            
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02)
    if type(save_dir) == str: 
        if type(years[0]) == list:
        # list the years
            year_list = [str(bud_month.dates[bud_month.yrpd[y,mstart]].year) for y in years[0]]
        else:
            year_list = [str(bud_month.dates[bud_month.yrpd[y,mstart]].year) for y in years]
        fig.savefig(save_dir+plot_label+'-'.join(year_list)+'_climatology'+plot_type,bbox_inches='tight')
#     plt.subplots_adjust(wspace=0.05, hspace=0.05)
    if show_plot: plt.show()

def budget_clim_lines(bud,mstart,nmonths,G,years=[],line_colors=[],t_set=[],line_labels = [],save_dir = False,show_plot = True,plot_type='.pdf',plot_label = '',mlabel_step = 2,extend_div= False):
    lines = []
    if extend_div:
        ### add another two rows to the plots
        nrows = 3
        ydim = 12
        var_list = budget_load_default[:4]+['positive_divergence','negative_divergence']
    else:
        nrows = 2
        ydim = 8
        var_list = budget_load_default[:4]
    if bud.new_files:
        pscale = 1e-12
    else:
        pscale = 1e-12*60*60*24
    mlist = [dt.datetime(2000,1,1) + relativedelta(months = mn) 
         for mn in range(mstart,mstart+nmonths)]
    mstr = [md.strftime('%b') for n,md in enumerate(mlist)]
    for i in range(nmonths):
        if np.mod(i,mlabel_step) != 0:
            mstr[i] = ''
    f= plt.figure(figsize=(8,ydim))
    for n, var in enumerate(var_list):
        ax = f.add_subplot(nrows,2,n+1)
        lines = []
        bud.set_var(var)
        bud.data *= np.array(bud.days_per_slice)[:,None,None]/30 ## convert to monthly growth 
        plot_temp = bud.clim_mean(first_p=mstart,time_set=t_set,method='ysum',mask=True,
                                        mult_array=G.xdist*G.ydist*pscale,full_time = True)
        l, = ax.plot(plot_temp,'--k',linewidth = 3)#,
#                      label = 'Mean')

        lines.append(l)
        #### loop over all the years
        for y in range(bud.nyrs):
            if not bud.yrpd.mask[y,mstart]:
                t_p = bud.yrpd[y,mstart]
                plot_temp = np.nansum(bud[t_p:t_p+nmonths]*G.xdist[None,:,:]*G.ydist[None,:,:]*pscale,axis=(1,2))
#                 l, = 
                ax.plot(plot_temp,'k',alpha= 0.5,linewidth=0.5)
        for y_care,ls in zip(years,line_colors):
            if type(y_care)==list:
                ### accumlate averages for those in list
                plot_temp = []
                for y in y_care:
                    t_p = bud.yrpd[y,mstart]
                    plot_temp.append(np.nansum(bud[t_p:t_p+nmonths]*G.xdist[None,:,:]*G.ydist[None,:,:]*pscale,axis=(1,2)))
                plot_temp = np.nanmean(plot_temp,axis = 0)
            else:
                t_p = bud.yrpd[y_care,mstart]
                plot_temp = np.nansum(bud[t_p:t_p+nmonths]*G.xdist[None,:,:]*G.ydist[None,:,:]*pscale,axis=(1,2))
            l, = ax.plot(plot_temp,color=ls,)
#                         label = bud.print_date(bud.yrpd[y_care,mstart],string='%Y'))
            lines.append(l)
#         ax.legend(lines)
        if len(line_labels) == 0 and np.any([type(y)==list for y in years]):
            #### make auto labels if needed and groups of years
            line_labels = ['set "line_labels"']*np.shape(years)[0]
        elif len(line_labels) == 0:
            line_labels = [gl.print_date(gl.yrpd[y_care,mstart],string='%Y') 
                           for y_care in years]
        ax.legend(lines,['Mean']+line_labels)
#         ax.legend(lines)
        ax.set_title(plot_label+var)
        ax.axes.set_xlim([0,nmonths-1])
        if n==0 or n==2:
            ax.set_ylabel('Volume change 10$^6$ km$^3$')
        if n==2 or n==3:
            ax.xaxis.set_ticks(range(0,nmonths))
            ax.xaxis.set_ticklabels(mstr)
        else:
            ax.xaxis.set_ticks([])
    if type(save_dir) == str: 
        f.savefig(save_dir+plot_label+'clim_lines'+plot_type,bbox_inches='tight')
    if show_plot: plt.show()


def gate_lines(gate_list,mstart,nmonths,years=[],line_colors=[],t_set=[],line_labels = [],save_dir = False,show_plot = True,plot_type='.pdf',plot_label = '',mlabel_step = 2,nrows=1):
    lines = []
    mlist = [dt.datetime(2000,1,1) + relativedelta(months = mn) 
         for mn in range(mstart,mstart+nmonths)]
    mstr = [md.strftime('%b') for n,md in enumerate(mlist)]
    for i in range(nmonths):
        if np.mod(i,mlabel_step) != 0:
            mstr[i] = ''
    n_gates = len(gate_list)
    ncols = int(np.ceil(n_gates/nrows))
    lines = []
    f= plt.figure(figsize=(ncols*3,nrows*3))
    for n,gl in enumerate(gate_list):
        if gl.new_files:
            pscale = 1e-12
        else:
            pscale = 1e-12*60*60*24
        ax = f.add_subplot(nrows,ncols,n+1)
        gl.set_var('transport_x')
    # bud.data *= np.array(bud.days_per_slice)[:,None,None]
        plot_temp = gl.clim_mean(first_p=mstart,method='ysum',full_time=True,time_set=t_set,
                                   mult_array=gl.G.xdist*pscale/2)
        l, = ax.plot(np.array(plot_temp),'--k',linewidth= 4)
        lines.append(l)

        #### loop over all the years
        for y in range(gl.nyrs):
            if not gl.yrpd.mask[y,mstart]:
                t_p = gl.yrpd[y,mstart]
                plot_temp = np.nansum(gl[t_p:t_p+nmonths]*gl.G.xdist[None,:,:]*pscale/2,axis=(1,2))
                l, = ax.plot(plot_temp,'k',alpha= 0.5)

        for y_care,ls in zip(years,line_colors):
            if type(y_care)==list:
                ### accumlate averages for those in list
                plot_temp = []
                for y in y_care:
                    t_p = gl.yrpd[y,mstart]
                    plot_temp.append(np.nansum(gl[t_p:t_p+nmonths]*gl.G.xdist[None,:,:]*pscale/2,axis=(1,2)))
                plot_temp = np.nanmean(plot_temp,axis = 0)
                
            else:
                t_p = gl.yrpd[y_care,mstart]
                plot_temp = np.nansum(gl[t_p:t_p+nmonths]*gl.G.xdist[None,:,:]*pscale/2,axis=(1,2))
            l, = ax.plot(plot_temp,color=ls,)
            lines.append(l)
        if len(line_labels) == 0 and np.any([type(y)==list for y in years]):
            #### make auto labels if needed and groups of years
            line_labels = ['set "line_labels"']*np.shape(y_care)[0]
        elif len(line_labels) == 0:
            line_labels = [gl.print_date(gl.yrpd[y_care,mstart],string='%Y') 
                           for y_care in years]
            
        ax.legend(lines,['Mean']+line_labels)
#         ax.legend(lines)
        ax.set_title(gl.name)
        ax.axes.set_xlim([0,nmonths-1])
        if np.mod(n,ncols) ==0:
            ax.set_ylabel('Volume change 10$^6$ km$^3$')
        if n>=(nrows-1)*ncols:
            ax.xaxis.set_ticks(range(0,nmonths))
            ax.xaxis.set_ticklabels(mstr)
        else:
            ax.xaxis.set_ticks([])
    f.subplots_adjust( wspace=0.2)
    if type(save_dir) == str: 
        f.savefig(save_dir+plot_label+'Gate_lines'+plot_type,bbox_inches='tight')
    if show_plot: plt.show()
        

def budget_clim_lines_norm(bud,mstart,nmonths,G,years=[],line_colors=[],t_set=[],line_labels = [],save_dir = False,show_plot = True,plot_type='.pdf',plot_label = '',mlabel_step = 2,extend_div= False,norm_conc = False):
    lines = []
    if extend_div:
        ### add another two rows to the plots
        nrows = 3
        ydim = 12
        var_list = budget_load_default[:4]+['positive_divergence','negative_divergence']
    else:
        nrows = 2
        ydim = 8
        var_list = budget_load_default[:4]
    if bud.new_files:
        pscale = 1e-12
    else:
        pscale = 1e-12*60*60*24
    mlist = [dt.datetime(2000,1,1) + relativedelta(months = mn) 
         for mn in range(mstart,mstart+nmonths)]
    mstr = [md.strftime('%b') for n,md in enumerate(mlist)]
    for i in range(nmonths):
        if np.mod(i,mlabel_step) != 0:
            mstr[i] = ''
    f= plt.figure(figsize=(8,ydim))
    for n, var in enumerate(var_list):
        ax = f.add_subplot(nrows,2,n+1)
        lines = []
        bud.set_var(var)
        bud.data *= np.array(bud.days_per_slice)[:,None,None]/30 ## convert to monthly growth 
        if norm_conc:
            # scale all data by colocated IC
            bud.data = bud.data/bud.concentration
        plot_temp = bud.clim_mean(first_p=mstart,time_set=t_set,
                                  method='mean',mask=True,full_time = True)
        l, = ax.plot(plot_temp,'--k',linewidth = 3)
        lines.append(l)
        #### loop over all the years
        for y in range(bud.nyrs):
            if not bud.yrpd.mask[y,mstart]:
                t_p = bud.yrpd[y,mstart]
                plot_temp = np.nanmean(bud[t_p:t_p+nmonths],axis=(1,2))
#                 l, = 
                ax.plot(plot_temp,'k',alpha= 0.5,linewidth=0.5)
        for y_care,ls in zip(years,line_colors):
            if type(y_care)==list:
                ### accumlate averages for those in list
                plot_temp = []
                for y in y_care:
                    t_p = bud.yrpd[y,mstart]
                    plot_temp.append(np.nanmean(bud[t_p:t_p+nmonths],axis=(1,2)))
                plot_temp = np.nanmean(plot_temp,axis = 0)
            else:
                t_p = bud.yrpd[y_care,mstart]
                plot_temp = np.nanmean(bud[t_p:t_p+nmonths],axis=(1,2))
#                 l, = 
            l, = ax.plot(plot_temp,color=ls,)
#                         label = bud.print_date(bud.yrpd[y_care,mstart],string='%Y'))
            lines.append(l)
#         ax.legend(lines)
        if len(line_labels) == 0 and np.any([type(y)==list for y in years]):
            #### make auto labels if needed and groups of years
            line_labels = ['set "line_labels"']*np.shape(y_care)[0]
        elif len(line_labels) == 0:
            line_labels = [gl.print_date(gl.yrpd[y_care,mstart],string='%Y') 
                           for y_care in years]
        ax.legend(lines,['Mean']+line_labels)
#         ax.legend(lines)
        ax.set_title(plot_label+var)
        ax.axes.set_xlim([0,nmonths-1])
        if n==0 or n==2:
            ax.set_ylabel('Thickness Change (m/unit area)')
        if n==2 or n==3:
            ax.xaxis.set_ticks(range(0,nmonths))
            ax.xaxis.set_ticklabels(mstr)
        else:
            ax.xaxis.set_ticks([])
    if type(save_dir) == str: 
        f.savefig(save_dir+plot_label+'clim_lines'+plot_type,bbox_inches='tight')
    if show_plot: plt.show()


def plot_single(data,dlevels,GS,cmap=plt.cm.RdYlBu,label='BLANK',
                vecx = None,vecy=None,rvGplot=None,rv=(15,0.2,0.01,1),
                save_dir = False,show_plot = True,plot_type='.png',
                plot_label = '',depth_nn=False,return_ax=False,give_ax=None):
    bs0 = dlevels[0]
    bs1 = dlevels[-1]
    if vecx is not None:
        rm = int(GS.m/rv[0])
        rn = int(GS.n/rv[0])
        ra = np.sqrt(rm+rn)
        ra=ra*rv[1]
        vec_plot=rvGplot.rg_vecs(vecx,vecy)
    if give_ax is None:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(1,1,1,projection=GS.ccrs)
    else:
        ax = give_ax
    ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
    if type(depth_nn) is not bool:
        llevels = np.arange(-2000,0,250) # check etopo.ravel().max()
        s = ax.contour(GS.xpts, GS.ypts, -depth_nn, vmin=-8000, vmax=0,levels=llevels,colors='k', linestyles='-',alpha=0.05)
            
    #### data contour        
    s = ax.contourf(GS.xpts, GS.ypts, data, vmin=bs0, vmax=bs1,levels=dlevels,cmap=cmap, extend='both')
    cax = plt.colorbar(s,  pad=-0.12,shrink= 0.85)
    cax.ax.yaxis.set_ticks_position('left')
    cax.set_label(label)
    #### vectors
    if vecx is not None:
        ax.quiver(rvGplot.mesh_new[0],rvGplot.mesh_new[1],
                              vec_plot[0],vec_plot[1],
                              scale = ra,width=rv[2],alpha=rv[3])
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    if return_ax: return ax
    if type(save_dir) == str: 
        fig.savefig(save_dir+plot_label+'_single'+plot_type,
                  bbox_inches='tight')
    if show_plot: plt.show()


def load_gate_grids(gate_grid_dir,m,old_list = False):
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
    if not old_list:
        ### CAA_west
        CAA_west = gs.grid_set(m)
        ### CAA_east
        CAA_east = gs.grid_set(m)
        Gates_list.extend([
            [CAA_east,'Gate_CAA_east/'],
            [CAA_west,'Gate_CAA_west/'],])

    for gl in Gates_list:
        grid_name = gl[1].split('/')[0]+'.npz'
        mask_name = gl[1].split('/')[0]+'_mask.npz'
        gl[0].load_grid(gate_grid_dir+grid_name)
        gl[0].load_mask(gate_grid_dir+mask_name)
    return Gates_list


def load_gates(gate_data_dir,Gates_list,year_limit = None,var_load = ['transport_y'],seasonal=False,daily=False,ignore_file = [],verbos= False):
    Gate_data = []
    for gl in Gates_list:
        save_dir = gate_data_dir+gl[1]
        gate_bud = budget_result(save_dir,seasonal=seasonal,daily=daily,
                            var_load = var_load,year_limit = year_limit,ignore_file = ignore_file,verbos=verbos)
        gate_bud.G = gl[0]
        gate_bud.name = gl[1].split('Gate_')[1].split('/')[0]
        Gate_data.append(gate_bud)
    return Gate_data