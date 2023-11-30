
import numpy as np
import datetime as dt
import struct
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from os.path import exists
import glob
from scipy.interpolate import griddata
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree


class Pathfinder():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,grid=False,hemi='north'):
        self.name = 'Pathfinder'
        self.path = ppath
        self.hemi = hemi
        self.vyear_load = 0
        self.erryear_load = 0
        self.vels_loaded = False
        self.err_loaded = False
        if type(grid) == bool:
            self.check_grid = False
        else:
            self.check_grid = True
            self.grid = grid
        
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        dates =[]
        d0 = dt.datetime(1970,1,1)
        n_yrs = (time_end.year - time_start.year)+1
        for y in range(n_yrs):
            yu = time_start.year + y
            if self.hemi == 'north':
                f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
            elif self.hemi == 'south':
                ### icemotion_daily_sh_25km_20100101_20101231_v4.1.nc
                f_name = 'icemotion_daily_sh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
            if exists(self.path+f_name):
                f_nc = Dataset(self.path+f_name)
                [dates.append(d0 + relativedelta(days = d))
                     for d in f_nc['time'][:]]
                f_nc.close()
        self.dates = dates
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')

    # daily points in yearly files

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self,dates_u,verbos=False):
        d0 = dt.datetime(1970,1,1)
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.vyear_load != yu) or (not self.vels_loaded)):
                print('loading new year of data: '+str(yu))
                if self.hemi == 'north':
                    f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                elif self.hemi == 'south':
                    ### icemotion_daily_sh_25km_20100101_20101231_v4.1.nc
                    f_name = 'icemotion_daily_sh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
#                 f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                f_nc = Dataset(self.path+f_name)
        #         print(p0,p1)
                self.u = f_nc['u'][:]
                self.v = f_nc['v'][:]
#                 self.u[self.u> 1000] = np.nan
#                 self.u[self.u<-1000] = np.nan
#                 self.u[self.v> 1000] = np.nan
#                 self.u[self.v<-1000] = np.nan
                self.u[self.u.mask] = np.nan
                self.v[self.v.mask] = np.nan
                f_nc.close()
                self.vyear_load = yu
                self.vels_loaded= True
            p0 = dates_u[ 0].timetuple().tm_yday -1
            p1 = dates_u[-1].timetuple().tm_yday 
            datau = self.u[p0:p1,:,:].transpose((0,2,1))/100
            datav = self.v[p0:p1,:,:].transpose((0,2,1))/100
            if self.check_grid:
                for n in range(np.shape(datau)[0]):
                    datau[n][self.grid.lats>88] = np.nan
                    datav[n][self.grid.lats>88] = np.nan
            return datau,datav
        
    def get_err(self,dates_u,verbos=False):
        ## errs need to be dimensional
        ### always make sure you call get_vels frist
        d0 = dt.datetime(1970,1,1)
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.erryear_load != yu) or (not self.err_loaded)):
                print('loading new year of data: '+str(yu))
#                 f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                if self.hemi == 'north':
                    f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                elif self.hemi == 'south':
                    ### icemotion_daily_sh_25km_20100101_20101231_v4.1.nc
                    f_name = 'icemotion_daily_sh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                f_nc = Dataset(self.path+f_name)
        #         print(p0,p1)
                self.err = f_nc['icemotion_error_estimate'][:]
#                 self.err = self.err/100.0 #percent to fraction
                self.err= np.sqrt(self.err)/100#### sqrt then/100 to m/s
                self.err[self.err.mask] = np.nan
                self.err[self.err<0] = np.nan
#                 self.uerr = self.err*self.u
#                 self.verr = self.err*self.v
                f_nc.close()
                self.erryear_load = yu
                self.err_loaded= True
            p0 = dates_u[ 0].timetuple().tm_yday -1
            p1 = dates_u[-1].timetuple().tm_yday 
            data_err = self.err[p0:p1,:,:].transpose((0,2,1))
#             data_uerr = self.uerr[p0:p1,:,:].transpose((0,2,1))
#             data_verr = self.verr[p0:p1,:,:].transpose((0,2,1))
            if self.check_grid:
                for n in range(np.shape(datau)[0]):
                    data_err[n][self.grid.lats>88] = np.nanmean
#                     data_uerr[n][self.grid.lats>88] = np.nanmean
#                     data_verr[n][self.grid.lats>88] = np.nanmean
            return data_err
#             return data_uerr,data_verr

class Pathfinder_weekly():
    """
    forcing class for the budget
    lets the forcing load efficiently
    assumes a single NRT nc file for all the Pathfinder
    """
    def __init__(self,ppath,grid=False):
        self.name = 'Pathfinder_weekly'
        self.path = ppath
        self.files = glob.glob(ppath+'*.nc')
        self.files.sort()
        self.vels_loaded = False
        self.vyear_load = 0
        
        if type(grid) == bool:
            self.check_grid = False
        else:
            self.check_grid = True
            self.grid = grid
        
        
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        self.ds = []
        self.de = []
        self.all_dates = []
        self.dates = []
        d0 = dt.datetime(1970,1,1)
        for f in self.files:
            Pnc = Dataset(f)
            ds = d0 + relativedelta(days =  int(Pnc.variables['time'][0]))
            de = d0 + relativedelta(days =  int(Pnc.variables['time'][-1]))
#             print(ds.year,time_start.year,de.year,time_end.year)
            if ds.year>=time_start.year and de.year<=time_end.year:
                self.ds.append(ds)
                self.de.append(de)
                self.all_dates.append([d0 + relativedelta(days = int(d)) for d in Pnc.variables['time']])
                self.dates.extend([t for t in self.all_dates[-1] if t>=time_start and t<=time_end])
            Pnc.close()
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self,dates_u,verbos=False):
        set_var_now = False
        join_data = False
        
        u_arrays = []
        v_arrays = []
        ### loop the dates
        for d in dates_u:
            iload = np.searchsorted(np.array(self.ds),d,side='right')-1
            if (d.year != self.vyear_load) or (not self.vels_loaded):
                self.load_vels([d],verbos=verbos)
            ### fill the out arrays
            load_point = self.all_dates[iload].index(d) 
            u_arrays.append(self.u[load_point])
            v_arrays.append(self.v[load_point])
        

        datau = np.ma.array(u_arrays)
        datav = np.ma.array(v_arrays)
        datau[datau.mask] = np.nan
        datav[datav.mask] = np.nan
        if self.check_grid:
            for n in range(np.shape(datau)[0]):
                datau[n][self.grid.lats>88] = np.nan
                datav[n][self.grid.lats>88] = np.nan
        return datau,datav

    def load_vels(self,dates_u,verbos=False):
        iload = np.searchsorted(np.array(self.ds),dates_u[0],side='right')-1
        Pnc = Dataset(self.files[iload])
        if verbos: print(self.files[iload])
        self.u = Pnc.variables['u'][:].transpose((0,2,1))/100
        self.v = Pnc.variables['v'][:].transpose((0,2,1))/100
        self.vels_loaded = True
        self.vyear_load = self.ds[iload].year
        Pnc.close()
        print('Filling buffer from '+self.ds[iload].strftime('%Y-%m-%d'))
        

class PIOMAS():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,time_smooth = 0):
        self.name = 'PIOMAS'
        self.check_years = []
        self.path = ppath
        self.vyear_load = 0
        self.hyear_load = 0
        self.ayear_load = 0
        self.vels_loaded = False
        self.hi_loaded = False
        self.aice_loaded = False
        self.tsmth = time_smooth
        
    def check_dates_in_file(self,data_f,ndata=1):
        ### return the no. of days for files listed
        with open(data_f, mode='rb') as file:
            filecontent = file.read()
            data = struct.unpack("f" * (len(filecontent)// 4), filecontent)
        n_points=np.asarray(data).shape[0]
        return int(n_points/120/360/ndata)
        
        print('loading new year of data: '+str(yu))
        
    def get_dates(self,time_start,time_end,check_years = []):
        """
        returns the all encompassing date list for use with the forcing object
        PIOMAS is a standardised list so we can just build a list
        Now adding a check years option - this allows us to list specific years to check file size
        For these particular years the BI will open an aice file to see how many days are in it and adjust accordingly
        """
        self.check_years = check_years
        dates =[]
        n_yrs = (time_end.year - time_start.year)-1
        if n_yrs>-1:
            y0 = dt.datetime(time_start.year,1,1)
            ye = dt.datetime(time_start.year,12,31)
            data_f = self.path+'aiday.H'+y0.strftime('%Y')
            if exists(data_f):
                for d in range(time_start.timetuple().tm_yday-1,
                               ye.timetuple().tm_yday):
                    dates.append(y0 +  relativedelta(days = d))
            for y in range(n_yrs):
                y0 += relativedelta(years=1)
                ye += relativedelta(years=1)
                data_f = self.path+'aiday.H'+y0.strftime('%Y')
                if exists(data_f) and y0.year in check_years:
                    n_days = self.check_dates_in_file(data_f)
                    for d in range(n_days):
                        dates.append(y0 +  relativedelta(days = d))
                elif exists(data_f):
                    for d in range(ye.timetuple().tm_yday):
                        dates.append(y0 +  relativedelta(days = d))
            y0 += relativedelta(years=1)
            ye  = time_end
            data_f = self.path+'aiday.H'+y0.strftime('%Y')
            if exists(data_f) and y0.year in check_years:
                n_days = self.check_dates_in_file(data_f)
                for d in range(n_days):
                    dates.append(y0 +  relativedelta(days = d))
            elif exists(data_f):
                for d in range(ye.timetuple().tm_yday):
                    dates.append(y0 +  relativedelta(days = d))
        else:
            y0 = dt.datetime(time_start.year,1,1)
            data_f = self.path+'aiday.H'+y0.strftime('%Y')
            if exists(data_f) and y0.year in check_years:
                n_days = self.check_dates_in_file(data_f)
                for d in range(n_days):
                    dates.append(y0 +  relativedelta(days = d))
            if exists(data_f):
                for d in range(time_start.timetuple().tm_yday-1,
                               time_end.timetuple().tm_yday):
                    dates.append(y0 +  relativedelta(days = d))

        self.dates= dates
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')

    # daily points in yearly files

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.vyear_load != yu) or (not self.vels_loaded)):
                print('loading new year of data: '+str(yu))
                data_f = self.path+'uiday.H'+str(yu)
                
                with open(data_f, mode='rb') as file:

                    filecontent = file.read()
                    data = struct.unpack("f" * (len(filecontent)// 4), filecontent)
                data_load = np.asarray(data)
                ndays = int(data_load.shape[0]/2/120/360)
                self.vels=data_load.reshape(ndays,2,120,360)
                self.vyear_load = yu
                self.vels_loaded= True

            p0 = dates_u[ 0].timetuple().tm_yday -1
            p1 = dates_u[-1].timetuple().tm_yday 
    #         print(p0,p1)
            if self.tsmth < 1:
                datau = self.vels[p0:p1,0,:,:].transpose((0,2,1))
                datav = self.vels[p0:p1,1,:,:].transpose((0,2,1))
                return datau,datav
            elif np.shape(dates_u)[0]>1:
                print('time smoothing not compatible with multiple dates')
                datau = self.vels[p0:p1,0,:,:].transpose((0,2,1))
                datav = self.vels[p0:p1,1,:,:].transpose((0,2,1))
                return datau,datav
            else:
                #### each time slice is the mean of 2*tsmth+1
                p0 = np.maximum(p0-self.tsmth,0)
                datau = self.vels[p0:p1+self.tsmth,0,:,:].transpose((0,2,1))
                datav = self.vels[p0:p1+self.tsmth,1,:,:].transpose((0,2,1))
#                 print(np.shape(datau))
                datau2 = np.expand_dims(np.nanmean(datau,axis=0),0)
                datav2 = np.expand_dims(np.nanmean(datav,axis=0),0)
                return datau2,datav2

    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.hyear_load != yu) or (not self.hi_loaded)):
                print('loading new year of data: '+str(yu))
                data_f = self.path+'hiday.H'+str(yu)
                with open(data_f, mode='rb') as file:

                    fileContent = file.read()
                    data = struct.unpack("f" * (len(fileContent)// 4), fileContent)
                data_load = np.asarray(data)
                ndays = int(data_load.shape[0]/120/360)
                self.hi=data_load.reshape(ndays,120,360)
                self.hyear_load = yu
                self.hi_loaded  = True

            p0 = dates_u[ 0].timetuple().tm_yday -1
            p1 = dates_u[-1].timetuple().tm_yday 
    #         print(p0,p1)
            if self.tsmth < 1:
                data = self.hi[p0:p1,:,:].transpose((0,2,1))
                return data
            elif np.shape(dates_u)[0]>1:
                print('Time smoothing not compatible with multiple dates')
                data = self.hi[p0:p1,:,:].transpose((0,2,1))
                return data
            else:
                #### each time slice is the mean of 2*tsmth+1
                p0 = np.maximum(p0-self.tsmth,0)
                data = self.hi[p0:p1+self.tsmth,0,:,:].transpose((0,2,1))
                data2 = np.expand_dims(np.nanmean(data,axis=0),0)
                return data2

    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.ayear_load != yu) or (not self.aice_loaded)):
                print('loading new year of data: '+str(yu))
                data_f = self.path+'aiday.H'+str(yu)
                with open(data_f, mode='rb') as file:

                    fileContent = file.read()
                    data = struct.unpack("f" * (len(fileContent)// 4), fileContent)
                data_load = np.asarray(data)
                ndays = int(data_load.shape[0]/120/360)
                self.aice=data_load.reshape(ndays,120,360)
                self.ayear_load = yu
                self.aice_loaded= True

            p0 = dates_u[ 0].timetuple().tm_yday -1
            p1 = dates_u[-1].timetuple().tm_yday 
    #         print(p0,p1)
            if self.tsmth < 1:
                data = self.aice[p0:p1,:,:].transpose((0,2,1))
                return data
            elif np.shape(dates_u)[0]>1:
                print('Time smoothing not compatible with multiple dates')
                data = self.aice[p0:p1,:,:].transpose((0,2,1))
                return data
            else:
                #### each time slice is the mean of 2*tsmth+1
                p0 = np.maximum(p0-self.tsmth,0)
                data = self.aice[p0:p1+self.tsmth,0,:,:].transpose((0,2,1))
                data2 = np.expand_dims(np.nanmean(data,axis=0),0)
                return data2

class NSIDC_nt():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'NSIDC_n'
        self.path = ppath
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #daily files
        dimY = 304
        dimX = 448
        d_no = np.shape(dates_u)[0]
        data =  np.empty([d_no, dimX, dimY])
        for n,d in enumerate(dates_u):
#             infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year<2020:
# #                 infile = self.path+"/nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#                 infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year>2019:
# #                 infile = self.path+"/nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
#                 infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
#             if d>=dt.datetime(2020,11,1):
#                 infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
#             else:
#             if d.year>2019:
#                 infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year<2019:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year>2018:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_*.bin"
            flist  = glob.glob(infile)
            if len(flist) > 0:
                infile = flist[0]
            with open(infile, 'rb') as fr:
                hdr = fr.read(300)
                ice = np.fromfile(fr, dtype=np.uint8)

            ice = ice.reshape(dimX,dimY)
            ice = np.flipud(ice)
            data[n] = ice / 250.
        data[data>1.0] = np.nan
        return data

    def get_dates(self,time_start,time_end):
        # does dates_u cover one year or more
        #daily files
        dates_u = []
        d_no = (time_end-time_start).days +3 
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
#             if d>=dt.datetime(2020,11,1):
#                 infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
#             else:
#             if d.year>2019:
            infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_*.bin"
            flist  = glob.glob(infile)
            if len(flist) > 0:
                infile = flist[0]
#                 print(infile)
                # check infile exists 
#                 if exists(infile):
                dates_u.append(d)
            #if it does append dates_u
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')

        

class NSIDC_bt():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'NSIDC_b'
        self.path = ppath
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #daily files
        dimY = 304
        dimX = 448
        d_no = np.shape(dates_u)[0]
        data =  np.empty([d_no, dimX, dimY])
        for n,d in enumerate(dates_u):
            infile = self.path+d.strftime('/%Y/')+"bt_"+d.strftime('%Y%m%d')+"_f17_v3.1_n.bin"
#             if d.year<2019:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year>2018:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            with open(infile, 'rb') as fr:
                hdr = fr.read(0)
                ice = np.fromfile(fr, dtype="<i2")

            ice = ice.reshape(dimX,dimY)
            ice = np.flipud(ice)
            data[n] = ice / 1000.
        data[data>1.0] = np.nan
        return data

    def get_dates(self,time_start,time_end):
        # does dates_u cover one year or more
        #daily files
        dates_u = []
        d_no = (time_end-time_start).days +3 
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
#             if d.year<2019:
            infile = self.path+d.strftime('/%Y/')+"bt_"+d.strftime('%Y%m%d')+"_f17_v3.1_n.bin"
#             if d.year>2018:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            # check infile exists 
            if exists(infile):
                dates_u.append(d)
            #if it does append dates_u
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')



class AWI_SMOS_daily():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'AWI_SMOS_daily'
        self.path = ppath
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        blurb = 'W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_'
        dn = np.shape(dates_u)[0]
        d = dates_u[0]
        if d > dt.datetime(2020,6,1):
            blurb2 = '_r_v203_01_l4sit.nc'
        elif d > dt.datetime(2019,6,1):
            blurb2 = '_r_v202_02_l4sit.nc'
        else:
            blurb2 = '_r_v202_01_l4sit.nc'
        t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
        t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
        file = self.path+d.strftime('%Y/%m/')+blurb+t0+t1+blurb2
        f_nc = Dataset(file)
        hi = f_nc['analysis_sea_ice_thickness'][0]
        hi[hi.mask] = np.nan
        dx,dy = hi.shape 
        data = np.empty([dn,dx,dy])
        data[0] = hi
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            if d > dt.datetime(2020,6,1):
                blurb2 = '_r_v203_01_l4sit.nc'
            elif d > dt.datetime(2019,6,1):
                blurb2 = '_r_v202_02_l4sit.nc'
            else:
                blurb2 = '_r_v202_01_l4sit.nc'
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/%m/')+blurb+t0+t1+blurb2
            f_nc = Dataset(file)
            hi = f_nc['analysis_sea_ice_thickness'][0]
            hi[hi.mask] = np.nan
            data[n+1] = hi
            f_nc.close()
            #if it does append dates_u
        return data


    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        blurb = 'W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_'
        dn = np.shape(dates_u)[0]
        d = dates_u[0]
        if d > dt.datetime(2020,6,1):
            blurb2 = '_r_v203_01_l4sit.nc'
        elif d > dt.datetime(2019,6,1):
            blurb2 = '_r_v202_02_l4sit.nc'
        else:
            blurb2 = '_r_v202_01_l4sit.nc'
        t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
        t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
        file = self.path+d.strftime('/%Y/%m/')+blurb+t0+t1+blurb2
        f_nc = Dataset(file)
        aice = f_nc['sea_ice_concentration'][0]
        aice[aice.mask] = np.nan
        dx,dy = aice.shape 
        data = np.empty([dn,dx,dy])
        data[0] = aice/100
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            if d > dt.datetime(2020,6,1):
                blurb2 = '_r_v203_01_l4sit.nc'
            elif d > dt.datetime(2019,6,1):
                blurb2 = '_r_v202_02_l4sit.nc'
            else:
                blurb2 = '_r_v202_01_l4sit.nc'
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/%m/')+blurb+t0+t1+blurb2
            f_nc = Dataset(file)
            aice = f_nc['sea_ice_concentration'][0]
            aice[aice.mask] = np.nan
            data[n+1] = aice/100
            #if it does append dates_u
            f_nc.close()
        return data


    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_dates(self,time_start,time_end):
        # does dates_u cover one year or more
        #daily files
        blurb = 'W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_'
        dates_u = []
        dd = (time_end-time_start).days+1
        # make sure we get the bracket points
        for dn in range(dd):
            d = time_start+ relativedelta(days = dn)
            if d > dt.datetime(2020,6,1):
                blurb2 = '_r_v203_01_l4sit.nc'
            elif d > dt.datetime(2019,6,1):
                blurb2 = '_r_v202_02_l4sit.nc'
            else:
                blurb2 = '_r_v202_01_l4sit.nc'
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/%m/')+blurb+t0+t1+blurb2
            ## got a corrupted file on t0 = 20201109
            if t0 == '20201109_':
                print('AWI SMOS Avoiding 2020-11-09')
            elif exists(file):
                dates_u.append(d)
            #if it does append dates_u
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')

class AWI_weekly():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'AWI_weekly'
        self.path = ppath
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        blurb = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-'
        blurb2 = '-fv2p2.nc'
        dn = np.shape(dates_u)[0]
        d = dates_u[0]
        t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
        t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
        file = self.path+d.strftime('%Y/')+blurb+t0+t1+blurb2
        f_nc = Dataset(file)
        hi = f_nc['analysis_sea_ice_thickness'][0]
        hi[hi.mask] = np.nan
        dx,dy = hi.shape 
        data = np.empty([dn,dx,dy])
        data[0] = hi
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/')+blurb+t0+t1+blurb2
            f_nc = Dataset(file)
            hi = f_nc['analysis_sea_ice_thickness'][0]
            hi[hi.mask] = np.nan
            data[n+1] = hi
            f_nc.close()
            #if it does append dates_u
        return data


    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        blurb = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-'
        blurb2 = '-fv2p2.nc'
        dn = np.shape(dates_u)[0]
        d = dates_u[0]
        t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
        t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
        file = self.path+d.strftime('/%Y/')+blurb+t0+t1+blurb2
        f_nc = Dataset(file)
        aice = f_nc['sea_ice_concentration'][0]
        aice[aice.mask] = np.nan
        dx,dy = aice.shape 
        data = np.empty([dn,dx,dy])
        data[0] = aice/100
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/')+blurb+t0+t1+blurb2
            f_nc = Dataset(file)
            aice = f_nc['sea_ice_concentration'][0]
            aice[aice.mask] = np.nan
            data[n+1] = aice/100
            #if it does append dates_u
            f_nc.close()
        return data


    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_dates(self,time_start,time_end):
        # does dates_u cover one year or more
        # awi weekly - first ever week is 20101101 to 20101107
        # midpoint of 20101104
        w0 = dt.datetime(2010,11,4)
        # find first time point before time start
        w1 = int((time_start-w0).days/7)
        # find first time point after time end
        w2 = int((time_end-w0).days/7)+1
        print(w1,w2)
        #daily files
        blurb = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-'
        blurb2 = '-fv2p2.nc'
        dates_u = []
        # make sure we get the bracket points
        for w in range(w1,w2):
            d = w0+ relativedelta(days = w*7)
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/')+blurb+t0+t1+blurb2
            if exists(file):
                dates_u.append(d)
            #if it does append dates_u
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')



class AWI_monthly():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'AWI_monthly'
        self.path = ppath
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        blurb = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-'
        dn = np.shape(dates_u)[0]
        d = dates_u[0]
        dload = d.replace(day=1)
        file = self.path+dload.strftime('/%Y/')+blurb+dload.strftime('%Y%m')+'-fv2p2.nc'
        f_nc = Dataset(file)
        hi = f_nc['sea_ice_thickness'][0]
        dx,dy = hi.shape 
        data = np.empty([dn,dx,dy])
        data[0] = hi
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            dload = d.replace(day=1)
            file = self.path+dload.strftime('/%Y/')+blurb+dload.strftime('%Y%m')+'-fv2p2.nc'
            f_nc = Dataset(file)
            hi = f_nc['sea_ice_thickness'][0]
            data[n+1] = hi
            f_nc.close()
            #if it does append dates_u
        return data


    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        blurb = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-'
        dn = np.shape(dates_u)[0]
        d = dates_u[0]
        dload = d.replace(day=1)
        file = self.path+dload.strftime('/%Y/')+blurb+dload.strftime('%Y%m')+'-fv2p2.nc'
        f_nc = Dataset(file)
        aice = f_nc['sea_ice_concentration'][0]
        dx,dy = aice.shape 
        data = np.empty([dn,dx,dy])
        data[0] = aice/100
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            dload = d.replace(day=1)
            file = self.path+dload.strftime('/%Y/')+blurb+dload.strftime('%Y%m')+'-fv2p2.nc'
            f_nc = Dataset(file)
            aice = f_nc['sea_ice_concentration'][0]
            data[n+1] = aice/100
            #if it does append dates_u
            f_nc.close()
        return data


    # next function will take a list of dates and return an appropriately orientated arrays
    def get_dates(self,time_start,time_end,fill_end_months=False):
        # does dates_u cover one year or more
        #daily files
        blurb = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2-'
        dates_u = []
        dy = time_end.year-time_start.year
        dm = time_end.month-time_start.month
        m_no = dy*12 + dm +2
        # make sure we get the bracket points
        ts_m = dt.datetime(time_start.year,time_start.month,1)
        for mn in range(m_no):
            d = ts_m+ relativedelta(months = mn )
            file = self.path+d.strftime('/%Y/')+blurb+d.strftime('%Y%m')+'-fv2p2.nc'
            if exists(file):
                if d.month==2:
                    mid_day = 13
                else:
                    mid_day = 14
                dates_u.append(d + relativedelta(days=mid_day))
        ### now work over the date and adjust for summer
        ### remove all months = [5,6,8,7,9]
        ### also need hole end points to be at month end not mid
        self.dates= []
        month_keep = [1,2,3,4,10,11,12]
        for d in dates_u:
            if d.month in month_keep:
                if fill_end_months and d.month == 4:
                    self.dates.append(d)
                    d_end = d.replace(day=30)
                    self.dates.append(d_end)
                elif fill_end_months and d.month == 10:
                    d_start = d.replace(day=1)
                    self.dates.append(d_start)
                    self.dates.append(d)
                else:
                    self.dates.append(d)
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')

class OSISAF():
    """
    forcing class for the budget
    lets the forcing load efficiently
    version is 405 for older, 455 for the new data (2022)
    
    """
    def __init__(self,ppath,version='405',hemi='north'):
        self.path = ppath
        self.name = 'osisaf'
        self.version = version
        self.hemi = hemi
        if   self.version =='405': 
            self.read = self.read_405
            self.fname = self.fname_405
        elif self.version =='455': 
            self.read = self.read_455
            self.fname = self.fname_455
            
    def fname_405(self,d):
        t1 = (d - relativedelta(hours=12)).strftime('%Y%m%d%H00')
        t2 = (d + relativedelta(hours=36)).strftime('%Y%m%d%H00')
        # check if were on the last day of the month
        next_month = d.replace(day=28) + dt.timedelta(days=4)  # this will never fail
        next_month.replace(day = 1)
        if (next_month - d).days>1:
            td = d.strftime('%Y/%m/')
        else:
            td = (d + relativedelta(months=1)).strftime('%Y/%m/')
        if self.hemi == 'north':
            f_name = td+'ice_drift_nh_polstere-625_multi-oi_'+t1+'-'+t2+'.nc'
        elif self.hemi == 'south':
            f_name = td+'ice_drift_sh_polstere-625_multi-oi_'+t1+'-'+t2+'.nc'
        return f_name
    
    def fname_455(self,d):
        ### ice_drift_nh_ease2-750_cdr-v1p0_24h-200203041200.nc'
#         t1 = (d - relativedelta(hours=12)).strftime('%Y%m%d%00')
        t2 = (d + relativedelta(hours=12)).strftime('%Y%m%d%H00')
        # check if were on the last day of the month
        next_month = d.replace(day=28) + dt.timedelta(days=4)  # this will never fail
        next_month.replace(day = 1)
        if (next_month - d).days>0:
            td = d.strftime('%Y/%m/')
        else:
            td = (d + relativedelta(months=1)).strftime('%Y/%m/')
        if self.hemi == 'north':
            f_name = td+'ice_drift_nh_ease2-750_cdr-v1p0_24h-'+t2+'.nc'
        elif self.hemi == 'south':
            f_name = td+'ice_drift_sh_ease2-750_cdr-v1p0_24h-'+t2+'.nc'
        return f_name 
    
    def read_405(self,u,v,d):
        vmult = -1.0
        if d > dt.datetime(2015,9,12): vmult = 1.0
        unow = np.fliplr(u.T)
        vnow = np.fliplr(v.T)
        unow[unow.mask] = np.nan
        vnow[vnow.mask] = np.nan
        # convert km/48hrs to m/s
        u =  unow*1e3/60/60/48
        v = vmult*vnow*1e3/60/60/48
        
        return u,v
    
    def read_455(self,u,v,d):
        # convert km/24hrs to m/s
        u[u.mask] = np.nan
        v[v.mask] = np.nan
        unow = -v*1e3/60/60/24
        vnow =  u*1e3/60/60/24
        
        return unow,vnow
            
# give a 
# next function will take a list of dates and return an appropriately orientated arrays
    def get_dates(self,time_start,time_end):
        # individual drift for each day slice
        dates_u = []
        d_no = (time_end-time_start).days +3 
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
            f_name = self.fname(d)
            if exists(self.path+f_name):
                dates_u.append(d)
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')



    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self,dates_u,verbos=False):
        # individual drift for each day slice
        d_no = np.shape(dates_u)[0]
        d = dates_u[0]
        f_name = self.fname(d)
        f_nc = Dataset(self.path+f_name)
        unow = f_nc['dX'][0]
        vnow = f_nc['dY'][0]
        u,v = self.read(unow,vnow,d)
        dx,dy = np.shape(unow)
        data_u = np.zeros([d_no,dx,dy])
        data_v = np.zeros([d_no,dx,dy])
        # convert km/48hrs to m/s
        data_u[0] = u
        data_v[0] = v
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            f_name = self.fname(d)
            f_nc = Dataset(self.path+f_name)
            unow = f_nc['dX'][0]
            vnow = f_nc['dY'][0]
            u,v = self.read(unow,vnow,d)
            data_u[n+1] = u
            data_v[n+1] = v
            f_nc.close()
        return data_u,data_v
    
    def get_err(self,dates_u,verbos=False):
        ### only for 455 at the moment
        ## errs need to be dimensional
        ### always make sure you call get_vels frist
        d_no = np.shape(dates_u)[0]
        d = dates_u[0]
        f_name = self.fname(d)
        f_nc = Dataset(self.path+f_name)
        u_err = f_nc['uncert_dX_and_dY'][0]
        u_err[u_err.mask] = np.nan
        dx,dy = np.shape(u_err)
        data_err = np.zeros([d_no,dx,dy])
        # convert km/24hrs to m/s
        data_err[0] = u_err*1e3/60/60/24
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            f_name = self.fname(d)
            f_nc = Dataset(self.path+f_name)
            u_err = f_nc['uncert_dX_and_dY'][0]
            u_err[u_err.mask] = np.nan
            data_err[n+1] = u_err*1e3/60/60/24
            f_nc.close()
        return data_err


class CPOM_hi:
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,G):
        self.name = 'CPOM_monthly'
        self.path = ppath
        ### we need to regrid as we load
        ### for efficiency we will buffer two months
        self.hi1 = 0.0
        self.dl1 = None
        self.hi2 = 0.0
        self.dl2 = None
        self.G = G
        self.THRESH = np.hypot(np.mean(np.diff(self.G.ypts)),
                               np.mean(np.diff(self.G.xpts)))/2
#         self.edges_x = self.G.xpts[:,int(G.n/2)] - self.G.dxRes/2
#         self.edges_y = self.G.ypts[int(G.m/2),:] - self.G.dyRes/2
#         self.edges_x = np.append(self.edges_x,2*self.G.xpts[-1,int(G.n/2)] - self.G.xpts[-2,int(G.n/2)])
#         self.edges_y = np.append(self.edges_y,2*self.G.ypts[int(G.m/2),-1] - self.G.ypts[int(G.m/2),-2])
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #dmonthly
        reload = [False,False]
        d_no = np.shape(dates_u)[0]
        if d_no>2: 
            print('Warning CPOM hi not compatible with tw over 1 month')
            return None
        ### first check if any of the dates in dates_u
        ### match the preloaded, ddL
        ### if the first date matches the second
        data_u = np.zeros([d_no,self.G.m,self.G.n])
        for n,d in enumerate(dates_u):
            if d == self.dl1:
                ### return buffered
                data_u[n] = self.hi1
            elif d == self.dl2:
                ### return buffered
                data_u[n] = self.hi2
            else:
                ### load and bin
                ### check dates
                ### normal month, data on the first
                dload = d.replace(day=1)
                file = self.path+d.strftime('%Y%m_')+'Thick.map'
                if verbos: print(file)
                f = np.genfromtxt(file)
                hi = f[:,2]
                lon = f[:,1]
                lat = f[:,0]
#                 xpts,ypts  = self.G.mplot(lon,lat)
#                 print('Binning awkward CPOM data')
# #              
#                 ret = binned_statistic_2d( xpts, ypts, hi,
#                     statistic='mean', bins=[self.edges_x,self.edges_y]) 
#                 data_u[n] = ret.statistic
#                 ret = binned_statistic_2d( xpts, ypts, [],
#                     statistic='count', bins=[self.edges_x,self.edges_y]) 
#                 data_u[n][ret.statistic<4] = np.nan
                print('Regridding awkward CPOM data')
                xpts,ypts = self.G.mplot(lon,lat)
                xy = np.vstack([xpts,ypts]).T
                data_u[n] = griddata(xy, hi, (self.G.xpts, self.G.ypts),method='nearest')

                # Construct kd-tree, functionality copied from scipy.interpolate
                ### this finds empty values on the grid
                tree = cKDTree(xy)
                xi = _ndim_coords_from_arrays((self.G.xpts, self.G.ypts))
                dists, indexes = tree.query(xi)
                # Copy original result but mask missing values with NaNs
                data_u[n][dists > self.THRESH] = np.nan
                reload[n] = True
        ### update buffered thickness
        if d_no == 2:
            self.dl1 = dates_u[0]
            self.hi1 = data_u[0]
            self.dl2 = dates_u[1]
            self.hi2 = data_u[1]
        else:
            self.dl2 = dates_u[0]
            self.hi2 = data_u[0]
        return data_u

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_dates(self,time_start,time_end,fill_end_months=False):
        # does dates_u cover one year or more
        #daily files
        dates_u = []
        dy = time_end.year-time_start.year
        dm = time_end.month-time_start.month
        m_no = dy*12 + dm +2
        # make sure we get the bracket points
        ts_m = dt.datetime(time_start.year,time_start.month,1)
        for mn in range(m_no):
            d = ts_m+ relativedelta(months = mn )
            file = self.path+d.strftime('%Y%m_')+'Thick.map'
            if exists(file):
                if d.month==2:
                    mid_day = 13
                else:
                    mid_day = 14
                dates_u.append(d + relativedelta(days=mid_day))
            #if it does append dates_u
        ### now work over the date and adjust for summer
        ### remove all months = [5,6,8,7,9]
        ### also need hole end points to be at month end not mid
        self.dates= []
        month_keep = [1,2,3,4,10,11,12]
        for d in dates_u:
            if d.month in month_keep:
                if fill_end_months and d.month == 4:
                    self.dates.append(d)
                    d_end = d.replace(day=30)
                    self.dates.append(d_end)
                elif fill_end_months and d.month == 10:
                    d_start = d.replace(day=1)
                    self.dates.append(d_start)
                    self.dates.append(d)
                else:
                    self.dates.append(d)
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')



class Kimura():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath):
        self.name = 'Kimura_drift'
        self.path = ppath
        
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        dates =[]
        d_no = (time_end-time_start).days +3 
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
            infile = self.path+d.strftime('%y%m%d')+".amsr36i"
            if exists(infile):
                dates.append(d)
            else:
                infile = self.path+d.strftime('%y%m%d')+".amsr18i"
                if exists(infile):
                    dates.append(d)
        self.dates = dates
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')

    # daily points in yearly files

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self,dates_u,verbos=False):
        dimY = 145
        dimX = 145
        d_no = np.shape(dates_u)[0]
        data_u = np.zeros([d_no,dimX,dimY])
        data_v = np.zeros([d_no,dimX,dimY])
        for n,d in enumerate(dates_u):
            infile = self.path+d.strftime('%y%m%d')+".amsr36i"
            if not exists(infile):
                infile = self.path+d.strftime('%y%m%d')+".amsr18i"
            if verbos: print(infile)
#             print(d.strftime('%y%m%d'))
            with open(infile, 'rb') as fr:
                hdr = fr.read(0)
                ice = np.fromfile(fr, dtype=np.float32)
            ice = ice/100
            ice[ice>9] = np.nan
            iceG = ice.reshape(dimX,dimY,2)
            data_u[n] =-iceG[:,:,1]
            data_v[n] = iceG[:,:,0]
        return data_u,data_v

class CICE_jas:
    def __init__(self,ppath):
        self.name = 'CICE_1deg'
        self.path = ppath        
        self.vyear_load = 0
        self.hyear_load = 0
        self.ayear_load = 0
        self.vels_loaded = False
        self.hi_loaded = False
        self.aice_loaded = False

    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        dates =[]
        n_yrs = (time_end.year - time_start.year)+1
        for y in range(n_yrs):
            yu = time_start.year + y
            d0 = dt.datetime(yu,1,1)
            f_name = 'cice_daily_'+str(yu)+'.nc'
            if exists(self.path+f_name):
                f_nc = Dataset(self.path+f_name)
                [dates.append(d0 + relativedelta(days = d)) 
                     for d in range(f_nc['time'].shape[0])]
                f_nc.close()
        self.dates = dates
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')

    def get_vels(self,dates_u,verbos=False):
        d0 = dt.datetime(1970,1,1)
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            dates = dates_u
            single_year = True
        else: ## multiple years
            dates = [d for d in dates_u if d.year == dates_u[0].year]
            single_year = False
            # one year, one file
        yu = dates_u[0].year
        if ((self.vyear_load != yu) or (not self.vels_loaded)):
            print('loading new year of data: '+str(yu))
            f_name = 'cice_daily_'+str(yu)+'.nc'
            f_nc = Dataset(self.path+f_name)
    #         print(p0,p1)
            self.u = f_nc['uvel_d'][:]
            self.v = f_nc['vvel_d'][:]
            self.u[self.u.mask] = np.nan
            self.v[self.v.mask] = np.nan
            f_nc.close()
            self.vyear_load = yu
            self.vels_loaded= True
        p0 = dates[ 0].timetuple().tm_yday -1
        p1 = dates[-1].timetuple().tm_yday 
        datau = self.u[p0:p1,:,:].transpose((0,2,1))
        datav = self.v[p0:p1,:,:].transpose((0,2,1))
        if single_year:
            return datau,datav
        else: ## multiple years
            data1 = datau
            dates = [d for d in dates_u if d.year == dates_u[-1].year]
            yu = dates[-1].year
            print('loading new year of data: '+str(yu))
            f_name = 'cice_daily_'+str(yu)+'.nc'
            f_nc = Dataset(self.path+f_name)
    #         print(p0,p1)
            self.u = f_nc['uvel_d'][:]
            self.v = f_nc['vvel_d'][:]
            self.u[self.u.mask] = np.nan
            self.v[self.v.mask] = np.nan
            f_nc.close()
            self.hyear_load = yu
            self.hi_loaded= True
            p0 = dates[ 0].timetuple().tm_yday -1
            p1 = dates[-1].timetuple().tm_yday 
            datau2= self.u[p0:p1,:,:].transpose((0,2,1))
            datav2= self.v[p0:p1,:,:].transpose((0,2,1))
            datau = np.vstack([datau,datau2])
            datav = np.vstack([datav,datav2])
            return datau,datav

    def get_aice(self,dates_u,verbos=False):
        d0 = dt.datetime(1970,1,1)
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            dates = dates_u
            single_year = True
        else: ## multiple years
            dates = [d for d in dates_u if d.year == dates_u[0].year]
            single_year = False
            # one year, one file
        yu = dates_u[0].year
        if ((self.ayear_load != yu) or (not self.aice_loaded)):
            print('loading new year of data: '+str(yu))
            f_name = 'cice_daily_'+str(yu)+'.nc'
            f_nc = Dataset(self.path+f_name)
    #         print(p0,p1)
            self.a = f_nc['aice_d'][:]
            self.a[self.a.mask] = np.nan
            f_nc.close()
            self.ayear_load = yu
            self.aice_loaded= True
        p0 = dates[ 0].timetuple().tm_yday -1
        p1 = dates[-1].timetuple().tm_yday 
        datau = self.a[p0:p1,:,:].transpose((0,2,1))
        if single_year:
            return datau
        else: ## multiple years
            data1 = datau
            dates = [d for d in dates_u if d.year == dates_u[-1].year]
            yu = dates[-1].year
            print('loading new year of data: '+str(yu))
            f_name = 'cice_daily_'+str(yu)+'.nc'
            f_nc = Dataset(self.path+f_name)
    #         print(p0,p1)
            self.a = f_nc['aice_d'][:]
            self.a[self.a.mask] = np.nan
            f_nc.close()
            self.ayear_load = yu
            self.aice_loaded= True
            p0 = dates[ 0].timetuple().tm_yday -1
            p1 = dates[-1].timetuple().tm_yday 
            data2 = self.a[p0:p1,:,:].transpose((0,2,1))
            datau = np.vstack([data1,data2])
            return datau

    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) == 0:
            dates = dates_u
            single_year = True
        else: ## multiple years
            dates = [d for d in dates_u if d.year == dates_u[0].year]
            single_year = False
            # one year, one file
        yu = dates_u[0].year
        if ((self.hyear_load != yu) or (not self.hi_loaded)):
            print('loading new year of data: '+str(yu))
            f_name = 'cice_daily_'+str(yu)+'.nc'
            f_nc = Dataset(self.path+f_name)
    #         print(p0,p1)
            self.h = f_nc['hi_d'][:]
            self.h[self.h.mask] = np.nan
            f_nc.close()
            self.hyear_load = yu
            self.hi_loaded= True
        p0 = dates[ 0].timetuple().tm_yday -1
        p1 = dates[-1].timetuple().tm_yday 
        datau = self.h[p0:p1,:,:].transpose((0,2,1))
        if single_year:
            return datau
        else: ## multiple years
            data1 = datau
            dates = [d for d in dates_u if d.year == dates_u[-1].year]
            yu = dates[-1].year
            print('loading new year of data: '+str(yu))
            f_name = 'cice_daily_'+str(yu)+'.nc'
            f_nc = Dataset(self.path+f_name)
    #         print(p0,p1)
            self.h = f_nc['hi_d'][:]
            self.h[self.h.mask] = np.nan
            f_nc.close()
            self.hyear_load = yu
            self.hi_loaded= True
            p0 = dates[ 0].timetuple().tm_yday -1
            p1 = dates[-1].timetuple().tm_yday 
            data2 = self.h[p0:p1,:,:].transpose((0,2,1))
            datau = np.vstack([data1,data2])
            return datau
        

class Bristol_thickness:
    def __init__(self,ppath,var='Sea_Ice_Thickness_incSMOS'):
        """
        var can be any of
        'Sea_Ice_Thickness'
        'Sea_Ice_Thickness_W99'
        'Sea_Ice_Thickness_incSMOS'
        'Sea_Ice_Volume'
        'Sea_Ice_Volume_W99'
        'Sea_Ice_Volume_incSMOS'
        """
        self.name = 'Bristol_hi'
        self.path = ppath        
        self.hi1 = 0.0
        self.dl1 = None
        self.hi2 = 0.0
        self.dl2 = None
        self.var = var

    def get_dates(self,time_start,time_end,fill_end_months=False):
        """
        returns the all encompassing date list for use with the forcing object
        """
        blurb ='ubristol_cryosat2_seaicethickness_nh25km_'
        dy = time_end.year-time_start.year
        dm = time_end.month-time_start.month
        m_no = dy*12 + dm +2
        # make sure we get the bracket points
        ts_m = dt.datetime(time_start.year,time_start.month,1)
        dates_u = []
        for mn in range(m_no):
            d = ts_m+ relativedelta(months = mn )
            file = self.path+blurb+d.strftime('%Y_%m_')+'v1.nc'
            if exists(file):
                if d.month==2:
                    mid_day = 13
                else:
                    mid_day = 14
                dates_u.append(d + relativedelta(days=mid_day))
            #if it does append dates_u
        ### now work over the date and adjust for summer
        ### remove all months = [5,6,8,7,9]
        ### also need hole end points to be at month end not mid
        self.dates= []
        month_keep = [1,2,3,4,10,11,12]
        for d in dates_u:
            if d.month in month_keep:
                if fill_end_months and d.month == 4:
                    self.dates.append(d)
                    d_end = d.replace(day=30)
                    self.dates.append(d_end)
                elif fill_end_months and d.month == 10:
                    d_start = d.replace(day=1)
                    self.dates.append(d_start)
                    self.dates.append(d)
                else:
                    self.dates.append(d)
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')


    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        ### check dates
        ### normal month, data on the first
        blurb ='ubristol_cryosat2_seaicethickness_nh25km_'
        d  = dates_u[0]
        dn = np.shape(dates_u)[0]
        file = self.path+blurb+d.strftime('%Y_%m_')+'v1.nc'
        f_nc = Dataset(file)
        hi = f_nc[self.var][:]
        dx,dy = hi.shape 
        data = np.empty([dn,dx,dy])
        data[0] = hi
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            dload = d.replace(day=1)
            file = self.path+blurb+d.strftime('%Y_%m_')+'v1.nc'
            f_nc = Dataset(file)
            hi = f_nc[self.var][:]
            data[n+1] = hi
            f_nc.close()
            #if it does append dates_u
        return data

class ESA_CCI:
    def __init__(self,ppath,satellite):
        """
        satellite can be:
        'CRYOSAT2'
        'ENVISAT'
        """
        self.name = 'ESA_CCI_'+satellite
        self.hi1 = 0.0
        self.dl1 = None
        self.hi2 = 0.0
        self.dl2 = None
        if satellite == 'CRYOSAT2':
            self.sat_blurb = 'ESACCI-SEAICE-L3C-SITHICK-SIRAL_CRYOSAT2-NH25KMEASE2-'
            self.path = ppath+'CS2_CCI/'
        elif satellite == 'ENVISAT':
            self.sat_blurb = 'ESACCI-SEAICE-L3C-SITHICK-RA2_ENVISAT-NH25KMEASE2-'
            self.path = ppath+'Envisat_CCI/'

    def get_dates(self,time_start,time_end,fill_end_months=False):
        """
        returns the all encompassing date list for use with the forcing object
        """
        dy = time_end.year-time_start.year
        dm = time_end.month-time_start.month
        m_no = dy*12 + dm +2
        # make sure we get the bracket points
        ts_m = dt.datetime(time_start.year,time_start.month,1)
        dates_u = []
        for mn in range(m_no):
            d = ts_m+ relativedelta(months = mn )
            file = self.path+d.strftime('%Y/')+self.sat_blurb+d.strftime('%Y%m')+'-fv2.0.nc'
            if exists(file):
                if d.month==2:
                    mid_day = 13
                else:
                    mid_day = 14
                dates_u.append(d + relativedelta(days=mid_day))
            #if it does append dates_u
        ### now work over the date and adjust for summer
        ### remove all months = [5,6,8,7,9]
        ### also need hole end points to be at month end not mid
        self.dates= []
        month_keep = [1,2,3,4,10,11,12]
        for d in dates_u:
            if d.month in month_keep:
                if fill_end_months and d.month == 4:
                    self.dates.append(d)
                    d_end = d.replace(day=30)
                    self.dates.append(d_end)
                elif fill_end_months and d.month == 10:
                    d_start = d.replace(day=1)
                    self.dates.append(d_start)
                    self.dates.append(d)
                else:
                    self.dates.append(d)
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')


    def get_hi(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        ### check dates
        ### normal month, data on the first
        d  = dates_u[0]
        dn = np.shape(dates_u)[0]
        file = self.path+d.strftime('%Y/')+self.sat_blurb+d.strftime('%Y%m')+'-fv2.0.nc'
        f_nc = Dataset(file)
        hi = f_nc['sea_ice_thickness'][0].T
        dx,dy = hi.shape 
        data = np.empty([dn,dx,dy])
        data[0] = hi
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            dload = d.replace(day=1)
            file = self.path+d.strftime('%Y/')+self.sat_blurb+d.strftime('%Y%m')+'-fv2.0.nc'
            f_nc = Dataset(file)
            hi = f_nc['sea_ice_thickness'][0].T
            data[n+1] = hi
            f_nc.close()
            #if it does append dates_u
        return data


class Bristol_thickness_seasonal:
    def __init__(self,ppath,var='Sea_Ice_Thickness',
        infile = 'ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc'
                ):
        """
        var is
        'Sea_Ice_Thickness'
        """
        self.name = 'Bristol_hi_seasonal'
        self.path = ppath        
        self.var = var
        self.file = self.path+infile

    def get_dates(self,time_start,time_end,fill_end_months=False):
        """
        use the Time variable, single file sothis allows for access later 
        """
        
        self.f_nc = Dataset(self.file)
        self.time_vec = self.f_nc.variables['Time'][:]
        self.dates = [dt.datetime(1,1,1)+relativedelta(days=t) for t in self.time_vec]
        self.dates = [t+relativedelta(years = -1) for t in self.dates]
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')

    def get_hi(self,dates_u,verbos=False):
        """
        use the self.dates to find the index of the time in question
        """
        ### find the indices
        idx = [np.argwhere(np.array([d == du for d  in self.dates] ))[0,0] 
                                             for du in dates_u]
        if verbos:
            for i,du in zip(idx,dates_u):
                dcheck = dt.datetime(1,1,1)
                dcheck = dcheck+relativedelta(days=self.time_vec[i])
                dcheck = dcheck+relativedelta(years=-1)
                print(du.strftime('%Y%m%d-')+dcheck.strftime('%Y%m%d'))
        hi = self.f_nc.variables[self.var][idx]
        return hi.data

    def get_err(self,dates_u,verbos=False):
        ## errs need to be dimensional
        """
        use the self.dates to find the index of the time in question
        """
        ### find the indices
        idx = [np.argwhere(np.array([d == du for d  in self.dates] ))[0,0] 
                                             for du in dates_u]
        if verbos:
            for i,du in zip(idx,dates_u):
                dcheck = dt.datetime(1,1,1)
                dcheck = dcheck+relativedelta(days=self.time_vec[i])
                dcheck = dcheck+relativedelta(years=-1)
                print(du.strftime('%Y%m%d-')+dcheck.strftime('%Y%m%d'))
        hi = self.f_nc.variables['Sea_Ice_Thickness_Uncertainty'][idx]
        return hi.data


class ECCC_drift:
    def __init__(self,ppath):
        self.name = 'ECCC_Drift'
        self.path = ppath
        self.files = glob.glob(ppath+'*/*.nc')
        self.files.sort()
    
    def get_dates(self,time_start,time_end):
        dates =[]
        ##### date in the filename
        ### 2020-03/ECCC_SIM_PREMELT_20200301_20200303_v1.0.nc
        #### list the files
        for f in self.files:
            if 'PREMELT_' in f:
                #### extract the dates
                f_date_str= f.split('PREMELT_')[1].split('_')[0]
                f_date = dt.datetime.strptime(f_date_str,'%Y%m%d')+relativedelta(days=1)
                if f_date >=time_start and f_date <=time_end:
                    dates.append(f_date)
        self.dates = dates
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')
                
            
    def get_vels(self,dates_u,verbos=False):
        ##### date in the filename
        ### ECCC_SIM_PREMELT_20200301_20200303_v1.0.nc
        ### 2020-03/ECCC_SIM_PREMELT_20200301_20200303_v1.0.nc
        xvel = []
        yvel = []
        for d in dates_u:
            t1 = (d - relativedelta(days=1)).strftime('%Y%m%d')
            t2 = (d + relativedelta(days=1)).strftime('_%Y%m%d')
            f_name = (d - relativedelta(days=1)).strftime('%Y-%m/')+'ECCC_SIM_PREMELT_'+t1+t2+'_v1.0.nc'
            f_nc = Dataset(self.path+f_name)
            x =-f_nc.variables['v'][0]
            y =-f_nc.variables['u'][0]
            #### convert km/day to m/s
            xvel.append(x*1e3/(24*60*60))
            yvel.append(y*1e3/(24*60*60))
        return np.array(xvel),np.array(yvel)
            

class Bremmen_IC:
    def __init__(self,ppath,reduceECCC=False):
        self.name = 'Bremmen_IC'
        self.path = ppath
        self.reduceECCC = reduceECCC
    
    def get_dates(self,time_start,time_end):
        dates =[]
        ##### date in the filename
        ### sic_modis-aqua_amsr2-gcom-w1_merged_nh_1000m_20200116.nc
        #### list the dirs
        #### first get the years we care for
        ystart = time_start.year
        yend   = time_end.year
        for y in range(ystart,yend+1):
            files = glob.glob(self.path+str(y)+'/sic_modis*.nc')
            files.sort()
            for f in files:
                if 'sic_modis-aqua' in f:
                    #### extract the dates
                    f_date_str= f.split('merged_nh_1000m_')[1].split('.')[0]
                    f_date = dt.datetime.strptime(f_date_str,'%Y%m%d')
                    if f_date >=time_start and f_date <=time_end:
                        dates.append(f_date)
        self.dates = dates
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')
                
            
    def get_aice(self,dates_u,verbos=False):
        ##### date in the filename
        ### sic_modis-aqua_amsr2-gcom-w1_merged_nh_1000m_20200116.nc
        IC_list = []
        for d in dates_u:
            t1 = d.strftime('%Y%m%d')
            f_name = 'sic_modis-aqua_amsr2-gcom-w1_merged_nh_1000m_'+t1+'.nc'
            f_name = d.strftime('%Y/')+f_name
            f_nc = Dataset(self.path+f_name)
            if self.reduceECCC:
                GBx = slice(750,2700)
                GBy = slice(1750,4000)
                try:
                    aice = f_nc.variables['mersic'][GBx,GBy]
                    aice.mask[:] = False
                except KeyError:
                    aice = f_nc.variables['sic_merged']
                    aice = np.flipud(aice)
                    aice = aice[GBx,GBy]
            else:
                try:
                    aice = f_nc.variables['mersic'][:]
                    aice.mask[:] = False
                except KeyError:
                    aice = f_nc.variables['sic_merged'][:]
                    aice = np.flipud(aice)
            IC = aice/100
            IC[aice>100] = np.nan
            IC_list.append(IC)
        if len(IC_list) == 1:
            return  np.expand_dims(IC_list[0],axis=0)
        else:
            return  np.array(IC_list)

class Bristol_CPOM_daily:
    def __init__(self,ppath):
        self.name = 'Bristol_CPOM_daily_thickness'
        self.path = ppath
        self.hi_loaded = False
        self.file_start = None
    
    def get_dates(self,time_start,time_end):
        dates =[]
        ##### date in the filename
        #### cryosat2_cpomOIrfb_seaicethickness_nh_daily_25km_2019-2020_v1.nc
        blurb = 'cryosat2_cpomOIrfb_seaicethickness_nh_daily_25km_'
        #### list the dirs
        #### first get the years we care for
        ystart = time_start.year
        yend   = time_end.year
        for y in range(ystart,yend+1):
            # file per season, get time_start year
            files = glob.glob(self.path+blurb+str(y)+'*.nc')
            if len(files) == 1:
                ## open file and get Time array
                fnc = Dataset(files[0])
                t_vec = fnc.variables['Time'][:].data
                dt_vec = [dt.datetime(1,1,1)+relativedelta(days=t) for t in t_vec]
                dt_vec = [t+relativedelta(years = -1) for t in dt_vec]
                dt_vec = [t for t in dt_vec if t>=time_start and t <=time_end]
                dates.extend(dt_vec)
        self.dates = dates
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')
                
            
    def get_hi(self,dates_u,verbos=False):
        ##### date in the filename
        hi_list = []
        for d in dates_u:
            if self.file_start is not None:
                ddiff = (d-self.file_start).days
            if self.hi_loaded and ddiff < 300:
                ### get the data from the array
                hi_list.append(self.hi[ddiff])
            else:
                #### load a new array
                if d.month > 7: y = d.year
                else: y = d.year - 1
                print('loading new year of data: '+str(y))
                blurb = 'cryosat2_cpomOIrfb_seaicethickness_nh_daily_25km_'
                files = glob.glob(self.path+blurb+str(y)+'*.nc')
                fnc = Dataset(files[0])
                ### get hi array
                hi = fnc.variables['Sea_Ice_Thickness'][:]
                hi[hi>100] = np.nan
                self.hi = hi
                ### get file_start
                t = fnc.variables['Time'][0]
                tdt = dt.datetime(1,1,1)+relativedelta(days=t)
                tdt = tdt+relativedelta(years = -1)
                self.file_start = tdt
                self.hi_loaded = True
                ddiff = (d-self.file_start).days
                hi_list.append(self.hi[ddiff])
        return  np.array(hi_list)





class CPOM_FB_interp:
    def __init__(self,ppath):
        self.name = 'CPOM_FB_interp'
        self.path = ppath
    
    def get_dates(self,time_start,time_end):
        dates =[]
        ##### date in the filename
        ### FB_interp_2011-2012_50km_20111105.npy ## beginnig of the season
        ### FB_interp_2011-2012_50km_20120328.npy ## end of season
        blurb = 'FB_interp_'
        #### first get the files
        dates = []
        files = glob.glob(self.path+blurb+'*.npy')
        for f in files:
            ftime = f.split('50km_')[1].split('.npy')[0]
            t = dt.datetime.strptime(ftime,'%Y%m%d')
            ## open file and get Time array
            if t>=time_start and t <=time_end: dates.append(t)
        dates.sort()
        self.dates = dates
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')
                
            
    def get_hi(self,dates_u,verbos=False):
        ##### date in the filename
        blurb = 'FB_interp_'
        ### FB_interp_2011-2012_50km_20111105.npy ## beginnig of the season
        ### FB_interp_2011-2012_50km_20120328.npy ## end of season
        hi_list = []
        for d in dates_u:
            ### get the season numbers
            if d.month > 6: ### end of year
                season_str = str(d.year)+'-'+str(d.year+1)
            elif d.month < 6: ### beginning of year
                season_str = str(d.year-1)+'-'+str(d.year)
            fname = self.path+blurb+season_str+'_50km_'+d.strftime('%Y%m%d')+'.npy'
            hi = np.load(fname)
            hi_list.append(hi)
        return  np.array(hi_list)


class AVHRR:
    def __init__(self,ppath):
        """
        var can be any of
        cdr_surface_temperature
        cdr_surface_albedo
         cdr_surface_downwelling_shortwave_flux
         cdr_surface_downwelling_longwave_flux
         cdr_surface_upwelling_shortwave_flux
         cdr_surface_upwelling_longwave_flux
         cdr_sea_ice_thickness
         cdr_toa_net_downward_shortwave_flux
         cdr_toa_outgoing_shortwave_flux
         cdr_toa_outgoing_longwave_flux
         surface_type
          cloud_particle_phase
          cloud_particle_radius
           cloud_optical_depth
           cloud_top_pressure
           cloud_top_temperature
           cloud_type
           surface_shortwave_cloud_radiative_forcing
            surface_longwave_cloud_radiative_forcing
        SET in get_var('variable_name',dates)
        """
        self.name = 'Pathfinder_extended'
        self.path = ppath        
    
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        ### file per day, dir per year, set this first then search
        yrng = range(time_start.year,time_end.year+1)
        dates_u = []
        for y in yrng:
            ### all files for that year
            flist = glob.glob(self.path+str(y)+'/Polar-APP-X_*.nc')
            flist.sort()
            for f in flist:
                ### extract date from file name
                ftstr = f.split('/'+str(y)+'/')[1].split('_d')[1].split('_')[0]
                fdt = dt.datetime.strptime(ftstr,'%Y%m%d')
                if fdt>=time_start and fdt<=time_end:
                    dates_u.append(fdt)
        self.dates=dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')
    
    def get_var(self,var,dates_u,verbos=False):
        ### contruct filename from dates
        a_out = []
        for d in dates_u:
            file_search = self.path+d.strftime('%Y/Polar-APP-X_*%Y%m%d*')
            file = glob.glob(file_search)[0]
            if verbos: print(file)
            f_nc = Dataset(file)
            array = np.fliplr(f_nc[var][0].T)
            array[array>9998] = np.nan
            a_out.append(array)
            f_nc.close()
        return np.array(a_out)
    
    

class SM_LG:
    def __init__(self,ppath,forcing = 'ERA5'):
        """
        Lagrangian SNOW model - set the forcing as well as path
        
        forcing can be 'ERA5' or 'MERRA2'
        """
        self.name = '_'.join(['SM-LG',forcing])
        self.path = ppath        
        self.forcing = forcing        
        self.file = glob.glob(ppath+'*'+forcing+'*.nc')[0]
        self.dt0 = dt.datetime(1980,1,1)
        self.t0 = 17352696 ## this is to be true to filename
        ### self.t0 = 17352696 ## - 48 ### this is to be true to attribute. 
        
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        ### one file - time is hours since 1-1-1 00:00:00 (yikes) - 2days out
        ### or if the filename is to be believed  t-17352696 = hours since 1980-8-1
        f_nc = Dataset(self.file)
        self.t_vec = f_nc['time'][:]-self.t0
        self.dtime =[self.dt0+relativedelta(hours=t) for t in self.t_vec] 
        self.dates = [d for d in self.dtime if d>=time_start and d<=time_end]
        f_nc.close()
        print(self.name+' Found '+str(np.shape(self.dates)[0])+' dates')


    def get_snod(self,dates_u,verbos=False):
        ### contruct filename from dates
        if verbos: print([d.strftime('%Y-%m-%d') for d in dates_u])
        tindx = [(d-self.dt0).days*24 for d in dates_u ]
        indx = [np.abs(self.t_vec - t).argmin() for t in tindx]
        if verbos: print(indx)
        f_nc = Dataset(self.file)
        snod = f_nc['snod'][indx].transpose(0,2,1)
        snod[snod.mask] = np.nan
        return snod.data
        f_nc.close()