
import numpy as np
import datetime as dt
import struct
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from os.path import exists
from scipy.interpolate import griddata
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree


class Pathfinder():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,grid=False):
        self.name = 'Pathfinder'
        self.path = ppath
        self.vyear_load = 0
        self.vels_loaded = False
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
            f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
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
                f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                f_nc = Dataset(self.path+f_name)
        #         print(p0,p1)
                self.u = f_nc['u'][:]
                self.v = f_nc['v'][:]
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
                    datau[n][self.grid.lats>88] = np.nanmean
                    datav[n][self.grid.lats>88] = np.nanmean
            return datau,datav
    
    
    
class PIOMAS():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,time_smooth = 0):
        self.name = 'PIOMAS'
        self.path = ppath
        self.vyear_load = 0
        self.hyear_load = 0
        self.ayear_load = 0
        self.vels_loaded = False
        self.hi_loaded = False
        self.aice_loaded = False
        self.tsmth = time_smooth
        
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        PIOMAS is a standardised list so we can just build a list
        """
        dates =[]
        n_yrs = (time_end.year - time_start.year)-1
        if n_yrs>-1:
            y0 = dt.datetime(time_start.year,1,1)
            ye = dt.datetime(time_start.year,12,31)
            data_f = self.path+'uiday.H'+y0.strftime('%Y')
            if exists(data_f):
                for d in range(time_start.timetuple().tm_yday-1,
                               ye.timetuple().tm_yday):
                    dates.append(y0 +  relativedelta(days = d))
            for y in range(n_yrs):
                y0 += relativedelta(years=1)
                ye += relativedelta(years=1)
                data_f = self.path+'uiday.H'+y0.strftime('%Y')
                if exists(data_f):
                    for d in range(ye.timetuple().tm_yday):
                        dates.append(y0 +  relativedelta(days = d))
            y0 += relativedelta(years=1)
            ye  = time_end
            data_f = self.path+'uiday.H'+y0.strftime('%Y')
            if exists(data_f):
                for d in range(ye.timetuple().tm_yday):
                    dates.append(y0 +  relativedelta(days = d))
        else:
            y0 = dt.datetime(time_start.year,1,1)
            data_f = self.path+'uiday.H'+y0.strftime('%Y')
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

                self.vels=np.asarray(data).reshape(365,2,120,360)
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

                self.hi=np.asarray(data).reshape(365,120,360)
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

                self.aice=np.asarray(data).reshape(365,120,360)
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
            if d.year<2020:
#                 infile = self.path+"/nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
            if d.year>2019:
#                 infile = self.path+"/nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
#             if d.year<2019:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             if d.year>2018:
#                 infile = self.path+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
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
            if d.year<2020:
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
            if d.year>2019:
                infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            # check infile exists 
            if exists(infile):
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
        if d > dt.datetime(2019,6,1):
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
            if d > dt.datetime(2019,6,1):
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
        if d > dt.datetime(2019,6,1):
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
            if d > dt.datetime(2019,6,1):
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
            if d > dt.datetime(2019,6,1):
                blurb2 = '_r_v202_02_l4sit.nc'
            else:
                blurb2 = '_r_v202_01_l4sit.nc'
            t0 = (d - relativedelta(days = 3)).strftime('%Y%m%d_')
            t1 = (d + relativedelta(days = 3)).strftime('%Y%m%d')
            file = self.path+d.strftime('%Y/%m/')+blurb+t0+t1+blurb2
            if exists(file):
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
    
    """
    def __init__(self,ppath):
        self.path = ppath
        self.name = 'osisaf'
# next function will take a list of dates and return an appropriately orientated arrays
# give a 
    def get_dates(self,time_start,time_end):
        # individual drift for each day slice
        dates_u = []
        d_no = (time_end-time_start).days +3 
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
            t1 = (d - relativedelta(hours=12)).strftime('%Y%m%d%H00')
            t2 = (d + relativedelta(hours=36)).strftime('%Y%m%d%H00')
            # check if were on the last day of the month
            next_month = d.replace(day=28) + dt.timedelta(days=4)  # this will never fail
            next_month.replace(day = 1)
            if (next_month - d).days>1:
                td = d.strftime('%Y/%m/')
            else:
                td = (d + relativedelta(months=1)).strftime('%Y/%m/')
            f_name = td+'ice_drift_nh_polstere-625_multi-oi_'+t1+'-'+t2+'.nc'
    #         print(f_name)
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
        # one year, one file
        t1 = (d - relativedelta(hours=12)).strftime('%Y%m%d%H00')
        t2 = (d + relativedelta(hours=36)).strftime('%Y%m%d%H00')
            # check if were on the last day of the month
        next_month = d.replace(day=28) + dt.timedelta(days=4) 
        next_month.replace(day = 1)
        if (next_month - d).days>1:
            td = d.strftime('%Y/%m/')
        else:
            td = (d + relativedelta(months=1)).strftime('%Y/%m/')
        ## extra date check for v component
        vmult = -1.0
        if d > dt.datetime(2015,9,12): vmult = 1.0
        f_name = td+'ice_drift_nh_polstere-625_multi-oi_'+t1+'-'+t2+'.nc'
        f_nc = Dataset(self.path+f_name)
        unow = np.fliplr(f_nc['dX'][0].T)
        vnow = np.fliplr(f_nc['dY'][0].T)
        unow[unow.mask] = np.nan
        vnow[vnow.mask] = np.nan
        dx,dy = np.shape(unow)
        data_u = np.zeros([d_no,dx,dy])
        data_v = np.zeros([d_no,dx,dy])
        # convert km/48hrs to m/s
        data_u[0] =  unow*1e3/60/60/48
        data_v[0] = vmult*vnow*1e3/60/60/48
        f_nc.close()
        for n,d in enumerate(dates_u[1:]):
            # one year, one file
            t1 = (d - relativedelta(hours=12)).strftime('%Y%m%d%H00')
            t2 = (d + relativedelta(hours=36)).strftime('%Y%m%d%H00')
            # check if were on the last day of the month
            next_month = d.replace(day=28) + dt.timedelta(days=4) 
            next_month.replace(day = 1)
            if (next_month - d).days>1:
                td = d.strftime('%Y/%m/')
            else:
                td = (d + relativedelta(months=1)).strftime('%Y/%m/')
            f_name = td+'ice_drift_nh_polstere-625_multi-oi_'+t1+'-'+t2+'.nc'
            f_nc = Dataset(self.path+f_name)
            unow = np.fliplr(f_nc['dX'][0].T)
            vnow = np.fliplr(f_nc['dY'][0].T)
            unow[unow.mask] = np.nan
            vnow[vnow.mask] = np.nan
            data_u[n+1] =  unow*1e3/60/60/48
            data_v[n+1] = vmult*vnow*1e3/60/60/48
            f_nc.close()
        return data_u,data_v


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
                xpts,ypts  = self.G.mplot(lon,lat)
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
