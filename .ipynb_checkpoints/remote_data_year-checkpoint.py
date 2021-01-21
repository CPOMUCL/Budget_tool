# converts a remote bunch of data files into a data year
import numpy as np
import datetime as dt
import copy
from scipy import stats
from dateutil.relativedelta import relativedelta
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from sys import path
path.insert(0, '/Users/H/WAVES/geo_data_group/')
import data_year as dy

class remote_data_year(dy.data_year):
    # data year inits using
    ### data, dates, periods
    
    ### the remote data files we use are
    ### get_dates
    ### get_arrays
    ### initalise using a tstart and twindow
    ### and an inputs_objects
    
    
    def __init__(self,INPUT,periods):
        
        ### need a dates list - get_dates and and coverter
        ### periods can be asscessed from this
        
        ### need to convert INPUT object to have an __get_items__ methods to align with dates
        ### see below with rdy_input
        ### need to define the get_data method for INPUT
        
        rdyi = rdy_input(INPUT)
        
        ### get periods from dates list
        ## actually don't bother now
        
        super().__init__(rdyi,INPUT.dates,periods)
        

class remote_vec_data_year(dy.vec_data_year):

    
    def __init__(self,INPUT,periods):
        
        ### need a dates list - get_dates and and coverter
        ### periods can be asscessed from this
        
        ### need to convert INPUT object to have an __get_items__ methods to align with dates
        ### see below with rdy_input
        ### need to define the get_data method for INPUT
        
        rdyx = rvdy_xinput(INPUT)
        rdyy = rvdy_yinput(INPUT)
        
        ### get periods from dates list
        ## actually don't bother now
        
        super().__init__(rdyx,rdyy,INPUT.dates,periods)
        

class rdy_input():
    def __init__(self,Dobj):
        ###
        self.Dobj = Dobj
        d= [self.Dobj.dates[0]]
        array = self.Dobj.get_data(d)
        t_p,m,n = np.shape(array)
        t_p = np.shape(self.Dobj.dates)[0]
        self.shape = (t_p,m,n)
    

        
    def __getitem__(self,indx):
        if (type(indx) == int) or (type(indx) == np.int64):
            t_p = indx
            m = slice(None)
            n = slice(None)
        else:
            t_p,m,n = indx
        # can we convert an a:b version of t_p
        # to a list of dates??
        d = self.Dobj.dates[t_p]
        if type(d) == dt.datetime:
            d= [d]
        ## call some sort of data collection on d 
        array = self.Dobj.get_data(d)
        if np.shape(d)[0] ==1:
            return array[0,m,n]
        else:
            return array[:,m,n]
    

class rvdy_input():
    def __init__(self,Dobj):
        ###
        self.Dobj = Dobj
        d= [self.Dobj.dates[0]]
        array,y = self.Dobj.get_vector(d)
        t_p,m,n = np.shape(array)
        t_p = np.shape(self.Dobj.dates)[0]
        self.shape = (t_p,m,n)
    
        
    def __getitem__(self,indx):
        if (type(indx) == int) or (type(indx) == np.int64):
            t_p = indx
            m = slice(None)
            n = slice(None)
        else:
            t_p,m,n = indx
        # can we convert an a:b version of t_p
        # to a list of dates??
        d = self.Dobj.dates[t_p]
        if type(d) == dt.datetime:
            d= [d]
        ## call some sort of data collection on d 
        xarray,yarray = self.Dobj.get_vector(d)
        if np.shape(d)[0] ==1:
            return xarray[0,m,n], yarray[0,m,n]
        else:
            return xarray[:,m,n], yarray[:,m,n]
    
class rvdy_xinput(rvdy_input):
    def __init__(self,Dobj):
        ###
        super().__init__(Dobj)
        
    def __getitem__(self,indx):
        x,y = super().__getitem__(indx)
        return x
    
class rvdy_yinput(rvdy_input):
    def __init__(self,Dobj):
        ###
        super().__init__(Dobj)
        
    def __getitem__(self,indx):
        x,y = super().__getitem__(indx)
        return y 
