# here is the class that holds all the data days/months
# it has all the gridding scripts needed
# it will save load all the data/days/months as needed
# rewriting so all grids and in i-j coords

import numpy as np
import pandas as pd
import datetime
import copy
from netCDF4 import Dataset
from numba import jit
from scipy import stats
import data_year as dy
from dateutil.relativedelta import relativedelta
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from scipy import sparse

class grid_set:
# will make one of these at a time point (as a datetime) defined by timestart

    def __init__(self,mplot):
        self.mplot = mplot
        self.proj = True
        self.files = False
        self.saved = False
        self.grid = False
        self.gridinfo = False
        self.masked = False
        self.data = False
  
        
#     def import_regrid
#     # takes another data, on a seperate lon/lat, and regrids into a data/day/month
    

#     def import_regrid_nc
        
#     def import_regrid_vec
#     # takes vector input, a seperate lon/lat, and regrids into a data/day/month
#     # makes use of rotate/regrid/rotate methods to get it all pretty

#     def import_regrid_vec_nc
        


    def set_proj(self,mplot):
         # puts in a projection mplot too
#        # make make sure you keep hold of regridding projs
        self.mplot = mplot
        self.proj = True
    
    def reproject(self,mplot):
            self.xpts, self.ypts = mplot(self.lons,self.lats)
            self.mplot = mplot
            for a in dir(self):
                if a == 'xptp':
                    self.get_ptp()
                    break
        
    def set_grid_lon_lat(self,lons,lats,grid_list = False):
       # creates a grid depending on wanted resolution 
        if self.proj:
            xpts, ypts = self.mplot(lons,lats)
            self.lons = lons
            self.lats = lats
            self.xpts = xpts
            self.ypts = ypts
            if grid_list:
                print("Linear grid list. Following grid_set methods won't apply, though Gs2Gs regridding will")
                print("Zero values set for saving")
                self.dxRes = 1.0
                self.dyRes = 1.0
                self.m = 1
                self.n = 1
                self.ang_c = 1.0
                self.ang_s = 1.0
                self.xdist = 1.0
                self.ydist = 1.0
                self.gridinfo = True
                self.grid = True
            else:
                self.dxRes = np.mean(np.diff(xpts[0,:]))
                self.dyRes = np.mean(np.diff(ypts[:,0]))
                self.m,self.n = np.shape(lons)
                self.grid = True
                print("Got a grid res = ",self.m," x ",self.n)
                print("Note that all grid info is in nx x ny grids, whilst data is in nx x ny")
        else: print("Projection not defined yet, do that first")
            
    def get_ptp(self):
        """
        Generates pts arrays for pcolor and pcolormesh - midpoitns for grid areas
        """
        if self.grid:
            # extend longitude by 2
            xpt_pad = np.pad(self.xpts, ((0,0),(1,0)), 'edge')
            ypt_pad = np.pad(self.ypts, ((1,0),(0,0)), 'edge')

            self.xptp = xpt_pad[:,:-1]+0.5*(np.diff(xpt_pad,axis=1))
            self.yptp = ypt_pad[:-1,:]+0.5*(np.diff(ypt_pad,axis=0))
            

    def set_grid_dxy(self,dxRes,dyRes):
       # creates a grid depending on wanted resolution 
        if self.proj:
            nx = int((self.mplot.xmax-self.mplot.xmin)/dxRes)+1
            ny = int((self.mplot.ymax-self.mplot.ymin)/dyRes)+1
            lons, lats, xpts, ypts = self.mplot.makegrid(ny, nx, returnxy=True)
            self.lons = lons
            self.lats = lats
            self.xpts = xpts
            self.ypts = ypts
            self.dxRes = dxRes
            self.dyRes = dyRes
            self.grid = True
            self.m = nx
            self.n = ny
            print("Got a grid res = ",nx," x ",ny)
            print("Note that all grid info is in nx x ny grids, whilst data is in nx x ny")
        else: print("Projection not defined yet, do that first")

    def set_grid_mn(self,nx,ny):
       # creates a grid depending on wanted no. of points 
        if self.proj:
            lons, lats, xpts, ypts = self.mplot.makegrid(ny, nx, returnxy=True)
            self.lons = lons
            self.lats = lats
            self.xpts = xpts
            self.ypts = ypts
            self.grid = True
            self.dxRes = (self.mplot.xmax-self.mplot.xmin)/(nx - 1)
            self.dyRes = (self.mplot.ymax-self.mplot.ymin)/(ny - 1)
            self.m = nx
            self.n = ny
            print("Got a grid res = ",nx," x ",ny)
        else: print("Projection not defined yet, do that first")

    def get_grid_info(self,av_dist = True, av_ang = True):
       # creates a grid depending on wanted no. of points 
        # print( self.grid and (not self.gridinfo))
        if self.grid and (not self.gridinfo):
            #iterate over the grid to get dimensions and angles
            # first iterate all x dimensions - m-1/n array
            # then  iterate all y dimensions - m/n-1 array
            xdims = np.empty([self.m-1,self.n])
            ydims = np.empty([self.m,self.n-1])
            self.xdist = np.empty([self.m,self.n])
            self.ydist = np.empty([self.m,self.n])
            self.ang_c = np.empty([self.m,self.n])
            self.ang_s = np.empty([self.m,self.n])
            for i in range(self.m):
                for j in range(self.n-1):
                    ydims[i,j] = ellipsoidal_distance(
                        self.lons[i,j  ],self.lats[i,j  ],
                        self.lons[i,j+1],self.lats[i,j+1],deg=True)
            for i in range(self.m-1):
                for j in range(self.n):
                    xdims[i,j] = ellipsoidal_distance(
                        self.lons[i  ,j],self.lats[i  ,j],
                        self.lons[i+1,j],self.lats[i+1,j],deg=True)

            # then average the available distances i-1,i j-1,j
            if av_dist:
                for i in range(self.m):
                    for j in range(self.n):
                        self.xdist[i,j] = np.nanmean(xdims[:i+1,j][-2:])
                        self.ydist[i,j] = np.nanmean(ydims[i,:j+1][-2:])
            else:
                self.xdist[:-1,:] = xdims
                self.xdist[-1,:]  = xdims[-1,:]
                self.ydist[:,:-1] = ydims
                self.ydist[:,-1]  = ydims[:,-1]
            print("Grid distances calculated: ",np.nanmean(self.xdist)," x ",np.nanmean(self.ydist))
                     
            # then  iterate all angles - this is all points plus the extra possible angles
            # pad the lon lat arrays for iteration
            lon_pad = np.pad(self.lons, (1,1), 'linear_ramp', end_values=(np.nan))
            lat_pad = np.pad(self.lats, (1,1), 'linear_ramp', end_values=(np.nan))
            for j in range(self.m):
                for i in range(self.n):
                    # i + angle
                    yPlus_c,yPlus_s = lon_lat_angle(lon_pad[j+1,i+1],lat_pad[j+1,i+1],
                                                    lon_pad[j+2,i+1],lat_pad[j+2,i+1],
                                                    return_trig = True,deg=True)
                    if av_ang:
                        xPlus_c,xPlus_s = lon_lat_angle(lon_pad[j+1,i+1],lat_pad[j+1,i+1],
                                                        lon_pad[j+1,i+2],lat_pad[j+1,i+2],
                                                        return_trig = True,deg=True)
                        xMins_c,xMins_s = lon_lat_angle(lon_pad[j+1,i+1],lat_pad[j+1,i+1],
                                                        lon_pad[j+1,i  ],lat_pad[j+1,i  ],
                                                        return_trig = True,deg=True)
                        yMins_c,yMins_s = lon_lat_angle(lon_pad[j+1,i+1],lat_pad[j+1,i+1],
                                                        lon_pad[j  ,i+1],lat_pad[j  ,i+1],
                                                        return_trig = True,deg=True)
                    # average all the components first checking the orientation
                    # if j == 20 and i ==12:
                        # print([xPlus_c,xMins_c,yPlus_c,yMins_c])
                        # print([xPlus_s,xMins_s,yPlus_s,yMins_s])
                    if av_ang:
                        self.ang_c[j,i] = np.nanmean([-xPlus_s, xMins_s, yPlus_c,-yMins_c])
                        self.ang_s[j,i] = np.nanmean([ xPlus_c,-xMins_c, yPlus_s,-yMins_s])
                        mag = np.hypot(self.ang_c[j,i],self.ang_s[j,i])
                        self.ang_c[j,i] /= mag
                        self.ang_s[j,i] /= mag
                    else:
                        self.ang_c[j,i] =  yPlus_c
                        self.ang_s[j,i] =  yPlus_s
            print('Angles calculated')
            self.gridinfo = True
        else: print("Grid not defined yet, do that first")

    def get_square_points(self):
        """
        makes the xsq,ysq fields that will let you plot on a square grid
        uses np.meshgrid to make location arrasy statring lower left at (0,0)
        """
        self.xsq,self.ysq = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n),indexing = 'ij')
            
    def check_angles(self,point=False,scale=1.0,project = False):
        # return np.hypot of con/sin, min/max and mean
        check_ang = np.hypot(self.ang_c,self.ang_s)**2
        print('mean ='+str(np.nanmean(check_ang)))
        print('max  ='+str(np.nanmax(check_ang)))
        print('min  ='+str(np.nanmin(check_ang)))
        # if a point is given return a vector to north and x positive
        # so it can be plotted on projection
        if (type(point) == list and project):
            # do it using the projection
            i = point[0]
            j = point[1]
            # vector is due up (0,1)
            Out1 = (self.xpts[i,j],self.ypts[i,j])
            # due north (easy)
            xrot = np.array(0.0) #-self.ang_s[i,j]
            yrot = np.array(1.0) # self.ang_c[i,j]
            u,v = self.mplot.rotate_vector(xrot,yrot,self.lons[i,j],self.lats[i,j])
            # vertical on grid (0,1)
            xrot = -self.ang_c[i,j]
            yrot = -self.ang_s[i,j]
#             # horizontal on grid (1,0)
#             xrot = -self.ang_s[i,j]
#             yrot =  self.ang_c[i,j] 
            u1,v1 = self.mplot.rotate_vector(xrot,yrot,self.lons[i,j],self.lats[i,j])
            return u,v,u1,v1,Out1[0],Out1[1]
        elif type(point) == list:
            # returns two normalised vectors
            i = point[0]
            j = point[1]
            # line1 starts at point
            # goes in direction to j+1 (+ve x)
            xvec = self.xpts[i,j+1] - self.xpts[i,j]
            yvec = self.ypts[i,j+1] - self.ypts[i,j]
    #         print(xvec,yvec)
            # angles are between positive x and due north clockwise
#             xrot =  self.ang_c[i,j]*xvec + self.ang_s[i,j]*yvec
#             yrot =  self.ang_c[i,j]*yvec - self.ang_s[i,j]*xvec
            # rotation is -pi/4 + rotation xrot = xyvec, yrot = -xvec
            xrot =  self.ang_c[i,j]*yvec - self.ang_s[i,j]*xvec
            yrot = -self.ang_c[i,j]*xvec - self.ang_s[i,j]*yvec
    #         print(xrot,yrot)
            print(np.rad2deg(np.arctan2(self.ang_s[i,j],self.ang_c[i,j])))
            Out1 = (self.xpts[i,j],self.ypts[i,j])
            Out2 = (Out1[0] + xvec*scale,Out1[1] + yvec*scale)
            Out3 = (Out1[0] + xrot*scale,Out1[1] + yrot*scale)
            # return the list of x,y's needed for plot
            return ([Out1[0],Out2[0]],
                    [Out1[1],Out2[1]]),([Out1[0],Out3[0]],
                    [Out1[1],Out3[1]])
            
            # line2 starts at point 
            # goes in direction - j+1 plus rotation
    def rotate_vectors_to_plot(self,xvec,yvec):
        """
        utilises the ang_c and ang_s arrays along with the associated projection
        """
        # ur,vr will be in lon/lat
#         ur = xvec*self.ang_c + yvec*self.ang_s
#         vr = yvec*self.ang_c - xvec*self.ang_s
        # test
        ur = -yvec*self.ang_c - xvec*self.ang_s
        vr =  xvec*self.ang_c - yvec*self.ang_s
        
        urr,vrr = self.mplot.rotate_vector(ur,vr,self.lons,self.lats)
        return urr,vrr
        
    def save_grid(self,file):
        if self.grid and self.gridinfo:
            # save lat/lon pts 
            np.savez(file,
                lats = self.lats,
                lons = self.lons,
                xpts = self.xpts,
                ypts = self.ypts,
                dxRes = self.dxRes,
                dyRes = self.dyRes,
                m = self.m,
                n = self.n,
                ang_c = self.ang_c,
                ang_s = self.ang_s,
                xdist = self.xdist,
                ydist = self.ydist)
            print("Grid saved in "+file)
        else:
            print("No grid to save - run get_grid_info")

    def save_grid_nc(self,file,notes=''):
        if self.grid and self.gridinfo:
            # save lat/lon pts 
            NC_f = Dataset(file, 'w', format='NETCDF4')
            NC_f.description = 'python grid_set grid file'+notes

            
            # dimensions
            NC_f.createDimension('time', None)
            NC_f.createDimension('x', self.m)
            NC_f.createDimension('y', self.n)

            # variables
#             time = NC_f.createVariable('time', 'f8', ('time',))
            x = NC_f.createVariable('x', 'f4', ('x',))
            y = NC_f.createVariable('y', 'f4', ('y',))
            lons  = NC_f.createVariable('lons', 'f8', ('y', 'x',))
            lats  = NC_f.createVariable('lats', 'f8', ('y', 'x',))
            xpts  = NC_f.createVariable('xpts', 'f8', ('y', 'x',))
            ypts  = NC_f.createVariable('ypts', 'f8', ('y', 'x',))
            ang_c = NC_f.createVariable('ang_c', 'f8',('y', 'x',))
            ang_s = NC_f.createVariable('ang_s', 'f8',('y', 'x',))
            xdist = NC_f.createVariable('xdist', 'f8',('y', 'x',))
            ydist = NC_f.createVariable('ydist', 'f8',('y', 'x',))
            
            NC_f.setncattr_string('dxRes',self.dxRes)
            NC_f.setncattr_string('dyRes',self.dyRes)

            
            x[:] =self.xpts[0,:]
            y[:] =self.ypts[:,0]
            lons[:,:] = self.lons
            lats[:,:] = self.lats
            xpts[:,:] = self.xpts
            ypts[:,:] = self.ypts
            ang_c[:,:] = self.ang_c
            ang_s[:,:] = self.ang_s
            xdist[:,:] = self.xdist
            ydist[:,:] = self.ydist


            NC_f.close()
#             np.savez(file,
#                 lats = self.lats,
#                 lons = self.lons,
#                 xpts = self.xpts,
#                 ypts = self.ypts,
#                 dxRes = self.dxRes,
#                 dyRes = self.dyRes,
#                 m = self.m,
#                 n = self.n,
#                 ang_c = self.ang_c,
#                 ang_s = self.ang_s,
#                 xdist = self.xdist,
#                 ydist = self.ydist)
            print("Grid saved in "+file)
        else:
            print("No grid to save - run get_grid_info")

    def load_grid(self,file,grid_list=False):
        with np.load(file) as npzfile:
            self.lats = npzfile["lats"]
            self.lons = npzfile["lons"]
            self.xpts = npzfile["xpts"]
            self.ypts = npzfile["ypts"]
            self.dxRes = npzfile["dxRes"] 
            self.dyRes = npzfile["dyRes"] 
            self.m = npzfile["m"][()] 
            self.n = npzfile["n"][()]
            self.ang_c = npzfile["ang_c"] 
            self.ang_s = npzfile["ang_s"] 
            self.xdist = npzfile["xdist"] 
            self.ydist = npzfile["ydist"] 
        self.grid = True
        self.gridinfo = True
        if not grid_list:
            self.get_ptp()
        print("Loaded a grid: "+file)

    def check_grid(self):
        # makes sure the projection and loaded grid are consistent
        if self.proj and self.grid and self.gridinfo:
            proj_dim = self.mplot.xmax - self.mplot.xmin
            proj_dim = proj_dim/self.m
            print("Projection av xdim = ",proj_dim)
            print("dxRes              = ",self.dxRes)
            print("xdist av           = ",np.mean(self.xdist))


    def get_grid_mask(self,inflate = 0.0):
        # makes a land mask for each point then inflates by a distance m
        # makes a land mask for each point then inflates by a distance m
        if self.masked:
            print("Already masked, do it again? set mask = False first")
        else:
            self.mask = np.ones([self.m,self.n])
            for i in range(self.m):
                for j in range(self.n):
                    if self.mplot.is_land(self.xpts[i,j],self.ypts[i,j]):
                         self.mask[i,j] = np.nan
            inf_mask = np.ones([self.m,self.n])
            if (inflate>0.0) and self.gridinfo:
                self.mask_inflate = inflate
                for i in range(self.m):
                    for j in range(self.n):
                        if np.isnan(self.mask[i,j]):
                            inf_p = int(inflate/np.hypot(self.xdist[i,j],self.ydist[i,j]))
                            inf_mask[i-inf_p:i+inf_p+1,j-inf_p:j+inf_p+1] = np.nan
                self.mask = inf_mask
            elif self.gridinfo:
                self.mask_inflate = inflate
        self.masked = True
        

    def inflate_mask(self,inflate = 0.0):
        # makes a land mask for each point then inflates by a distance m
        # makes a land mask for each point then inflates by a distance m
        if self.masked and self.gridinfo:
            inf_mask = np.ones([self.m,self.n])
            if (inflate>0.0) and self.gridinfo:
                self.mask_inflate = inflate
                for i in range(self.m):
                    for j in range(self.n):
                        if np.isnan(self.mask[i,j]):
                            inf_p = int(inflate/np.hypot(self.xdist[i,j],self.ydist[i,j]))
                            inf_mask[i-inf_p:i+inf_p+1,j-inf_p:j+inf_p+1] = np.nan
                self.mask = inf_mask
            elif self.gridinfo:
                self.mask_inflate = inflate
        else:
            print("Not masked so can't inflate")
        

    
    def mask_point(self,lon,lat,inflate = 0):
        x,y = np.unravel_index(np.argmin(
        np.abs(self.lons - lon) + 
        np.abs(self.lats - lat)),
        np.shape(self.lons))
        if (inflate>0.0) and self.gridinfo:
            inf_p = int(inflate/np.hypot(self.xdist[x,y],self.ydist[x,y]))
            self.mask[x-inf_p:x+inf_p+1,y-inf_p:y+inf_p+1] = np.nan
        else:
            self.mask[x,y] = np.nan
            

    def save_mask(self,file):
        if self.masked:
            # save lat/lon pts 
            np.savez(file,
                mask = self.mask,
                mask_inflate = self.mask_inflate,
                m = self.m,
                n = self.n)
            print("Mask saved in "+file)
        else:
            print("No mask to save - run get_grid_mask")


    def load_mask(self,file):
        if self.masked:
            print("Masked already!")
        elif self.gridinfo:
            # save lat/lon pts 
            with np.load(file) as npzfile:
                self.mask = npzfile["mask"]
                self.mask_inflate = npzfile["mask_inflate"]
                m_check = npzfile["m"] 
                n_check = npzfile["n"] 
            if (m_check == self.m)&(n_check == self.n):
                print("Loaded mask, ",m_check," x ",n_check," inflated by ",self.mask_inflate)
                self.masked = True
            else: 
                print("Gird and mask dimensins inconsistent, check them") 
                print("Mask",m_check," x ",n_check," Grid, ",self.m," x ",self.n)
                
    def generate_mask_lonlat(self,lon_r,lat_r,add_mask = True,out='bool'):
        """
        give a lon_r = [l1,l2] range of lons to keep within the mask
        give a lat_r = [l1,l2] range of lats to keep within the mask
        makes a np array that keeps the given range unmaksed.
        add_mask keeps the new mask true to the original GS mask, ie, keeps a land mask
        out = 'bool' makes the out array a logical, T = unmaskes, F = masked
        out = 'float' makes the out array a float, 1.0 unmasked, np.nan = masked
        """
        new_mask = np.ones_like(self.mask)
        new_mask[self.lats.T<lat_r[0]] =np.nan
        new_mask[self.lats.T>lat_r[1]] =np.nan
        new_mask[self.lons.T<lon_r[0]] =np.nan
        new_mask[self.lons.T>lon_r[1]] =np.nan
        if add_mask:
            new_mask[np.isnan(self.mask)] =np.nan
        if out == 'Float':
            return new_mask
        elif out == 'bool':
            out_mask = np.ones_like(self.mask,dtype=bool)
            out_mask[np.isnan(new_mask)] = False
            return out_mask
        
def read_nc_single(ncfile,grid_set,lonlatk,valk,fill_lonlat = False):
    """
    # read and grids, then regrids a single data slice netcdf
    # data array
    # slightly flexible,
    # lonlatk = ['x','y'], say, gives the variable names for lon/lat
    # valkk = ['data'] gives the variable name for the data you want 
    # fill_lonlat = True, allows the code to deal with an even lon/lat grid
    # where the lon lat dimensions are given by a single dimension array
    """
    data_nc = Dataset(ncfile)
    lons = data_nc.variables[lonlatk[0]][:]
    lats = data_nc.variables[lonlatk[1]][:]
    d_array = data_nc.variables[valk[0]][:]
    # if the lat_lons need filling - do it
    if fill_lonlat:
        lon_a,lat_a = np.meshgrid(lons,lats)
    else:
        lon_a = lons
        lat_a = lats
    # regrid depending upon m and grid
    x_nc, y_nc = grid_set.mplot(lon_a.data, lat_a.data)
    new_d_array = griddata((x_nc[~d_array.mask].ravel(), y_nc[~d_array.mask].ravel()),
                d_array[~d_array.mask].ravel(), (grid_set.xpts, grid_set.ypts),
                method='linear')
    return new_d_array


def geo_gradient(array,grid_set):
    """
    gradient function that will take the grid info from the 
    grid_set type class to get gradients 
    the array has to be consistent with the grid set class so it can access the x/ydist parameters
    """
    # check if grid_set has grid info
    if not grid_set.gridinfo:
        print("No grid_set geo grid info - no result")
        return False, False
    in_mn = np.shape(array)
    if in_mn[0]!=grid_set.m or in_mn[1]!=grid_set.n :
        print("input array or geo grid_set not consistently shaped")
        return False, False
    else:
        out_Dax = np.empty_like(array) 
        out_Day = np.empty_like(array) 
        # np gradient can't do an uneven array
        # so we iterate along the columns, then the rows 
        # taking the gradient eachmtime
        # 1 . columns
        for i in range(grid_set.m):
            temp_space = [np.sum(grid_set.ydist[0:j+1,i])
                          for j in range(grid_set.n)]
            out_Day[i,:] = np.gradient(
            array[i,:],temp_space)
        # 2 . rows
        for j in range(grid_set.n):
            temp_space = [np.sum(grid_set.xdist[j,0:i+1])
                          for i in range(grid_set.m)]
            out_Dax[:,j] = np.gradient(
            array[:,j],temp_space)
        return out_Dax,out_Day

def geo_curl(u,v,grid_set):
    """
    curl function that will take the grid info from the 
    grid_set type class to get gradients 
    the array has to be consistent with the grid set class so it can access the x/ydist parameters
    """
    # check if grid_set has grid info
    if not grid_set.gridinfo:
        print("No grid_set geo grid info - no result")
        return False
    in_mn = np.shape(u)
    if in_mn[0]!=grid_set.m or in_mn[1]!=grid_set.n :
        print("input array or geo grid_set not consistently shaped")
        return False
    else:
        
        Dvdx = geo_gradient(v,grid_set)[1]
        Dudy = geo_gradient(u,grid_set)[0]

        zeta = Dvdx - Dudy

        return zeta

    
    
def de_ripple(array1,array2,rip_filt_std = 1,filt_ring_sig = 5,force_zero = False):
    from scipy.ndimage.filters import gaussian_filter
    # find the ripples by subtracting the arrays
    ripples = array1 - array2
    # fast fourier transform the difference
    rip_spec  = np.fft.fft2(np.double(ripples))
    rip_spec2 = np.fft.fftshift(rip_spec)
    # find the ring sprectrum the contains the ripples
    filt_ring = np.ones_like(array1)
    spec_r = np.mean(rip_spec2) + rip_filt_std*np.std(rip_spec2)
    filt_ring[rip_spec2>spec_r] = 0.0
    filt_ring = gaussian_filter(filt_ring,sigma = filt_ring_sig)
    if not type(force_zero) == bool:
        filt_ring[rip_spec2>spec_r] = filt_ring[rip_spec2>spec_r]*force_zero  
    # use this filter ring to remove the array1 fft spectrum
    a1_spec  = np.fft.fft2(np.double(array1))
    a1_spec2 = np.fft.fftshift(a1_spec)
    a1_spec2 = a1_spec2*filt_ring
    back = np.real(np.fft.ifft2(np.fft.ifftshift(a1_spec2)))
    return back

def geo_filter(array,grid_set,distance,mask = False):
    from scipy.ndimage.filters import gaussian_filter
    """
    filter function that will take the grid info from the 
    grid_set type class to get filter distances
    the array has to be consistent with the grid set class so it can access the x/ydist parameters
    """
    # takes the DOT and filters out the geoid harmonics
    # hopefully can implement variable gradient using 
    # grid info
    # can dx/dyres if needed
    # check if grid_set has grid info
    if type(mask)==bool:
        if mask:
            mask = grid_set.mask
        else:
            mask = np.ones_like(array)
    elif (np.shape(mask)[0] != grid_set.m 
           |np.shape(mask)[1] != grid_set.n):# check mask dimension)
        print("Mask array incorrect shape, ignoring it")
        mask = np.ones([grid_set.m,grid_set.n])
    if not grid_set.gridinfo:
        print("No grid_set geo grid info - no result")
        return False
    in_mn = np.shape(array)
    if in_mn[0]!=grid_set.m or in_mn[1]!=grid_set.n :
        print("input array or geo grid_set not consistently shaped")
        return False
    else:
        V = np.empty_like(array) 
        W = np.empty_like(array) 
        out_array = np.empty_like(array) 
        f_sig =[distance/d for d in [grid_set.dxRes,grid_set.dyRes]] # some function of the radius given..
        V[:,:]=array*mask
        V[np.isnan(V)]=0
        VV=gaussian_filter(V,sigma=f_sig)

        W[:,:]=0*array+1
        W = W*mask
        W[np.isnan(W)]=0
        WW=gaussian_filter(W,sigma=f_sig)

        out_array[:,:]=VV/WW
        out_array[np.isnan(array)] = np.nan
        
        return out_array*mask
    
def geo_convolve(array,grid_set,distance,limits,mask = False,set_kernel = False):
    from astropy.convolution import convolve, Gaussian2DKernel
    """
    filter function that will take the grid info from the 
    grid_set type class to get filter distances
    the array has to be consistent with the grid set class so it can access the x/ydist parameters
    mask is the same size as the array
    """
    # takes the DOT and filters out the geoid harmonics
    # hopefully can implement variable gradient using 
    # grid info
    # can dx/dyres if needed
    # check if grid_set has grid info
    if type(mask)==bool:
        if mask:
            mask = grid_set.mask
        else:
            mask = np.ones_like(array)
    elif (np.shape(mask)[0] != grid_set.m 
           or np.shape(mask)[1] != grid_set.n):# check mask dimension)
        print(np.shape(mask)[0],grid_set.m )
        print(np.shape(mask)[1],grid_set.n )
        print("Mask array incorrect shape, ignoring it")
        mask = np.ones([grid_set.m,grid_set.n])
    if not grid_set.gridinfo:
        print("No grid_set geo grid info - no result")
        return False
    if type(set_kernel)==bool:
        # set normal guassian kernel as a function of grid general dim
        f_sig =np.mean([distance/d for d in [grid_set.dxRes,grid_set.dyRes]]) 
        kernel = Gaussian2DKernel(f_sig)
    else: kernel = set_kernel
    in_mn = np.shape(array)
    if in_mn[0]!=grid_set.m or in_mn[1]!=grid_set.n :
        print("input array or geo grid_set not consistently shaped")
        return False
    else:
        # some function of the radius given..
        array_2 = copy.copy(array)
        array_2[array<limits[0]] = np.nan
        array_2[array>limits[1]] = np.nan 
        array_2[np.isnan(mask)] = np.nan 
        
        out_array = convolve(array_2,kernel,boundary = 'extend')
        out_array[np.isnan(mask)] = np.nan 
        
        return out_array

# takes generic data and regrids it into a data_year
def regrid_data(data,dates,lons,lats,grid_set,periods,
                fill_lonlat = False):
    """
    makes a data year object, nicely regridded on D_Class grid
    time dimension of data is default 0 
    currently setup to access list of lists, or arrays
    first list access is the time point
    retains old netcdf option to fill lat lon arrays from singular 
    axis arrays
    otherwise lon/lat need to be of the same shape as the data time slice
    periods is the number of time slices per year, ie. 12 for monthlies
    """
    n_t = np.shape(data)[0]
    
    new_d_array = np.empty([n_t,grid_set.m,grid_set.n])
    # if the lat_lons need filling - do it
    if fill_lonlat:
        lon_a,lat_a = np.meshgrid(lons,lats)
    else:
        lon_a = lons
        lat_a = lats
    # regrid depending upon m and grid
    x_d, y_d = grid_set.mplot(lon_a, lat_a)
    for tt in range(n_t):
        new_d_array[tt,:,:] = griddata((x_d.ravel(), y_d.ravel()),
                data[tt][:].ravel(), (grid_set.xpts.T, grid_set.ypts.T),
                method='linear')
    return dy.data_year(new_d_array,dates,periods)

# takes generic data and regrids it into a data_year
def regrid_vectors(x,y,dates,lons,lats,grid_set,periods,
                fill_lonlat = False,vector_angles = False):
    """
    makes a vector data year object, nicely regridded on D_Class grid
    time dimension of data is default 0 
    currently setup to access list of lists, or arrays
    first list access is the time point
    retains old netcdf option to fill lat lon arrays from singular 
    axis arrays
    otherwise lon/lat need to be of the same shape as the data time slice
    periods is the number of time slices per year, ie. 12 for monthlies
    # the original vectors may need to be rotated back to be square to
    # lon lat so they can be regridded
    # if vector_angles = false then they are already square ie on an x/y frame sqaure to lon/lat
    # otherwise vector_angles is the same shape as lon/lats etc 
    #  and is angle positive from gridded data x/y to lon/lat  
    # ie positive rotational angle from local y positive dimension to true north
    # so angle array is consistent to gridinfo method on a grid_set - so you can use that.
    """
    n_t = np.shape(x)[0]
    
    new_x_array = np.empty([n_t,grid_set.m,grid_set.n])
    new_y_array = np.empty([n_t,grid_set.m,grid_set.n])
    # if the lat_lons need filling - do it
    if fill_lonlat:
        lon_a,lat_a = np.meshgrid(lons,lats)
    else:
        lon_a = lons
        lat_a = lats
    if type(vector_angles) == bool:
        orig_c = np.ones_like(lon_a)
        orig_s = np.zeros_like(lon_a)
    else: 
        orig_c = np.cos(np.deg2rad(vector_angles.T))
        orig_s = np.sin(np.deg2rad(vector_angles.T))
    # regrid depending upon mplot and grid
    x_d, y_d = grid_set.mplot(lon_a, lat_a)
    for tt in range(n_t):
        # rotating back to lon lat
        orig_x = x*orig_c - y*orig_s
        orig_y = y*orig_c + x*orig_s
        # regridding
        temp_x = griddata((x_d.ravel(), y_d.ravel()),
                orig_x.ravel(), (grid_set.xpts.T, grid_set.ypts.T),
                method='linear')
        temp_y = griddata((x_d.ravel(), y_d.ravel()),
                orig_y.ravel(), (grid_set.xpts.T, grid_set.ypts.T),
                method='linear')
        
        # rotating to the new grid
        new_x_array[tt] = temp_x*grid_set.ang_c.T - temp_y*grid_set.ang_s.T 
        new_y_array[tt] = temp_y*grid_set.ang_c.T + temp_x*grid_set.ang_s.T 
        
    return dy.vec_data_year(new_x_array,new_y_array,dates,periods)


@jit
def ellipsoidal_distance(long1, lat1, long2, lat2,deg=False,eps=1e-10):
    """
    (long1, lat1, long2, lat2) all in radians
    outputs a distance in m
    """
    if np.isnan([long1, lat1, long2, lat2]).any():
        return np.nan
    else:
        if deg:
            long1 = np.deg2rad(long1)
            lat1  = np.deg2rad(lat1)
            long2 = np.deg2rad(long2)
            lat2  = np.deg2rad(lat2)

        a = 6378137.0 # equatorial radius in meters
        f = 1/298.257223563 # ellipsoid flattening
        b = (1 - f)*a
        tolerance = eps # to stop iteration

        phi1, phi2 = lat1, lat2
        U1 = np.arctan((1-f)*np.tan(phi1))
        U2 = np.arctan((1-f)*np.tan(phi2))
        L1, L2 = long1, long2
        L = L2 - L1
        i = 0

        lambda_old = L + 0

        while True:

            t =  (np.cos(U2)*np.sin(lambda_old))**2
            t += (np.cos(U1)*np.sin(U2) - np.sin(U1)*np.cos(U2)*np.cos(lambda_old))**2
            sin_sigma = t**0.5
            cos_sigma = np.sin(U1)*np.sin(U2) + np.cos(U1)*np.cos(U2)*np.cos(lambda_old)
            sigma     = np.arctan2(sin_sigma, cos_sigma)

            sin_alpha    = np.cos(U1)*np.cos(U2)*np.sin(lambda_old) / (sin_sigma)
            cos_sq_alpha = 1 - sin_alpha**2
            cos_2sigma_m = cos_sigma - 2*np.sin(U1)*np.sin(U2)/(cos_sq_alpha+1e-12)
            C            = f*cos_sq_alpha*(4 + f*(4-3*cos_sq_alpha))/16

            t          = sigma + C*sin_sigma*(cos_2sigma_m + C*cos_sigma*(-1 + 2*cos_2sigma_m**2))
            lambda_new = L + (1 - C)*f*sin_alpha*t
            if np.abs(lambda_new - lambda_old) <= tolerance:
                break
            elif i > 1000:
                return np.nan
                break
            else:
                lambda_old = lambda_new
                i += 1

        u2 = cos_sq_alpha*((a**2 - b**2)/b**2)
        A  = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320 - 175*u2)))
        B  = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
        t  = cos_2sigma_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2sigma_m**2))
        t -= (B/6)*cos_2sigma_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2sigma_m**2)
        delta_sigma = B * sin_sigma * t
        s = b*A*(sigma - delta_sigma)

        return s


@jit
def lon_lat_angle( lon1,lat1,lon2,lat2,deg=False,return_trig = False ):
    """
    #LAT_LON_ANGLE finds the geodesic angle from point 1 to point 2 
    #(lat lon in radians)
    #   This done by a series of axes rotations to get the 1st point at 0,0
    #   keeping the second the same relative distance apart. The roataion
    #   needed to get the 2nd point directly north of the 1st is the geodesic
    #   angle.
    """
    if np.isnan([lon1,lat1,lon2,lat2]).any():
        if return_trig: return np.nan, np.nan
        else: return np.nan
    else:
        if deg:
            lon1 = np.deg2rad(lon1)
            lat1 = np.deg2rad(lat1)
            lon2 = np.deg2rad(lon2)
            lat2 = np.deg2rad(lat2)

        C_lat=np.cos(lat2);
        S_lat=np.sin(lat2);
        C_lon=np.cos(lon2);
        S_lon=np.sin(lon2);

        C_1=np.cos(-lon1);
        S_1=np.sin(-lon1);
        C_2=np.cos(-lat1);
        S_2=np.sin(-lat1);

        A1=[[C_1, -S_1, 0],
            [S_1,  C_1, 0],
            [0,    0,   1]]

        A2=[[C_2, 0, -S_2],
            [0,   1,  0  ],
            [S_2, 0,  C_2]]

        Borig=[C_lat*C_lon,
               C_lat*S_lon,
               S_lat      ];

        B=np.matmul(A2,A1)
        B=np.matmul(B,Borig)
        # print(B)

        if return_trig:
            scale=np.hypot(B[1],B[2])

            angle_sin=-B[1]/scale
            angle_cos= B[2]/scale

            return angle_cos, angle_sin
        else:
            angle=np.arctan2(-B[1],B[2])

            return angle

def nearest_xy(lon,lat,grid_set): 
    x,y = np.unravel_index(np.argmin(
        np.abs(grid_set.lons - lon) + 
        np.abs(grid_set.lats - lat)),
        np.shape(grid_set.lons))
    return x,y
    

def nearest_dy(lon,lat,t,gs,dy,tr = [0,0],box = [0,0],time_vec = False,space_array = False):
    """
    give this a dy object and a gs object,
    the nearest point to the supplied lon lat will be returned
    tr is a time range option [time points previous, after]
    if tr > 0 time_vec=True will return a rs/H/WAVES/SWH/swh_arrays/SWH_10vector of the time point, False is nanmean
    box is a spatial range option [range in x, range in y]
    if there is a box, space_array=True returns the whole box, False is nanmean
    """
    y,x = nearest_xy(lon,lat,gs)
    
    out_array = dy[t-tr[0]:t+tr[1]+1,x-box[0]:x+box[0]+1,y-box[1]:y+box[1]+1]
    if time_vec and space_array:
        return out_array
    elif time_vec:
        return np.nanmean(out_array,axis = (1,2))
    elif space_array:
        return np.nanmean(out_array,axis = 0)
    else:
        return np.nanmean(out_array)


def nearest_interp(lon,lat,t,gs,dy,tr = [0,0],time_vec = False):
    """
    give this a dy object and a gs object,
    the nearest point to the supplied lon lat will be returned
    tr is a time range option [time points previous, after]
    if tr > 0 time_vec=True will return a vector of the time point, False is nanmean
    box is a spatial range option [range in x, range in y]
    if there is a box, space_array=True returns the whole box, False is nanmean
    """
    x,y = nearest_xy(lon,lat,gs)
    # find 4 point weighting via x,y or lon,lat?
    xpts,ypts = gs.mplot(lon,lat)
    
    if xpts<gs.xpts[x,y]:
        try:
            xw = 1.0 - (gs.xpts[x,y]-xpts)/(gs.xpts[x,y] - gs.xpts[x-1,y])
        except IndexError:
            xw = 1.0
    else:
        try:
            xw = 1.0 - (xpts - gs.xpts[x,y])/(gs.xpts[x+1,y] - gs.xpts[x,y])
        except IndexError:
            xw = 1.0
    if ypts<gs.ypts[x,y]:
        try:
            yw = 1.0 - (gs.ypts[x,y]-ypts)/(gs.ypts[x,y] - gs.ypts[x-1,y])
        except IndexError:
            yw = 1.0
    else:
        try:
            yw = 1.0 - (ypts - gs.ypts[x,y])/(gs.ypts[x+1,y] - gs.ypts[x,y])
        except IndexError:
            yw = 1.0
        
#         print(sigi,sigw,tori,torw)
        

    try:
        pwf1 =      xw *     yw *gs[t-tr[0]:t+tr[1]+1,x  ,y]
    except IndexError:
        pwf1 = np.nan
    try:
        pwf2 = (1.0-xw)*     yw *gs[t-tr[0]:t+tr[1]+1,x+1,y]
    except IndexError:
        pwf2 = np.nan
    try:
        pwf3 =      xw *(1.0-yw)*gs[t-tr[0]:t+tr[1]+1,x  ,y+1]
    except IndexError:
        pwf3 = np.nan
    try:
        pwf4 = (1.0-xw)*(1.0-yw)*gs[t-tr[0]:t+tr[1]+1,x+1,y+1]
    except IndexError:
        pwf4 = np.nan
    
    if time_vec:
        return np.nansum([pwf1,pwf2,pwf3,pwf4],axis = 1)
    else:
        return np.nansum([pwf1,pwf2,pwf3,pwf4])
    
#     out_array = dy[t-tr[0]:t+tr[1]+1,x-box[0]:x+box[0]+1,y-box[1]:y+box[1]+1]
#     if time_vec and space_array:
#         return out_array
#     elif time_vec:
#         return np.nanmean(out_array,axis = (1,2))
#     elif space_array:
#         return np.nanmean(out_array,axis = 0)
#     else:
#         return np.nanmean(out_array)


# def gen_from_lonlat():
#     """
#     Creates a lon/lat square grid from a 1d vector of each
#     Returns a grid set with points defined 
#     """

class Gs2Gs:
    """
    For regridding
    give a new projection, this builds the internals of griddata, but for repeat usage
    feed it two grid_sets and it'll be a function
    regridding from one grid_set to the other
    """
    def __init__(self,gs_native,gs_new,vectors = False,vectors_plot = False):
        """
        gs_native is the grid set the data is defined on
        gs_new is where you want it to go
        set vectors = True if you want to regrid vectors too
        this will require all the correct angular grid info
        """
        from scipy.spatial import Delaunay
        # get the tri angulation
        self.vectors = vectors
        self.vectors_plot = vectors_plot
        if self.vectors:
            self.vectors_plot = True
        xorig,yorig = gs_new.mplot(gs_native.lons,gs_native.lats)
        xyorig = np.vstack((xorig.ravel(),yorig.ravel())).T
        self.tri = Delaunay(xyorig)  # Compute the triangulation
        # destination mesh
        self.mesh_new = (gs_new.xpts,gs_new.ypts)
        if vectors:
            # record the neccesary angles to de-rotate the input vectors 
            # and re-rotate the output vectors
            self.in_ang_c = gs_native.ang_c
            self.in_ang_s = gs_native.ang_s
            self.out_ang_c = gs_new.ang_c
            self.out_ang_s = gs_new.ang_s
            self.new_mplot = gs_new.mplot
            self.new_lons = gs_new.lons
            self.new_lats = gs_new.lats
        if vectors_plot:
            # record the neccesary angles to de-rotate the input vectors 
            # and re-rotate the output vectors
            self.in_ang_c = gs_native.ang_c
            self.in_ang_s = gs_native.ang_s
            self.new_mplot = gs_new.mplot
            self.new_lons = gs_new.lons
            self.new_lats = gs_new.lats

    def rg_array(self,arr):
        """
        the regridding function
        feed it the array defined on gs_native
        out pops a new array on gs_new
        """
        from scipy.interpolate import LinearNDInterpolator
        # define the function
        interpolator = LinearNDInterpolator(self.tri, arr.ravel())
        return interpolator(self.mesh_new)


    def rg_vecs_to_plot(self,x,y):
        """
        the regridding function, this just allows plotting in the new GS. Plotting vectors are not neccesarily equal to gridding - non-square to projection girds for example
        feed it the x,y comps defined on gs_native
        out pops new arrays on gs_new
        """
        from scipy.interpolate import LinearNDInterpolator
        if self.vectors_plot:
            # de-rotate the input vecs (back to lon lat square)
            xr = -y*self.in_ang_c - x*self.in_ang_s
            yr =  x*self.in_ang_c - y*self.in_ang_s
#             xr = x*self.in_ang_c - y*self.in_ang_s
#             yr = y*self.in_ang_c + x*self.in_ang_s
            # define the function
            interpolator = LinearNDInterpolator(self.tri, xr.ravel())
            # use it
            xrr = interpolator(self.mesh_new)
            # define the function
            interpolator = LinearNDInterpolator(self.tri, yr.ravel())
            # use it
            yrr = interpolator(self.mesh_new)
            return self.new_mplot.rotate_vector(xrr,yrr,
                      self.new_lons,self.new_lats)
        else:
            print('Gs2Gs not defined for vectors, re-initialise')
            
    
    def rg_vecs(self,x,y):
        """
        the regridding function
        feed it the x,y comps defined on gs_native
        out pops a new array on gs_new
        """
        from scipy.interpolate import LinearNDInterpolator
        if self.vectors:
            # de-rotate the input vecs (back to lon lat square)
            xr = -y*self.in_ang_c - x*self.in_ang_s
            yr =  x*self.in_ang_c - y*self.in_ang_s
#             xr = x*self.in_ang_c - y*self.in_ang_s
#             yr = y*self.in_ang_c + x*self.in_ang_s
            # define the function
            interpolator = LinearNDInterpolator(self.tri, xr.ravel())
            # use it
            xrr = interpolator(self.mesh_new)
            # define the function
            interpolator = LinearNDInterpolator(self.tri, yr.ravel())
            # use it
            yrr = interpolator(self.mesh_new)
            
            xrout = -yrr*self.out_ang_c + xrr*self.out_ang_s
            yrout =  xrr*self.out_ang_c + yrr*self.out_ang_s
            return xrout, yrout
        else:
            print('Gs2Gs not defined for vectors, re-initialise')
            
def border(arr): # Input array : arr
    alist=[arr[0,:-1], arr[:-1,-1], arr[-1,::-1], arr[-2:0:-1,0],[arr[0,0]]]
    
    return np.concatenate(alist)


class vary_smooth2d():
    """
    variable smooth object
    initalises with varc,varr which need to give local 
    number of grid cells to smooth over
    Initialising creats a large sparse array to use for later repeated smoohting operations
    Inputs
    varr: array of local row smoothing distances (no. of grid cells, float)
    varc: array of local column smoothing distancese (no. of grid cells, float)
    """
    def __init__(self,varr,varc,verbos = False):
        ii,jj = varc.shape
        if verbos:
            print('Bulding vary smoother, av cell dist = ['+
                  '{:0.1f}'.format(np.nanmean(varr))+', '+ 
                  '{:0.1f}'.format(np.nanmax(varr))+'], ['+ 
                  '{:0.1f}'.format(np.nanmean(varc))+', '+ 
                  '{:0.1f}'.format(np.nanmax(varc))+']')
        varc[np.isnan(varc)] = 1.0
        varc[np.isinf(varc)] = 1.0
        varr[np.isnan(varr)] = 1.0
        varr[np.isinf(varr)] = 1.0
        varc[varc<1.0] = 1.0
        varr[varr<1.0] = 1.0
        varc = varc.astype(int)
        varr = varr.astype(int)
        self.ii=ii
        self.jj=jj
        self.matrix = sparse.lil_matrix((ii*jj,ii*jj))
        for i in range(ii):
            for j in range(jj):
                ## i,j are the target indicies
                ij = np.ravel_multi_index((i,j),(ii,jj))
                ## i+-varr[i,j]
                ## j+-varc[i,j]
                ## are the target indice
                ## place in the sparse matrix is np.ravel_multi_index((i,j))
                vij = []
                for vi in range(i-varr[i,j]+1,i+varr[i,j]):
                    for vj in range(j-varc[i,j]+1,j+varc[i,j]):
                        if (vi in range(ii)) and (vj in range(jj)):
                            vij.append(np.ravel_multi_index((vi,vj),(ii,jj)))
                weight = 1/np.shape(vij)[0]
        #         if np.shape(vij)[0] > 1: print(weight)
                for vv in vij:
                    self.matrix[ij,vv] = weight
        #             matrix[vv,ij] = weight
        #             if vv!=ij:print(vv,ij)
    def smooth(self,array):
        """
        Takes the initalised object with local smoothing distances
        Smooths an array using the pre-given distances
        input array here is the smae size as the initialise varr,varc
        """
        A = np.isnan(array)
        array[A] = 0.0
        a_r = array.ravel()
        narray = (self.matrix*a_r).reshape((self.ii,self.jj))
        narray[A] = np.nan
        return narray
    
class geo_vary_smooth():
    """
        Uses the dimensions of a grid_set class to set appropriate
        Dims for a vary_smooth class
        init distance is the approximate distance we want represented 
        by the local box dimensions
    """
    def __init__(self,grid_set,distance,verbos=False):
        dist = 75e3
        varr = distance/grid_set.xdist
        varc = distance/grid_set.ydist
        self.Vsm2d = vary_smooth2d(varr,varc,verbos=verbos)
    
    def smooth(self,array):
        return self.Vsm2d.smooth(array)
        #som