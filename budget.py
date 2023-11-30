# functions for budget
import numpy as np
import pickle
import datetime as dt
from dateutil.relativedelta import relativedelta
from scipy.sparse import spdiags
from copy import copy
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def advection(U,V,Vol,G):
    advtemp=-V[1:-1,1:-1]*(
        (Vol[1:-1,2:]-Vol[1:-1,:-2])/(2*G.ydist[1:-1,1:-1])) \
           -U[1:-1,1:-1]*(
        (Vol[2:,1:-1]-Vol[:-2,1:-1])/(2*G.xdist[1:-1,1:-1]))
    return np.pad(advtemp,((1,1),(1,1)),'constant',constant_values = np.nan)

def divergence(U,V,Vol,G):
    divtemp=-Vol[1:-1,1:-1]*(
        (V[1:-1,2:]-V[1:-1,:-2])/(2*G.ydist[1:-1,1:-1])) \
           -Vol[1:-1,1:-1]*(
        (U[2:,1:-1]-U[:-2,1:-1])/(2*G.xdist[1:-1,1:-1]))
    return np.pad(divtemp,((1,1),(1,1)),'constant',constant_values = np.nan)


def intensification(Vols,nsecs):
    return (Vols[2,:,:]-Vols[0,:,:])/(2*nsecs)
# import grid_set as gs



# class input_slices():
#     ### holds all the input data and does all the loading
#     def __init__(self,G,get_func,In,smooth=False,Rgd=False,tsmth_rad=0,vectors=False)
#     """
#     setup with options for
#     size: use G for it
#     time smoothing: give int
#     give appropriate: get_''_array
#     give appropriate: input class
#     give appropriate: smoothing function
#     give appropriate: regridding object
    
#     defs area: 
#     intialise
#     pre-fill (before loop)
#     loop_fill (during loop)
#     """
#     self.G = G
#     self.get_func = get_func
#     self.In = In
#     self.vectors = vectors
#     ### allowing for smooth to pass options
#     if type(smooth) == function:
#         self.smooth = True
#         self.smoothF = smooth
#     else:
#         self.smooth = False
#     if type(Rdg) == gs.Gs2Gs:
#         self.Rdg = True
#         self.RdgF = Rdg
#     else:
#         self.Rdg = False
#     if tsmth_rad > 0:
#         ## then we're in the moving average pain cave
#         ## have to step into the future to retrieve the next value
#         ## how is that going to work with data holes?
#         self.time_smooth = True
#         if self.vectors:
#             self.TUslices = np.zeros([1+tsmth_rad,self.G.m,self.G.n])
#             self.TVslices = np.zeros([1+tsmth_rad,self.G.m,self.G.n])
#         else:
#             self.Tslices = np.zeros([1+tsmth_rad,self.G.m,self.G.n])
#     if self.vectors:
#         self.Uslices = np.zeros([3,self.G.m,self.G.n])
#         self.Vslices = np.zeros([3,self.G.m,self.G.n])
#     else:
#         self.slices = np.zeros([3,self.G.m,self.G.n])
    
    
#     def prefill(self,ts,tw):
#         """
#         we're now going to fill the slices for the beginnig
#         for no time smoothing - fill slice 0 from ts,tw
#         the copy to slice 1
#         """
#         if self.time_smooth:
#             ## update fill Tslice array 
            
#             ## use Tslice array to fill slice array
#         else:
#             ##  just fill slice array
#             if self.vectors
#                 Utemp,Vtemp = self.get_func(InputV,ts,tw)
#             if self.smooth:
#                 Uslices[0] = b.smooth2a(Utemp,velsmrad)
#     Vslices[0] = b.smooth2a(Vtemp,velsmrad)
        
        
        
#     def loopfill(self,ts,tw):

def smooth2a(inputM,N):
    # %Copied from Paul Holland 
    # %Smooths 2D array data.  Ignores NaN's.
    # %
    # % This function smooths the data in matrixIn using a mean filter over a
    # % rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
    # % element "i" by the mean of the rectange centered on "i".  Any NaN
    # % elements are ignored in the averaging.  If element "i" is a NaN, then it
    # % will be preserved as NaN in the output.  At the edges of the matrix,
    # % where you cannot build a full rectangle, as much of the rectangle that
    # % fits on your matrix is used (similar to the default on Matlab's builtin
    # % function "smooth").
    # % 
    # % "matrixIn": original matrix
    # % "Nr": number of points used to smooth rows
    # % "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
    # % 
    # % "matrixOut": smoothed version of original matrix
    # % 
    # % 
    # % 	Written by Greg Reeves, March 2009.
    # % 	Division of Biology
    # % 	Caltech
    # % 
    # % 	Inspired by "smooth2", written by Kelly Hilands, October 2004
    # % 	Applied Research Laboratory
    # % 	Penn State University
    # % 
    # % 	Developed from code written by Olof Liungman, 1997
    # % 	Dept. of Oceanography, Earth Sciences Centre
    # % 	Gï¿½teborg University, Sweden
    # % 	E-mail: olof.liungman@oce.gu.se

    matrixIn = copy(inputM)

    row, col = matrixIn.shape
    #print(row,col)
#     diags=np.array([-2, -1, 0, 1, 2])
    diags=np.linspace(-N, N, num=2*N+1)
    eL = spdiags(np.ones((2*N+1,row)),diags,row,row).toarray()
    eR = spdiags(np.ones((2*N+1,col)),diags,col,col).toarray()


    # %
    # % Setting all "NaN" elements of "matrixIn" to zero so that these will not
    # % affect the summation.  (If this isn't done, any sum that includes a NaN
    # % will also become NaN.)
    # %

    A = np.isnan(matrixIn)
    matrixIn[A]=0.0
    
    # %
    # % For each element, we have to count how many non-NaN elements went into
    # % the sums.  This is so we can divide by that number to get a mean.  We use
    # % the same matrices to do this (ie, "eL" and "eR").
    # %

    nrmlize = np.dot(np.dot(eL,np.logical_not(A)),eR)
    nrmlize[A]=np.nan

    # % Actually taking the mean.

    matrixOut = np.dot(np.dot(eL,matrixIn[:,:]),eR)
    matrixOut = matrixOut/nrmlize
    
    return matrixOut



# we now need to work on the dates
# got through entire list and return those within the givne point and range
# first checks finest resolution is a day
# frist and last should be outliers
def select_dates(dates,time_u,time_w,diag = False):
    window = False
    dates_u = []
    # extra check for aligning frst date
    # then we need an empty extended entry
    if (time_u - dates[0] ).days   == 0:
        dates_u.append([])
    for n,d in enumerate(dates):
        # go along list checking if we get passed time_u

        if ((time_u - d ).days -1 < 0) and not window:
            window = True
#             print(d.strftime('%Y%m%d'))
            # we also need the 'bracket slices'
            # either side of the date in interest
            # so to interpolate when there isn't any
            if n > 0: dates_u.append(dprev)  

        if window: dates_u.append(d)    

        if (time_u + time_w - d ).days -1 < 0:
#             print(d.strftime('%Y%m%d'))
            window = False
            break
        dprev = d
    # if we've reached the end and swe're still windowing, add an empty entry
    if window:
        dates_u.append([])
    if diag:
        [print(d.strftime('%Y%m%d')) for d in dates_u if type(d)==dt.datetime]
    if len(dates_u)==1:
        dnew = [[],dates_u[0],[]]
        dates_u = dnew
    return dates_u



def get_load_points(datesQ,time_u,time_w,diag=False):
    # need to return d_load
    w0 = 0.0
    w1 = 0.0
    load_I = False
    if diag:
        print('datesQ length = '+str(np.shape(datesQ)[0]))
        if np.shape(datesQ)[0] ==1:
            print(datesQ[0].strftime('%Y%m%d'))
    if np.shape(datesQ)[0]>2:
        # we don't need the end point so load from
        # between
        d_load = datesQ[1:-1]

    #### this next block wiil need adding if we want weird
    #### combos of days (say 5) that can cross month boundaries
#     elif np.shape(datesQ)[0]>3:
#         # we may still need to interp if the gap is 
#         # not perfectly aligned
#         d_load = datesQ[1:-1]

        
    elif time_w.months>0:
        # if we're doing monthlies and the first end point hits
        # the first day then take that
        # aslo if monthlies we do soemthing clever with monthly points
        if (time_u - datesQ[0]).days == 0:
            d_load = [datesQ[0]]
        if (time_u - datesQ[1]).days == 0:
            d_load = [datesQ[1]]

    elif (time_w.days==1) and ((time_u - datesQ[0]).days == 0):
        # if we're doing single days and our day hits the first point exactly do that
        d_load = [datesQ[0]]
    elif (time_w.days==1) and ((time_u - datesQ[1]).days == 0):
        d_load = [datesQ[1]]

    else: 
        # then we're doing some sort of days interp so use the endpoints
        load_I = True
        # find weights depending on range we're interested in
        dw = time_w.days
        dgap = (datesQ[1]- datesQ[0]).days
        dgap0 = (time_u- datesQ[0]).days
        w0 = np.sum([1-n/dgap 
                     for n in range(dgap0,dgap0+dw)])/dw
        w1 = np.sum([n/dgap 
                     for n in range(dgap0,dgap0+dw)])/dw
        d_load = datesQ
    return d_load,load_I,w0,w1



def get_hi_array(bfhi,time_u,time_w,diag=False,verbos=False):
    # bfhi has to be preloaded with dates
    duse = select_dates(bfhi.dates,time_u,time_w,diag=diag)
    if diag: print('duse length = '+str(np.shape(duse)[0]))
    dload,interp,w0,w1 = get_load_points(duse,time_u,time_w,diag=diag)
    if diag: print('dload length = '+str(np.shape(dload)[0]))
    hi_in = bfhi.get_hi(dload,verbos=verbos)
    if interp:
        hi_out = hi_in[0]*w0 + hi_in[1]*w1
        if verbos: print("interpolating "+bfhi.name+" thickness")
    else:
        hi_out = np.nanmean(hi_in,axis=0)
        sno = np.shape(dload)[0]
        if verbos: print("loading "+bfhi.name+" thickness, "+str(sno)+" slices")
        if verbos: print("mean hi: ",'{:.3}'.format(np.nanmean(hi_out)))
    return hi_out


def get_aice_array(bfaice,time_u,time_w,diag=False,verbos=False):
    # bfaice has to be preloaded with dates
    duse = select_dates(bfaice.dates,time_u,time_w,diag=diag)
    if diag: print('duse length = '+str(np.shape(duse)[0]))
    dload,interp,w0,w1 = get_load_points(duse,time_u,time_w,diag=diag)
    if diag: print('dload length = '+str(np.shape(dload)[0]))
    hi_in = bfaice.get_aice(dload,verbos=verbos)
    if interp:
        hi_out = hi_in[0]*w0 + hi_in[1]*w1
        if verbos: print("interpolating "+bfaice.name+" concentration")
    else:
        hi_out = np.nanmean(hi_in,axis=0)
        sno = np.shape(dload)[0]
        if verbos: print("loading "+bfaice.name+" concentration, "+str(sno)+" slices")
        if verbos: print("mean aice: ",'{:.3}'.format(np.nanmean(hi_out)))
    return hi_out



def get_vels_array(bfvel,time_u,time_w,diag=False,verbos=False):
    # bfvel has to be preloaded with dates
    duse = select_dates(bfvel.dates,time_u,time_w,diag=diag)
    if diag: print('duse length = '+str(np.shape(duse)[0]))
    dload,interp,w0,w1 = get_load_points(duse,time_u,time_w,diag=diag)
    if diag: print('dload length = '+str(np.shape(dload)[0]))
    u_in,v_in = bfvel.get_vels(dload,verbos=verbos)
    if interp:
        u_out = u_in[0]*w0 + u_in[1]*w1
        v_out = v_in[0]*w0 + v_in[1]*w1
        if verbos: print("interpolating "+bfvel.name+" velocity")
    else:
        u_out = np.nanmean(u_in,axis=0)
        v_out = np.nanmean(v_in,axis=0)
        sno = np.shape(dload)[0]
        if verbos: print("loading "+bfvel.name+" velocity, "+str(sno)+" slices")
        if verbos: print("mean vel: ",'{:.3}'.format(np.nanmean(np.hypot(u_out,v_out))))
    return u_out,v_out



def get_loop_list(dsearch,tw,include_last = False):
    dates_u = []
    #### NOTE this only works well with single days or months
    #### when operating on monthly date points
    # search date list using tw to find all possible
    # dates that can successfully access data

    # loop over dates list
    dlen = np.shape(dsearch)[0]
    for n,d in enumerate(dsearch):
        # dsearch[0] entry will be in it
        if n==0: dates_u.append(d)
        # we need to skip all entries that are within a window
        elif (dates_u[-1]+tw-d).days > 0: continue
        else: dates_u.append(d)
        # if this is the last entry then stop
        if n==dlen-1: break
        # now progress beyond d in steps of tw 
        search  = True
        # need to check if there are gaps
        # limit will be a gap of greater than a month
        # in which case don't search 
        if (dsearch[n+1]-d).days >32: search = False
        # dnow is the point of interest
        dnow = d + tw
        while search:
            # first check if dnow is equal or passed the next entry
            if (dnow-dsearch[n+1]).days > -1:
                search = False
            else:
                dates_u.append(dnow)
                dnow += tw
    # print the important info
    print('made loop list length '+str(np.shape(dates_u)[0]),flush=True)
    if np.shape(dates_u)[0]>0:
        print('Start : '+dates_u[ 0].strftime('%Y%m%d'),flush=True)
        print('Finish: '+dates_u[-1].strftime('%Y%m%d'),flush=True)
    if include_last:
        return dates_u
    else:
        return dates_u[:-1]


def comp_dates(dates1,dates2):
    dout = []
    for d in dates1:
        if d in dates2:
            dout.append(d)
    print('made loop list length '+str(np.shape(dout)[0]),flush=True)
    if np.shape(dout)[0]>0:
        print('Start : '+dout[ 0].strftime('%Y%m%d'),flush=True)
        print('Finish: '+dout[-1].strftime('%Y%m%d'),flush=True)
    return dout
    

def get_step_size(d,tw):
    if tw.months>0:
        # get no. of days in this month
        ndays = ((d+tw) - d).days
        step_size = ndays*86400
    elif tw.days>0:
        step_size = tw.days*86400
    return step_size


class budhist:
    """
    Contains all the history variables and accumulation
    Also has dumping options
    Also has budget history variables 
    Also has drift statistics 
    Adding opitons - we need dump window and dumpd0
    Output regridding - this will be a Gs2Gs option - output on a particular grid
        this is an option in dump_nc 
    """
    def __init__(self,G,dlist,tw,drift_stats = False,
                 dumptw = relativedelta(months=1),budday0 = 'default',
                 split_div=False,transport = False,tot_budget = False,save_input=True):
        ## set first date
        self.split_div = split_div
        self.tot_budget = tot_budget
        self.transport = transport
        self.save_input = save_input
        self.ddump = dlist[0]
        self.dsdump = dlist[0]
        self.npts = np.shape(dlist)[0]
        self.dlist = dlist
        self.budtw = tw
        if budday0 == 'default':
            self.budday0 = self.ddump.replace(day = 1)
        else:
            self.budday0 = budday0 
        self.dumptw = dumptw
        self.normalddump = self.budday0
        self.hist_dump_now = False
        self.just_dumped = False
        self.hist_count = 0
        self.hole_wait = 0
        self.hole = False
        self.drift_stats = drift_stats
        ### days in the future we're loading from
        ### typically 1 day for no time smoothing data
        self.hole_spec = np.arange(1,3)
        self.hole_width = 3
        #### history fields
        self.budcount = np.zeros([G.m,G.n],dtype = int)
        self.buddifhist = np.zeros([G.m,G.n])
        self.budadvhist = np.zeros([G.m,G.n])
        self.buddivhist = np.zeros([G.m,G.n])
        if self.split_div:
            ### two more arrays
            self.buddivplushist = np.zeros([G.m,G.n])
            self.buddivminshist = np.zeros([G.m,G.n])
        if self.transport:
            self.transporthist_x = np.zeros([G.m,G.n])
            self.transporthist_y = np.zeros([G.m,G.n])
        self.budreshist = np.zeros([G.m,G.n])
        self.buddynhist = np.zeros([G.m,G.n])
        if self.save_input:
            self.budthckhist = np.zeros([G.m,G.n])
            self.budconchist = np.zeros([G.m,G.n])
            self.budxvelhist = np.zeros([G.m,G.n])
            self.budyvelhist = np.zeros([G.m,G.n])
        if drift_stats:
            self.ddiv = np.zeros([G.m,G.n])
            self.dcrl = np.zeros([G.m,G.n])
            self.dshr = np.zeros([G.m,G.n])
            self.dshist_count = 0

    def check_dump(self,n,day,hole_verbos=False):
        ### check list for 
        self.hist_dump_now =False
        ### 1. end of list
        if n == self.npts - 2: 
            self.hist_dump_now = True
            return
        ### 2. holes (gap between data of more than a month
        if n<self.npts-3 and (self.dlist[n+2] - day).days > 32:
            self.hist_dump_now = True
            # if we hit a hole, we need to ingore the next n slices
            if hole_verbos: 
                print("hole found d1 "+day.strftime('%Y%m%d'))
                print("hole found d2 "+self.dlist[n+2].strftime('%Y%m%d'))
            self.hole = True
            self.hole_wait = 0
        ### 3. if we are a set no. of tw from dt0 
        ### change day.day == 1 to self.normalddump + self.dumptw == day
#         if day.day==1 and (not self.hole): 
        if self.normalddump + self.dumptw == day and (not self.hole): 
            self.hist_dump_now = True 
            self.normalddump = day
        
    def accumulate_budget(self,d,dif,adv,div,thck,conc,xvel,yvel,hole_verbos=False):
        ### every step
        ### accumulate history
        ### except holes
        if (self.hole_wait==self.hole_spec).any():
            print("------- avoiding hole "+d.strftime('%Y%m%d'))
            if hole_verbos: print('Waiting '+str(self.hole_wait))
        else:
            if hole_verbos: print('No hole '+d.strftime('%Y%m%d'))
            ### update to nan + 
            self.budcount      = np.sum(np.dstack((self.budcount,np.isfinite(dif))),2)
            self.buddifhist    = np.nansum(np.dstack((self.buddifhist,dif)),2)
            self.budadvhist    = np.nansum(np.dstack((self.budadvhist,adv)),2)
            self.buddivhist    = np.nansum(np.dstack((self.buddivhist,div)),2)
            if self.split_div:
                #### split the div to positive and negative parts
                div_split = copy(div)
                div_split[div_split<0] = 0.0
                self.buddivplushist    = np.nansum(np.dstack((self.buddivplushist,div_split)),2)
                div_split = copy(div)
                div_split[div_split>0] = 0.0
                self.buddivminshist    = np.nansum(np.dstack((self.buddivminshist,div_split)),2)
            if self.transport:
                #### calculate the transport
                transport_x = thck*conc*xvel
                transport_y = thck*conc*yvel
                self.transporthist_x = np.nansum(np.dstack((self.transporthist_x,transport_x)),2)
                self.transporthist_y = np.nansum(np.dstack((self.transporthist_y,transport_y)),2)
            self.budreshist    = np.nansum(np.dstack((self.budreshist,dif-adv-div)),2)
            self.buddynhist    = np.nansum(np.dstack((self.buddynhist,adv+div)),2)
            if self.save_input:
                self.budthckhist   = np.nansum(np.dstack((self.budthckhist,thck)),2)
                self.budconchist   = np.nansum(np.dstack((self.budconchist,conc)),2)
                self.budxvelhist   = np.nansum(np.dstack((self.budxvelhist,xvel)),2)
                self.budyvelhist   = np.nansum(np.dstack((self.budyvelhist,yvel)),2)
            self.hist_count += 1
        if self.hole: 
            self.hole_wait +=1
            if hole_verbos: print("in a hole "+d.strftime('%Y%m%d'))
        if self.hole_wait == self.hole_width:
            self.hole = False
            self.hole_wait = 0
            self.ddump = d
            self.dsdump = d
            ## update self.normalddump to be nearest previous timepoint
            nddump_search = True
            nextnddump = self.normalddump
            while nddump_search:
                if nextnddump + self.dumptw > d:
                    nddump_search = False
                    self.normalddump = nextnddump
                else:
                    nextnddump += self.dumptw
            print(" end of hole "+d.strftime('%Y%m%d'))
            
    def accumulate_dstats(self,ddiv,dcrl,dshr):
         ### every step
        ### accumulate history
        ### except holes
        if (self.hole_wait==np.array([1,2])).any():
#             print("------- avoiding hole "+d.strftime('%Y%m%d'))
            pass
        else:
            self.ddiv = np.nansum(np.dstack((self.ddiv,ddiv)),2)
            self.dcrl = np.nansum(np.dstack((self.dcrl,dcrl)),2)
            self.dshr = np.nansum(np.dstack((self.dshr,dshr)),2)
            self.dshist_count += 1
 
    def dump_pkl(self,spath,dn,G,dstats = False):
        # first day of the month
        nsecs = get_step_size(dn,self.budtw)
        print(' - - - - History dump: '+self.ddump.strftime('%Y%m%d'),flush=True)
        print('mean dif = '+str(np.nanmean(self.buddifhist*G.mask)*nsecs*self.hist_count),flush=True)
        print('mean adv = '+str(np.nanmean(self.budadvhist*G.mask)*nsecs*self.hist_count),flush=True)
        print('mean div = '+str(np.nanmean(self.buddivhist*G.mask)*nsecs*self.hist_count),flush=True)
        print('data points = '+str(np.sum(np.isfinite(self.buddivhist))),flush=True)
        #### save the budget
        budfields = {} #initialise dictionary
        budfields['budcount']     = self.budcount 
        budfields['buddifmonth']  = self.buddifhist 
        budfields['budadvmonth']  = self.budadvhist 
        budfields['buddivmonth']  = self.buddivhist 
        budfields['budresmonth']  = self.budreshist 
        budfields['buddynmonth']  = self.buddynhist 
        if self.save_input:
            budfields['budthckmonth'] = self.budthckhist/self.hist_count
            budfields['budconcmonth'] = self.budconchist/self.hist_count
            budfields['budxvelmonth'] = self.budxvelhist/self.hist_count
            budfields['budyvelmonth'] = self.budyvelhist/self.hist_count
        if dstats:
            budfields['driftdiv'] = self.ddiv/self.dshist_count
            budfields['driftcrl'] = self.dcrl/self.dshist_count
            budfields['driftshr'] = self.dshr/self.dshist_count

    #save
        f = open(spath+'budfields_'+self.ddump.strftime('%Y%m%d')+'.pkl','wb')
        pickle.dump(budfields,f)
        f.close()
        
        # empty fields
        self.budcount[:] = 0
        self.buddifhist[:] = 0.0
        self.budadvhist[:] = 0.0
        self.buddivhist[:] = 0.0
        self.budreshist[:] = 0.0
        self.buddynhist[:] = 0.0
        if self.save_input:
            self.budthckhist[:] = 0.0
            self.budconchist[:] = 0.0
            self.budxvelhist[:] = 0.0
            self.budyvelhist[:] = 0.0
        self.hist_count = 0
        if dstats:
            self.ddiv[:] = 0.0
            self.dcrl[:] = 0.0
            self.dshr[:] = 0.0
            self.dshist_count = 0
        self.ddump = dn
        

    def dump_nc(self,spath,dn,G,InputV,InputT,InputC,outrgd = False,accu_reset = True,extra_arrays = False,rgd_meth = 'nearest',print_info = False):
        """
        This is the history file write option.
        When this is called all the accumulated variables are written to file
        New options 2020-10:
            outrgd is a Gs2Gs option that will regrid the output onto a new grid
            outrdg = [Gs2Gs,'label',grid_set] where the 'label' is appended to the history file, and grid_set is the new grid
            accu_reset = True, this resests the accumulated variable to dump when writing to file. This is the default behaviour for a single type of history file
            Only set to False if you want to write to two types - say a regridded out too. The last written file will have to have accu_reset = True to stop madness
        New options 2021-2:
            extra_arrays is a list of arrays (in original shape) to add to the output (for example dump point thickness)
            each element in the the array is ['label',np.array]
            the label will be the Nc_file variable name
        """
        # first day of the month
        if type(outrgd) == bool:
            outrgd = [False,'']
        nsecs = get_step_size(dn,self.budtw)
        bud_scale  = 1.0
        if self.tot_budget: bud_scale = nsecs
        if print_info:
            print(' - - - - History dump: '+self.ddump.strftime('%Y%m%d'))
        if type(outrgd[0])==bool and print_info:
            print('mean dif = '+str(np.nanmean(self.buddifhist*G.mask)*nsecs),flush=True)
            print('mean adv = '+str(np.nanmean(self.budadvhist*G.mask)*nsecs),flush=True)
            print('mean div = '+str(np.nanmean(self.buddivhist*G.mask)*nsecs),flush=True)
            print('data points = '+str(np.sum(np.isfinite(self.buddivhist))),flush=True)
        elif print_info:
            print(outrgd[1]+' mean dif = '+str(np.nanmean(outrgd[0].rg_array(
                self.buddifhist,method=rgd_meth)*G.mask)*nsecs),flush=True)
            print(outrgd[1]+' mean adv = '+str(np.nanmean(outrgd[0].rg_array(
                self.budadvhist,method=rgd_meth)*G.mask)*nsecs),flush=True)
            print(outrgd[1]+' mean div = '+str(np.nanmean(outrgd[0].rg_array(
                self.buddivhist,method=rgd_meth)*G.mask)*nsecs),flush=True)
            
        file = spath+outrgd[1]+'budfields_'+self.ddump.strftime('%Y%m%d--')+dn.strftime('%Y%m%d')+'.nc'
        NC_f = Dataset(file, 'w', format='NETCDF4')
        NC_f.description = 'CPOM sea ice budget calculation outputs' 
        if type(outrgd[0])==bool:
            nc_m = G.m
            nc_n = G.n
        else:
            m,n = np.shape(outrgd[0].mesh_new[0])
            nc_m = m
            nc_n = n
            
#         NC_f.createDimension('time', 1)
        NC_f.createDimension('x', nc_m)
        NC_f.createDimension('y', nc_n)
        # to save:
        # A,swh,t0
        count=NC_f.createVariable('data_count', 'i4', ('x','y'))
        count.setncatts({
            'long_name':"no. of days of data present within each grid cell",
            'units':"count"}) 
        dif = NC_f.createVariable('intensification', 'f4', ('x','y'))
        dif.setncatts({'intensification:': "rate of change in sea ice thickness" ,
            'units':"meters" ,
            'comment':"rate of change in sea ice thickness (unit vol), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell." 
        })
        div = NC_f.createVariable('divergence', 'f4', ('x','y'))
        div.setncatts({
            'long_name':"change in sea ice thickness due to ice drift diveregence" ,
            'units':"meters" ,
            'comment':"divergence of the ice drift velocity field x ice thickness (unit volume), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell accountable from the diverging sea ice drift." 
            }        )
        if self.split_div:
            divplus = NC_f.createVariable('positive_divergence', 'f4', ('x','y'))
            divplus.setncatts({
            'long_name':"change in sea ice thickness due to ice drift diveregence, positive divergence only" ,
            'units':"meters" ,
            'comment':"divergence of the ice drift velocity field x ice thickness (unit volume), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell accountable from the positive diverging sea ice only sea ice drift." 
            }        )
            divmins = NC_f.createVariable('negative_divergence', 'f4', ('x','y'))
            divmins.setncatts({
            'long_name':"change in sea ice thickness due to ice drift diveregence, negative divergence only" ,
            'units':"meters" ,
            'comment':"divergence of the ice drift velocity field x ice thickness (unit volume), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell accountable from the negative diverging sea ice only sea ice drift." 
            }        )
        if self.transport:
            trnspt_x = NC_f.createVariable('transport_x', 'f4', ('x','y'))
            trnspt_x.setncatts({
            'long_name':"volume of sea ice drifting across a fixed line" ,
            'units':"meters" ,
            'comment':"Used only for gates, this grid is the value of concentration x thickness x (xvel) x no_of_seconds per history file period" 
            }        )
            trnspt_y = NC_f.createVariable('transport_y', 'f4', ('x','y'))
            trnspt_y.setncatts({
            'long_name':"volume of sea ice drifting across a fixed line" ,
            'units':"meters" ,
            'comment':"Used only for gates, this grid is the value of concentration x thickness x (yvel) x no_of_seconds per history file period" 
            }        )
        adv = NC_f.createVariable('advection', 'f4', ('x','y'))
        adv.setncatts({
            'long_name':"change in sea ice thickness due to advected sea ice" ,
            'units':"meters" ,
            'comment':"spatial variance in ice thickness x ice drift (unit volume), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell accountable due to the movement sea ice of variable thickness" 
            }        )
        dyn = NC_f.createVariable('dynamics', 'f4', ('x','y'))
        dyn.setncatts({
            'long_name':"change in sea ice thickness due to divergence and advection" ,
            'units':"meters" ,
            'comment':"the summed advection and divergence terms." 
            }        )
        res = NC_f.createVariable('residual', 'f4', ('x','y'))
        res.setncatts({
            'long_name':"Residual change in sea ice thickness (thermodynamic)" ,
            'units':"meters" ,
            'comment':"The rate of change is sea ice thickness (unit vol) that can not be accounted for from divergence and advection. intensificationc - dynmics. This value represents the total thermodynamic growth of sea within the 'Budget_period_length'." 
            }        )
        thk = NC_f.createVariable('thickness', 'f4', ('x','y'))
        thk.setncatts({
            'long_name':"averaged ice thickness on the destination grid" ,
            'units':"meters" ,
            'comment':"The averaged daily ice thickness used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " 
            }        )
        con = NC_f.createVariable('concentration', 'f4', ('x','y'))
        con.setncatts({
            'long_name':"averaged ice concentration on the destination grid" ,
            'units':"unit[0,1]" ,
            'comment':"The averaged daily ice concentration used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " 
            }        )
        xvl = NC_f.createVariable('ice_drift_x', 'f4', ('x','y'))
        xvl.setncatts({
            'long_name':"averaged ice drift x component" ,
            'units':"m/s" ,
            'comment':"The averaged daily ice drift x component used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " 
            }        )
        yvl = NC_f.createVariable('ice_drift_y', 'f4', ('x','y'))
        yvl.setncatts({
            'long_name':"averaged ice drift y component" ,
            'units':"m/s" ,
            'comment':"The averaged daily ice drift y component used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " 
            }        )
        # lat,lon,time(time in timestamp format)
        lat = NC_f.createVariable('lat', 'f4', ('x','y'))
        lat.setncatts({
            'long_name':"Latitude coordinate of grid cell" ,
            'units':"Degrees north" 
            }        )
        lon = NC_f.createVariable('lon', 'f4', ('x','y'))
        lon.setncatts({
            'long_name':"Longitude coordinate of grid cell" ,
            'units':"Degrees East" 
            }        )
        ### additional attributes
        if type(extra_arrays)==bool:
            pass
        else:
            nc_ea = []
            for ea in extra_arrays:
                nc_ea.append(NC_f.createVariable(ea[0], 'f4', ('x','y')))
                
        
#         NC_f.setncattr_string('Ice_drift_data_source',InputV.name)
#         NC_f.setncattr_string('Ice_thickness_data_source',InputT.name)
#         NC_f.setncattr_string('Ice_concentration_data_source',InputC.name)
#         NC_f.setncattr_string('Created_on',dt.date.today().strftime('%Y-%m-%d'))
#         NC_f.setncattr_string('Budget period start',self.ddump.strftime('%Y-%m-%d'))
#         NC_f.setncattr_string('Budget period end',dn.strftime('%Y-%m-%d'))
#         NC_f.setncattr_string('Budget period length (days)',(dn - self.ddump).days)
#         NC_f.setncattr_string('Budget time slices',self.hist_count)
        NC_f.setncatts({
            'description':"CPOM sea ice budget calculation outputs" ,
            'comment':"Budget processing code by H. Heorton, code available at https://github.com/CPOMUCL/Budget_tool/" ,
        'Ice_drift_data_source':InputV.name,
        'Ice_thickness_data_source':InputT.name,
        'Ice_concentration_data_source':InputC.name,
        'Created_on':dt.date.today().strftime('%Y-%m-%d'),
        'Budget period start':self.ddump.strftime('%Y-%m-%d'),
        'Budget period end':dn.strftime('%Y-%m-%d'),
        'Budget period length (days)':str((dn - self.ddump).days),
        'Budget time slices':self.hist_count,
            'publisher_name':"UCL_CPOM" ,
            'publisher_type':"institution" ,
            'publisher_email':"a.muir@ucl.ac.uk" ,
            'publisher_url':"http://www.cpom.ucl.ac.uk/csopr/" 
        })
        
        self.buddifhist[self.budcount==0] = np.nan
        self.buddivhist[self.budcount==0] = np.nan
        if self.split_div:
            self.buddivplushist[self.budcount==0] = np.nan
            self.buddivminshist[self.budcount==0] = np.nan
        if self.transport:
            self.transporthist_x[self.budcount==0] = np.nan
            self.transporthist_y[self.budcount==0] = np.nan
        self.budadvhist[self.budcount==0] = np.nan
        self.buddynhist[self.budcount==0] = np.nan
        self.budreshist[self.budcount==0] = np.nan
        if self.save_input:
            self.budthckhist[self.budcount==0] = np.nan
            self.budconchist[self.budcount==0] = np.nan
            self.budxvelhist[self.budcount==0] = np.nan
            self.budyvelhist[self.budcount==0] = np.nan
        if type(outrgd[0])==bool:
            count[:]=self.budcount
            dif[:] = self.buddifhist*bud_scale
            div[:] = self.buddivhist*bud_scale
            if self.split_div:
                divplus[:] = self.buddivplushist*bud_scale
                divmins[:] = self.buddivminshist*bud_scale
            if self.transport:
                trnspt_x[:] = self.transporthist_x*bud_scale
                trnspt_y[:] = self.transporthist_y*bud_scale
            adv[:] = self.budadvhist*bud_scale
            dyn[:] = self.buddynhist*bud_scale
            res[:] = self.budreshist*bud_scale
            if self.save_input:
                thk[:] = self.budthckhist/self.hist_count
                con[:] = self.budconchist/self.hist_count
                xvl[:] = self.budxvelhist/self.hist_count
                yvl[:] = self.budyvelhist/self.hist_count
            lon[:] = G.lons
            lat[:] = G.lats
            
        else:
            count[:]=outrgd[0].rg_array(self.budcount,method=rgd_meth)
            dif[:] = outrgd[0].rg_array(self.buddifhist*bud_scale,method=rgd_meth)
            div[:] = outrgd[0].rg_array(self.buddivhist*bud_scale,method=rgd_meth)
            if self.split_div:
                divplus[:] = outrgd[0].rg_array(self.buddivplushist*bud_scale,method=rgd_meth)
                divmins[:] = outrgd[0].rg_array(self.buddivminshist*bud_scale,method=rgd_meth)
            if self.transport:
                xtemp,ytemp = outrgd[0].rg_vecs(self.transporthist_x*bud_scale,
                                       self.transporthist_y*bud_scale,method=rgd_meth)
                trnspt_x[:] = xtemp
                trnspt_y[:] = ytemp
            adv[:] = outrgd[0].rg_array(self.budadvhist*bud_scale,method=rgd_meth)
            dyn[:] = outrgd[0].rg_array(self.buddynhist*bud_scale,method=rgd_meth)
            res[:] = outrgd[0].rg_array(self.budreshist*bud_scale,method=rgd_meth)
            if self.save_input:
                thk[:] = outrgd[0].rg_array(self.budthckhist/self.hist_count,method=rgd_meth)
                con[:] = outrgd[0].rg_array(self.budconchist/self.hist_count,method=rgd_meth)
                xtemp,ytemp = outrgd[0].rg_vecs(self.budxvelhist/self.hist_count,
                                           self.budyvelhist/self.hist_count,method=rgd_meth)
                xvl[:] = xtemp
                yvl[:] = ytemp
            lon[:] = outrgd[2].lons
            lat[:] = outrgd[2].lats
        if type(extra_arrays)==bool:
            pass
        else:
            for ea,nc_v in zip(extra_arrays,nc_ea):
                if type(outrgd[0])==bool:
                    nc_v[:] = ea[1]
                else:
                    nc_v[:] = outrgd[0].rg_array(ea[1],method=rgd_meth)
        NC_f.close()
        # empty fields
        if accu_reset:
            self.budcount[:] = 0
            self.buddifhist[:] = 0.0
            self.budadvhist[:] = 0.0
            self.buddivhist[:] = 0.0
            if self.split_div:
                self.buddivplushist[:] = 0.0
                self.buddivminshist[:] = 0.0
            if self.transport:
                self.transporthist_x[:] = 0.0
                self.transporthist_y[:] = 0.0
            self.budreshist[:] = 0.0
            self.buddynhist[:] = 0.0
            if self.save_input:
                self.budthckhist[:] = 0.0
                self.budconchist[:] = 0.0
                self.budxvelhist[:] = 0.0
                self.budyvelhist[:] = 0.0
            self.hist_count = 0
            self.ddump = dn
            
            
##### NETCDF budget outputs for help
# netcdf budfields_20210101 {
# dimensions:
#         time = 1 ;
#         x = 119 ;
#         y = 177 ;
# variables:
#         int data_count(time, x, y) ;
#                 data_count:long_name = "no. of days of data present within each grid cell" ;
#                 data_count:units = "count" ;
#         float intensification(time, x, y) ;
#                 intensification:long_name = "rate of change in sea ice thickness" ;
#                 intensification:units = "meters" ;
#                 intensification:comment = "rate of change in sea ice thickness (unit vol), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell." ;
#         float divergence(time, x, y) ;
#                 divergence:long_name = "change in sea ice thickness due to ice drift diveregence" ;
#                 divergence:units = "meters" ;
#                 divergence:comment = "divergence of the ice drift velocity field x ice thickness (unit volume), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell accountable from the diverging sea ice drift." ;
#         float advection(time, x, y) ;
#                 advection:long_name = "change in sea ice thickness due to advected sea ice" ;
#                 advection:units = "meters" ;
#                 advection:comment = "spatial variance in ice thickness x ice drift (unit volume), calculated on a daily timescale, accumulated over the 'Budget_period_length'. This value shows the difference in sea ice thickness (unit vol) within one static grid cell accountable due to the movement sea ice of variable thickness" ;
#         float dynamics(time, x, y) ;
#                 dynamics:long_name = "change in sea ice thickness due to divergence and advection" ;
#                 dynamics:units = "meters" ;
#                 dynamics:comment = "the summed advection and divergence terms." ;
#         float residual(time, x, y) ;
#                 residual:long_name = "Residual change in sea ice thickness (thermodynamic)" ;
#                 residual:units = "meters" ;

#                 residual:comment = "The rate of change is sea ice thickness (unit vol) that can not be accounted for from divergence and advection. intensificationc - dynmics. This value represents the total thermodynamic growth of sea within the 'Budget_period_length'." ;
#         float thickness(time, x, y) ;
#                 thickness:long_name = "averaged ice thickness on the destination grid" ;
#                 thickness:units = "meters" ;
#                 thickness:comment = "The averaged daily ice thickness used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " ;
#         float concentration(time, x, y) ;
#                 concentration:long_name = "averaged ice concentration on the destination grid" ;
#                 concentration:units = "unit[0,1]" ;
#                 concentration:comment = "The averaged daily ice concentration used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " ;
#         float ice_drift_x(time, x, y) ;
#                 ice_drift_x:long_name = "averaged ice drift x component" ;
#  float ice_drift_y(time, x, y) ;
#                 ice_drift_y:long_name = "averaged ice drift y component" ;
#                 ice_drift_y:units = "m/s" ;
#                 ice_drift_y:comment = "The averaged daily ice drift y component used within the budget calculation. Accumulated over the 'Budget_period_length' and averaged. " ;
#         float lat(time, x, y) ;
#                 lat:long_name = "Latitude coordinate of grid cell" ;
#                 lat:units = "Degrees north" ;
#         float lon(time, x, y) ;
#                 lon:long_name = "Longitude coordinate of grid cell" ;
#                 lon:units = "Degrees East" ;

# // global attributes:
#                 :description = "CPOM sea ice budget calculation outputs" ;
#                 :comment = "Budget processing code by H. Heorton, code available at https://github.com/CPOMUCL/Budget_tool/" ;
#                 :Budget_destination_grid = "osisaf" ;
#                 :Ice_drift_data_source = "osisaf" ;
#                 :Ice_thickness_data_source = "AWI_SMOS_daily" ;
#                 :Ice_concentration_data_source = "NSIDC_n" ;
#                 :Created_on = "2021-03-11" ;
#                 :Budget_period_start = "2021-01-01" ;
#                 :Budget_period_end = "2021-02-01" ;
#                 :Budget_period_length_days = 31_days ;
#                 :Budget_time_slices = 31 ;
#                 :publisher_name = "UCL_CPOM" ;
#                 :publisher_type = "institution" ;
#                 :publisher_email = "a.muir@ucl.ac.uk" ;
#                 :publisher_url = "http://www.cpom.ucl.ac.uk/csopr/" ;
# }}

    def dump_dstats_nc(self,spath,dn,G,InputV):
        # first day of the month
        nsecs = get_step_size(dn,self.budtw)
        print(' - - - - Dstats dump: '+self.dsdump.strftime('%Y%m%d'))
        print('mean div = '+str(np.nanmean(self.ddiv*G.mask)*nsecs*self.dshist_count),flush=True)
        print('mean crl = '+str(np.nanmean(self.dcrl*G.mask)*nsecs*self.dshist_count),flush=True)
        print('mean shr = '+str(np.nanmean(self.dshr*G.mask)*nsecs*self.dshist_count),flush=True)
        print('data points = '+str(np.sum(np.isfinite(self.ddiv))),flush=True)
        file = spath+InputV.name+'_stats_'+self.dsdump.strftime('%Y%m%d')+'.nc'
        NC_f = Dataset(file, 'w', format='NETCDF4')
        NC_f.description = 'CPOM sea ice budget calculation outputs, drift statistics' 
            
        NC_f.createDimension('time', 1)
        NC_f.createDimension('x', G.m)
        NC_f.createDimension('y', G.n)
        # to save:
        # A,swh,t0
        xvl = NC_f.createVariable('ice_drift_x', 'f4', ('time','x','y'))
        yvl = NC_f.createVariable('ice_drift_y', 'f4', ('time','x','y'))
        # lat,lon,time(time in timestamp format)
        lat = NC_f.createVariable('lat', 'f4', ('time','x','y'))
        lon = NC_f.createVariable('lon', 'f4', ('time','x','y'))
        ddv = NC_f.createVariable('ice_divergence', 'f4', ('time','x','y'))
        dcl = NC_f.createVariable('ice_curl', 'f4', ('time','x','y'))
        dsh = NC_f.createVariable('ice_shear', 'f4', ('time','x','y'))
        ### additional attributes
        
        NC_f.setncattr_string('Ice_drift_data_source',InputV.name)
        NC_f.setncattr_string('Created_on',dt.date.today().strftime('%Y-%m-%d'))
        xvl[:] = self.budxvelhist/self.dshist_count
        yvl[:] = self.budyvelhist/self.dshist_count
        ddv[:] = self.ddiv/self.dshist_count
        dcl[:] = self.dcrl/self.dshist_count
        dsh[:] = self.dshr/self.dshist_count
        lon[:] = G.lons
        lat[:] = G.lats
        NC_f.close()
        # empty fields
        self.ddiv[:] = 0.0
        self.dcrl[:] = 0.0
        self.dshr[:] = 0.0
        self.dshist_count = 0
        self.dsdump = dn
        
def plot_inputs(Gvel,Gthk,Gcon,
                InputV,InputT,InputC,
                Utemp,Vtemp,Tplot,Cplot,spath,ts):
    """
    Diagnostic figures for the input data - projected properly
    All on their native grids
    """
    f = plt.figure(figsize=[9,3])
    # dtarg = dt.datetime(2015,1,1)
    # t_p1 = B1.condy.get_index(dtarg)
    # t_p2 = B2.condy.get_index(dtarg)
    # t_p3 = B3.condy.get_index(dtarg)
    vlim = 0.4
    a_no = 30
#     a_sc = 3.9e-1
    a_sc = 2.5*np.nanpercentile(np.hypot(Utemp,Vtemp),[90])[0]
#     print(a_sc)
    m = Gvel.mplot

    plt.subplot(1,3,1)
    m.pcolormesh(Gcon.xptp,Gcon.yptp,Cplot,rasterized=True)
    m.colorbar(location='bottom')
    plt.clim([0.0,1.0])
    m.drawcoastlines()
    plt.ylabel(InputC.name)
    plt.title('Ice conc. '+ts.strftime('%Y%m%d'))

    plt.subplot(1,3,2)
    m.pcolormesh(Gthk.xptp,Gthk.yptp,Tplot,rasterized=True)
    m.colorbar(location='bottom')
    plt.clim([0.0,5.0])
    m.drawcoastlines()
    plt.ylabel(InputT.name)
    plt.title('Ice thick. '+ts.strftime('%Y%m%d'))

    plt.subplot(1,3,3)
    rm = int(Gvel.m/a_no)
    rn = int(Gvel.n/a_no)
    ra = np.sqrt(rm+rn)
    ra=ra*a_sc
    ur,vr = Gvel.rotate_vectors_to_plot(Utemp,Vtemp)
    
    m.pcolormesh(Gvel.xptp,Gvel.yptp,np.hypot(Utemp,Vtemp),rasterized=True)
    m.colorbar(location='bottom')
    plt.clim([0.0,0.3])
    m.drawcoastlines()
    m.quiver(Gvel.xpts[::rm,::rn],Gvel.ypts[::rm,::rn],
             ur[::rm,::rn],vr[::rm,::rn],scale = ra,width=0.005)
    plt.ylabel(InputV.name)
    plt.title('Ice drift. '+ts.strftime('%Y%m%d'))
    f.savefig(spath+'Inputs_'+ts.strftime('%Y%m%d')+'.pdf',
              bbox_inches='tight')
    print('Saving figure: '+spath+'Inputs_'+ts.strftime('%Y%m%d')+'.pdf')


def plot_budget(Gvel,buddif,buddiv3,budadv3,Utemp,Vtemp,spath,ts):
    """
    Diagnostic figures for the budget - projected properly
    """
    f = plt.figure(figsize=[9,9])
    # dtarg = dt.datetime(2015,1,1)
    # t_p1 = B1.condy.get_index(dtarg)
    # t_p2 = B2.condy.get_index(dtarg)
    # t_p3 = B3.condy.get_index(dtarg)
    vlim = 0.4
    a_no = 30
#     a_sc = 3.9e-1
    # a_sc = 2*np.nanmax(np.hypot(Utemp,Vtemp))
#     print(a_sc)
    a_sc = 2.5*np.nanpercentile(np.hypot(Utemp,Vtemp),[90])[0]
    m = Gvel.mplot
    p_rng = np.nanpercentile(buddif,[2,98])
    pr = np.max(np.abs(p_rng))
    p_rng = [-pr,pr]

    ### intensification
    plt.subplot(2,2,1)
    m.pcolormesh(Gvel.xptp,Gvel.yptp,buddif,cmap='RdBu',rasterized=True)
    m.colorbar(location='bottom')
    plt.clim(p_rng)
    m.drawcoastlines()
    plt.title('Intensification '+ts.strftime('%Y%m%d'))

    ### DIVERGENCE
    plt.subplot(2,2,2)
    rm = int(Gvel.m/a_no)
    rn = int(Gvel.n/a_no)
    ra = np.sqrt(rm+rn)
    ra=ra*a_sc
    m.pcolormesh(Gvel.xptp,Gvel.yptp,buddiv3,cmap='RdBu',rasterized=True)
    plt.clim(p_rng)
    m.colorbar(location='bottom')
    ur,vr = Gvel.rotate_vectors_to_plot(Utemp,Vtemp)
    m.quiver(Gvel.xpts[::rm,::rn],Gvel.ypts[::rm,::rn],
             ur[::rm,::rn],vr[::rm,::rn],scale = ra,width=0.005)
    m.drawcoastlines()
    plt.title('Divergence '+ts.strftime('%Y%m%d'))

    ### ADVECTION
    plt.subplot(2,2,3)
    rm = int(Gvel.m/a_no)
    rn = int(Gvel.n/a_no)
    ra = np.sqrt(rm+rn)
    ra=ra*a_sc
    m.pcolormesh(Gvel.xptp,Gvel.yptp,budadv3,cmap='RdBu',rasterized=True)
    plt.clim(p_rng)
    m.colorbar(location='bottom')
    m.quiver(Gvel.xpts[::rm,::rn],Gvel.ypts[::rm,::rn],
             ur[::rm,::rn],vr[::rm,::rn],scale = ra,width=0.005)
    m.drawcoastlines()
    plt.title('Advection '+ts.strftime('%Y%m%d'))
    
    ### intensification
    plt.subplot(2,2,4)
    m.pcolormesh(Gvel.xptp,Gvel.yptp,buddif-buddiv3-budadv3,cmap='RdBu',rasterized=True)
    plt.clim(p_rng)
    m.colorbar(location='bottom')
    m.drawcoastlines()
    plt.title('Residual '+ts.strftime('%Y%m%d'))

    
    f.savefig(spath+'Budget_components_'+ts.strftime('%Y%m%d')+'.pdf',
              bbox_inches='tight')
    print('Saving figure: '+spath+'Budget_components_'+ts.strftime('%Y%m%d')+'.pdf')
    


def plot_budget_square(Gvel,buddif,buddiv3,budadv3,Utemp,Vtemp,spath,ts):
    """
    Diagnostic figures for the budget - square grid
    """
    f = plt.figure(figsize=[9,9])
    # dtarg = dt.datetime(2015,1,1)
    # t_p1 = B1.condy.get_index(dtarg)
    # t_p2 = B2.condy.get_index(dtarg)
    # t_p3 = B3.condy.get_index(dtarg)
    vlim = 0.4
    a_no = 20
#     a_sc = 3.9e-1
#     a_sc = 2*np.nanmax(np.hypot(Utemp,Vtemp))
    a_sc = 2.5*np.nanpercentile(np.hypot(Utemp,Vtemp),[90])[0]
#     print(a_sc)
    p_rng = np.nanpercentile(buddif,[2,98])
    pr = np.max(np.abs(p_rng))
    p_rng = [-pr,pr]
    
    Gvel.get_square_points()

    ### intensification
    plt.subplot(2,2,1)
    plt.pcolormesh(Gvel.xsq,Gvel.ysq,buddif,cmap='RdBu',rasterized=True)
    plt.colorbar(orientation="horizontal")
    plt.clim(p_rng)
    plt.title('Intensification '+ts.strftime('%Y%m%d'))

    ### DIVERGENCE
    plt.subplot(2,2,2)
    rm = int(Gvel.m/a_no)
    rn = int(Gvel.n/a_no)
    ra = np.sqrt(rm+rn)
    ra=ra*a_sc
    plt.pcolormesh(Gvel.xsq,Gvel.ysq,buddiv3,cmap='RdBu',rasterized=True)
    plt.clim(p_rng)
    plt.colorbar(orientation="horizontal")
    plt.quiver(Gvel.xsq[::rm,::rn],Gvel.ysq[::rm,::rn],
             Utemp[::rm,::rn],Vtemp[::rm,::rn],scale = ra,width=0.005)
    plt.title('Divergence '+ts.strftime('%Y%m%d'))

    ### ADVECTION
    plt.subplot(2,2,3)
    plt.pcolormesh(Gvel.xsq,Gvel.ysq,budadv3,cmap='RdBu',rasterized=True)
    plt.clim(p_rng)
    plt.colorbar(orientation="horizontal")
    plt.quiver(Gvel.xsq[::rm,::rn],Gvel.ysq[::rm,::rn],
             Utemp[::rm,::rn],Vtemp[::rm,::rn],scale = ra,width=0.005)
    plt.title('Advection '+ts.strftime('%Y%m%d'))
    
    ### intensification
    plt.subplot(2,2,4)
    plt.pcolormesh(Gvel.xsq,Gvel.ysq,buddif-buddiv3-budadv3,cmap='RdBu',rasterized=True)
    plt.clim(p_rng)
    plt.colorbar(orientation="horizontal")
    plt.title('Residual '+ts.strftime('%Y%m%d'))

    
    f.savefig(spath+'Budget_components_square_'+ts.strftime('%Y%m%d')+'.pdf',
              bbox_inches='tight')
    print('Saving figure: '+spath+'Budget_components_square_'+ts.strftime('%Y%m%d')+'.pdf')


def plot_dstats(Gvel,Utemp,Vtemp,InputV,ddiv,dcrl,dshr,spath,ts):
    """
    Diagnostic figures for the drift statistics - projected properly
    """
    f = plt.figure(figsize=[9,3])
    # dtarg = dt.datetime(2015,1,1)
    # t_p1 = B1.condy.get_index(dtarg)
    # t_p2 = B2.condy.get_index(dtarg)
    # t_p3 = B3.condy.get_index(dtarg)
    vlim = 0.4
    a_no = 30
#     a_sc = 3.9e-1
#     a_sc = 2*np.nanmax(np.hypot(Utemp,Vtemp))
    a_sc = 2*np.nanpercentile(np.hypot(Utemp,Vtemp),[90])[0]
    rm = int(Gvel.m/a_no)
    rn = int(Gvel.n/a_no)
    ra = np.sqrt(rm+rn)
    ra=ra*a_sc
    m = Gvel.mplot

    plt.subplot(1,3,1)
    m.pcolormesh(Gvel.xptp,Gvel.yptp,ddiv,cmap='RdBu',rasterized=True)
    m.colorbar(location='bottom')
    p_rng = np.nanpercentile(ddiv,[40,60])
    pr = np.max(np.abs(p_rng))
    p_rng = [-pr,pr]
#     plt.clim(p_rng)
#     plt.clim([0.0,1.0])
    m.drawcoastlines()
    ur,vr = Gvel.rotate_vectors_to_plot(Utemp,Vtemp)
    m.quiver(Gvel.xpts[::rm,::rn],Gvel.ypts[::rm,::rn],
             ur[::rm,::rn],vr[::rm,::rn],scale = ra,width=0.005)
    plt.ylabel(InputV.name)
    plt.title('Drift div. '+ts.strftime('%Y%m%d'))

    plt.subplot(1,3,2)
    m.pcolormesh(Gvel.xptp,Gvel.yptp,dcrl,cmap='RdBu',rasterized=True)
    m.colorbar(location='bottom')
    p_rng = np.nanpercentile(dcrl,[8,92])
    pr = np.max(np.abs(p_rng))
    p_rng = [-pr,pr]
    plt.clim(p_rng)
#     plt.clim([0.0,5.0])
    m.drawcoastlines()
    m.quiver(Gvel.xpts[::rm,::rn],Gvel.ypts[::rm,::rn],
             ur[::rm,::rn],vr[::rm,::rn],scale = ra,width=0.005)
    plt.ylabel(InputV.name)
    plt.title('Drift curl '+ts.strftime('%Y%m%d'))

    plt.subplot(1,3,3)
    rm = int(Gvel.m/a_no)
    rn = int(Gvel.n/a_no)
    ra = np.sqrt(rm+rn)
    ra=ra*a_sc
    
    m.pcolormesh(Gvel.xptp,Gvel.yptp,dshr,cmap='YlGnBu',rasterized=True)
    m.colorbar(location='bottom')
    p_rng = np.nanpercentile(dshr,[0,87])
    plt.clim(p_rng)
#     plt.clim([0.0,0.3])
    m.drawcoastlines()
    m.quiver(Gvel.xpts[::rm,::rn],Gvel.ypts[::rm,::rn],
             ur[::rm,::rn],vr[::rm,::rn],scale = ra,width=0.005)
    plt.ylabel(InputV.name)
    plt.title('Drift shear '+ts.strftime('%Y%m%d'))
    f.savefig(spath+'Drift_statistics_'+ts.strftime('%Y%m%d')+'.pdf',
              bbox_inches='tight')
    print('Saving figure: '+spath+'Drift_statistics_'+ts.strftime('%Y%m%d')+'.pdf')
