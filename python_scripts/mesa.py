''' MESA output data loading and plotting

    v0.2, 15OCT2012: NuGrid collaboration (Sam Jones, Michael Bennett, Daniel Conti,
                     William Hillary, Falk Herwig, Christian Ritter)
    v0.1, 23JUN2010: Falk Herwig 
    

    mesa.py provides tools to get MESA stellar evolution data output
    into your favourite python session. In the LOGS directory MESA
    outputs two types of files: history.data or star.log is a time
    evolution output, printing one line per so many cycles (e.g. each
    cycle) of all sorts of things. profilennn.data or lognnn.data
    files are profile data files.  nnn is the number of profile.data
    or log.data files that is translated into model cycles in the
    profiles.index file.

    MESA allows users to freely define what should go into these two
    types of outputs, which means that column numbers can and do
    change. mesa.py reads in both types of files and present them (as
    well as any header attributes) as arrays that can be referenced by
    the actual column name as defined in the header section of the
    files. mesa.py then defines a (hopefully growing) set of standard
    plots that make use of the data just obtained.

    mesa.py is organised as a module that can be imported into any
    python or ipython session. It is related to nmuh5.py which is a
    similar module to deal with 'se' output, used by the NuGrid
    collaboration. mesa.py does not need se libraries. The 'se' output
    files that can be written with MESA can be read and processed with
    the nuh5.py tool.

    mesa.py is providing two class objects, mesa_profile and history_data.
    The first makes profile data available, the second reads and plots the
    history.dataor star.log file. Note that several instances of these can be
    initiated within one session and data from different instances
    (i.e. models, tracks etc) can be overplotted.

    Here is how a simple session could look like that is plotting an
    HRD (We prefer to load ipython with matplotlib and numpy support
    via the alias
    alias mpython='ipython -pylab -p numpy -editor emacsclient')

        vortex$ mpython
        Python 2.5.2 (r252:60911, Feb 22 2008, 07:57:53) 
        Type "copyright", "credits" or "license" for more information.

        IPython 0.9.1 -- An enhanced Interactive Python.
        ?         -> Introduction and overview of IPython's features.
        %quickref -> Quick reference.
        help      -> Python's own help system.
        object?   -> Details about 'object'. ?object also works, ?? prints more.

        IPython profile: numpy

          Welcome to pylab, a matplotlib-based Python environment.
          For more information, type 'help(pylab)'.

        In [1]: import mesa as ms

        In [2]: help ms
        ------> help(ms)

        In [4]: s=ms.history_data('.')

        In [5]: s.hrd()
    
     In order to find out what header attributes and columns are
     available in history.data or star.log use:
     
        In [6]: s.header_attr
        Out[6]: 
        {'burn_min1': 50.0,
         'burn_min2': 1000.0,
         'c12_boundary_limit': 0.0001,
         'h1_boundary_limit': 0.0001,
         'he4_boundary_limit': 0.0001,
         'initial_mass': 2.0,
         'initial_z': 0.01}

        In [7]: s.cols
        Out[7]: 
        {'center_c12': 38,
         'center_h1': 36,
         'center_he4': 37,
          ...  
          
    In order to read the profile data from the first profile.data file
    in profiles.index, and then get the mass and temperature out and
    finally plot them try. Typically you will have already a
    Kippenhahn diagram as a function of model number in front of you,
    and you want to access profile information for a given cycle
    number. Typically you do not have profiles for all cycle
    numbers. The best way to start a profile instance is with
    num_type='nearest_model' (check the docstring for other ways to
    select profiles for a profile instance):

        In [9]: a1=ms.mesa_profile('LOGS',59070,,num_type='nearest_model')
        2001 in profiles.index file ...
        Found and load nearest profile for cycle 59000
        reading LOGS/profile1801.data ...
        Closing profile tool ...

        In [10]: T=a1.get('temperature')

        In [11]: mass=a1.get('mmid')

        In [12]: plot(mass,T)
        Out[12]: [<matplotlib.lines.Line2D object at 0x8456ed0>]

    Or, you could have had it easier in the following way:
        In [13]: a1.plot('mass','c12',logy=True,shape='-',legend='$^{12}\mathrm{C}$')
    where the superclass plot method interprets data column headers
    correctly and does all the work for you.

    Of course, a1.cols etc are available here as well and many other
    things. E.g. a.model contains an array with all the models for
    which profile.data or log.data  are available. You may initiate a profile object
    with a model number:

        In [14]: a2=ms.mesa_profile('.',55000,num_type='model')
        100 in profiles.index file ...
        reading ./profile87.data ...

   a1.log_ind (for any profile instance) provides a map of model
   number to profile file number. 
   a1.cols and a1.header_attr gives the column names and header attributes.

   

'''
import numpy as np
from data_plot import *
import numpy as np
import matplotlib
from pylab import *
import pylab as pyl
import matplotlib.pyplot as pl
import os
import sys



class mesa_profile(DataPlot):
    ''' read profiles.index and prepare reading MESA profile files

    starts with reading profiles.index and creates hash array
    profile.data can then be accessed via prof_plot
    
    '''

    sldir = ''



    
    def __init__(self,sldir,num,num_type='nearest_model',prof_ind_name='profiles.index',profile_prefix='profile',data_suffix='.data'):
        '''read a profile.data profile file

        input:
        sldir       directory path of LOGS

        num         by default this is the i. profile file (profile.data or log.data) available
                    (e.g. num=1 is the 1. available profile file),
                    however if you give 

        num_type    'model' (exact) or 'nearest_model': get the profile
                    profile.data file for model (or cycle number) used
                    by the stellar evolution code

                    'profile_num': num will be interpreted as the
                    profile.data or log.data number profile_num
                    (profile_num is the number that appears in the
                    file names of type profile23.data or log23.data)

                    'profiles_i': the ith file in profiles.index file
        prof_ind_name    use this optional argument if the profiles.index 
                         file hasn an alternative name, for example, do 
                         superpro=ms.profile('LOGS',1,prof_ind_name='super.prof') 
        log_prefix  Depending on what mesa version you use, for the case of mesa
		    version before 4442 and no log.data file is found, 
		    log_prefix is internal changed to 'log' for using log#.data files
	data_suffix are optional arguments that allow you to change
                    the defaults for the profile.data profile files. '''

        self.prof_ind_name = prof_ind_name
        self.sldir         = sldir

        if num_type is 'nearest_model' or num_type is 'model':
            self.profiles_index()
        if num_type is 'nearest_model':
            amods=array(self.model)
            nearmods=[where(amods<=num)[0][-1],where(amods>=num)[0][0]]
            sometable={}
            for thing in nearmods:
                sometable[abs(self.model[thing]-num)]=thing
            nearest       = min(abs(self.model[nearmods[0]]-num),\
                                    abs(self.model[nearmods[1]]-num))
            num = self.model[sometable[nearest]] 
            print 'Found and load nearest profile for cycle '+str(num)
            num_type = 'model'            
        if num_type is 'model':
            try:
                log_num=self.log_ind[num]
            except KeyError:
                print 'There is no profile file for this model'
                print "You may retry with num_type='nearest_model'"
                return
        elif num_type is 'profiles_i':
            log_num=self.log_file_ind(num)
            if log_num == -1:
                print "Could not find a profile file with that number"
                return
        elif num_type is 'profile_num':
            log_num = num
        else:
            print 'unknown num_type'
            return

        filename=self.sldir+'/'+profile_prefix+str(log_num)+data_suffix
	if not os.path.exists(filename):
		profile_prefix='log'			
		filename=self.sldir+'/'+profile_prefix+str(log_num)+data_suffix
		if not os.path.exists(filename):
			print 'error: no profile.data file found in '+sldir
			print 'error: no log.data file found in '+sldir       
			
 
        print 'reading '+filename+' ...'
        header_attr = read_mesafile(filename,only='header_attr')
        num_zones=int(header_attr['num_zones'])
        header_attr,cols,data = read_mesafile(filename,data_rows=num_zones,only='all')

        self.cols        = cols
        self.header_attr = header_attr
        self.data        = data


    def __del__(self):
        print 'Closing profile tool ...'

    def profiles_index(self):
        ''' read profiles.index and make hash array

        log_ind     hash array that returns profile.data or log.data file number from model number
        model       the models for which profile.data or log.data is available'''

        prof_ind_name = self.prof_ind_name 

        f = open(self.sldir+'/'+prof_ind_name,'r')
        line = f.readline()
        numlines=int(line.split()[0])
        print str(numlines)+' in profiles.index file ...'

        model=[]
        log_file_num=[]
        for line in f:
            model.append(int(line.split()[0]))
            log_file_num.append(int(line.split()[2]))

        log_ind={}    # profile.data number from model
        for a,b in zip(model,log_file_num):
            log_ind[a] = b
            
        self.log_ind=log_ind
        self.model=model

# let's start with functions that aquire data

    def log_file_ind(self,inum):
        ''' information about available profile.data or log.data files
        
        inmu       attempt to get number of inum's profile.data file
        inum_max   max number of profile.data or log.data files available'''
        
        self.profiles_index()
        if inum <= 0:
            print "Smallest argument is 1"
            return

        inum_max = len(self.log_ind)
        inum -= 1
        
        if inum > inum_max:
            print 'There are only '+str(inum_max)+' profile file available.'
            log_data_number = -1
            return log_data_number
        else:
            log_data_number=self.log_ind[self.model[inum]]
            print 'The '+str(inum+1)+'. profile.data file is '+ \
                  str(log_data_number)
            return log_data_number

    def get(self,str_name):
        ''' return a column of data with the name str_name
        
        str_name is the name of the column as printed in the
        profilennn.data or lognnn.data file; get the available columns from self.cols
        (where you replace self with the name of your instance)'''

        column_array = self.data[:,self.cols[str_name]-1].astype('float')
        return column_array



        
class history_data(DataPlot):
    ''' read history.data or star.log MESA output and plot various things, including
    HRD, Kippenhahn etc
    
    sldir              - which LOGS directory
    slname='history.data'  - if star.log is available instead, star.log file is read,
			     this is an optional argument if history.data or star.log 
			     file has an alternative name,
    clean_starlog=True - request new cleaning of history.data or star.log, makes 
			 history.datasa or star.logsa which is the file that 
			 is actually read and plotted
    use like this: another=ms.history_data('LOGS',slname='anothername')
    '''

    sldir  = ''
    slname = ''
    header_attr = []
    cols = [] 
    
    def __init__(self,sldir,slname='history.data',clean_starlog=True):
        self.sldir  = sldir
        self.slname = slname
        self.clean_starlog  = clean_starlog
        if not os.path.exists(self.sldir+'/'+self.slname):
	    if not os.path.exists(self.sldir+'/'+'star.log'):
		print 'error: no history.data file found in '+sldir
		print 'error: no star.log file found in '+sldir
	    else:
	    	self.slname='star.log'
		self.read_starlog()
	else:
            self.read_starlog()

    def __del__(self):
        print 'Closing', self.slname,' tool ...'

# let's start with functions that aquire data
    def read_starlog(self):
        ''' read history.data or star.log file again'''

        sldir   = self.sldir
        slname  = self.slname
        slaname = slname+'sa'
        if self.clean_starlog and os.path.exists(sldir+'/'+slaname):
            jonesmod=str(raw_input("Clean "+self.slname+" to make new "+self.slname+"sa? (y/n)"))
            if jonesmod == 'y':
                os.remove(sldir+'/'+slaname)
            
        if not os.path.exists(sldir+'/'+slaname):
            print 'No '+self.slname+'sa file found, create new one from '+self.slname
            cleanstarlog(sldir+'/'+slname)
        else:
            print 'Using old '+self.slname+'sa file ...'
            
        cmd=os.popen('wc '+sldir+'/'+slaname)    
        cmd_out=cmd.readline()
        cnum_cycles=cmd_out.split()[0]
        num_cycles=int(cnum_cycles) - 6

        filename=sldir+'/'+slaname

        header_attr,cols,data = read_mesafile(filename,data_rows=num_cycles)

        self.cols        = cols
        self.header_attr = header_attr
        self.data        = data
        
    def get(self,str_name):
        ''' return a column of data with the name str_name
        
        str_name is the name of the column as printed in history.data or
	star.log get the available columns from self.cols (where you replace
        self with the name of your instance'''

        column_array = self.data[:,self.cols[str_name]-1].astype('float')
        return column_array
        
    def hrd(self):
		''' plot an HR diagram '''
	
		pyl.plot(self.data[:,self.cols['log_Teff']-1],self.data[:,self.cols['log_L']-1],label = "M="+str(self.header_attr['initial_mass'])+", Z="+str(self.header_attr['initial_z']))
		pyl.legend()
		pyl.xlabel('log Teff')
		pyl.ylabel('log L')
    
    def hrd_key(self,key_str):
		''' plot an HR diagram 
		
		key_str    a label string'''
	
		pyl.plot(self.data[:,self.cols['log_Teff']-1],self.data[:,self.cols['log_L']-1],label = key_str)
		pyl.legend()
		pyl.xlabel('log Teff')
		pyl.ylabel('log L')
    
    def kippenhahn_CO(self,num_frame,xax,t0_model=0,title='Kippenhahn diagram',\
                       tp_agb=0.):
		''' Kippenhahn plot as a function of time or model with CO ratio
		
		num_frame    number of frame to plot this plot into 
                xax          string that is either model or time to
                             indicate what is to be used on the x-axis

                t0_model     model for the zero point in time, for AGB
                             plots this would be usually the model of
                             the 1st TP, which can be found with the
                             Kippenhahn plot 
                title        figure title

                tp_agb       if >= 0 then 
                             ylim=[h1_min*1.-tp_agb/100 : h1_max*1.+tp_agb/100] 
                             with h1_min, h1_max the min and max H-free 
                             core mass coordinate
                '''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction.'+\
			  ' needs to be "time" or "model"'
		
                t0_mod=xaxisarray[t0_model]
	    
		h1_boundary_mass  = self.get('he_core_mass')
		he4_boundary_mass = self.get('c_core_mass')
		star_mass         = self.get('star_mass')
		mx1_bot           = self.get('mx1_bot')*star_mass
		mx1_top           = self.get('mx1_top')*star_mass
		mx2_bot           = self.get('mx2_bot')*star_mass
		mx2_top           = self.get('mx2_top')*star_mass
		surface_c12       = self.get('surface_c12')
		surface_o16       = self.get('surface_o16')
	
		COratio=(surface_c12*4.)/(surface_o16*3.)
	
		pyl.plot(xaxisarray[t0_model:]-t0_mod,COratio[t0_model:],'-k',label='CO ratio')
		pyl.ylabel('C/O ratio')
		pyl.legend(loc=4)

		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
	
		pyl.twinx()
		pyl.plot(xaxisarray[t0_model:]-t0_mod,h1_boundary_mass[t0_model:],label='h1_boundary_mass')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,he4_boundary_mass[t0_model:],label='he4_boundary_mass')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_bot[t0_model:],',r',label='conv bound')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_top[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_bot[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_top[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,star_mass[t0_model:],label='star_mass')
		pyl.ylabel('mass coordinate')
		pyl.legend(loc=2)
                if tp_agb > 0.:
                    h1_min = min(h1_boundary_mass[t0_model:])
                    h1_max = max(h1_boundary_mass[t0_model:])
                    h1_min = h1_min*(1.-tp_agb/100.)
                    h1_max = h1_max*(1.+tp_agb/100.)
                    print 'setting ylim to zoom in on H-burning:',h1_min,h1_max                    
                    pyl.ylim(h1_min,h1_max)

    def kippenhahn(self,num_frame,xax,t0_model=0,title='Kippenhahn diagram',\
                       tp_agb=0.,t_eps=5.e2,plot_star_mass=True,symbol_size=8,\
                       c12_bm=False,print_legend=True):
		''' Kippenhahn plot as a function of time or model
		
		num_frame    number of frame to plot this plot into, if <0 open no new figure
                xax          string that is either 'model', 'time' or 'logtimerev' to
                             indicate what is to be used on the x-axis
                t_eps        final time for logtimerev             
                t0_model     model for the zero point in time, for AGB
                             plots this would be usually the model of
                             the 1st TP, which can be found with the
                             Kippenhahn plot 
                title        figure title

                tp_agb       if >= 0 then 
                             ylim=[h1_min*1.-tp_agb/100 : h1_max*1.+tp_agb/100] 
                             with h1_min, h1_max the min and max H-free 
                             core mass coordinate
                plot_star_mass    True - then plot the stellar mass as a line as well
                symbol_size  size of convection boundary marker
                c12_bm       boolean if we plot c12_boundary_mass or not
                print_legend       boolean, show or do not show legend
                '''
	
                if num_frame >= 0:
                    pyl.figure(num_frame)

		t0_mod=[]

		if xax == 'time':
		    xaxisarray = self.get('star_age')
                    t0_mod=xaxisarray[t0_model]
                    print 'zero time is '+str(t0_mod)
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
                    #t0_mod=xaxisarray[t0_model]
                    t0_mod = 0
		elif xax == 'logtimerev':
                    xaxi    = self.get('star_age')
                    xaxisarray = np.log10(np.max(xaxi)+t_eps-xaxi)
                    t0_mod = 0.
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction.'+\
			  ' needs to be "time" or "model"'
		
	    
		h1_boundary_mass  = self.get('he_core_mass')
		he4_boundary_mass = self.get('c_core_mass')
                if c12_bm:
                    c12_boundary_mass = self.get('c12_boundary_mass')
		star_mass         = self.get('star_mass')
		mx1_bot           = self.get('mx1_bot')*star_mass
		mx1_top           = self.get('mx1_top')*star_mass
		mx2_bot           = self.get('mx2_bot')*star_mass
		mx2_top           = self.get('mx2_top')*star_mass
	

		if xax == 'time':
                    if t0_model>0:
                        pyl.xlabel('$t - t_0$ $\mathrm{[yr]}$')
                    else:
                        pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
		elif xax == 'logtimerev':
                    pyl.xlabel('$\log(t_{final} - t)$  $\mathrm{[yr]}$')
	
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_bot[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size,label='convection zones')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_top[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size)
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_bot[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size)
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_top[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size)
		pyl.plot(xaxisarray[t0_model:]-t0_mod,h1_boundary_mass[t0_model:],color='red',linewidth=2,label='H-free core')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,he4_boundary_mass[t0_model:],color='green',linewidth=2,linestyle='dashed',label='He-free core')
                if c12_bm:
                    pyl.plot(xaxisarray[t0_model:]-t0_mod,c12_boundary_mass[t0_model:],color='purple',linewidth=2,linestyle='dotted',label='C-free core')
                if plot_star_mass is True:
                    pyl.plot(xaxisarray[t0_model:]-t0_mod,star_mass[t0_model:],label='$M_\star$')
		pyl.ylabel('$m_\mathrm{r}/\mathrm{M}_\odot$')
		if print_legend:
                    pyl.legend(loc=2)
                if tp_agb > 0.:
                    h1_min = min(h1_boundary_mass[t0_model:])
                    h1_max = max(h1_boundary_mass[t0_model:])
                    h1_min = h1_min*(1.-tp_agb/100.)
                    h1_max = h1_max*(1.+tp_agb/100.)
                    print 'setting ylim to zoom in on H-burning:',h1_min,h1_max                    
                    pyl.ylim(h1_min,h1_max)

    def t_surfabu(self,num_frame,xax,t0_model=0,title='surface abundance',t_eps=1.e-3,plot_CO_ratio=False):
		''' t_surfabu plots surface abundance evolution as a function of time
		
		num_frame    number of frame to plot this plot into, if <0 don't open figure
                xax          string that is either model, time or logrevtime
                             to indicate what is to be used on the x-axis

                t0_model     model for the zero point in time, for AGB
                             plots this would be usually the model of
                             the 1st TP, which can be found with the
                             Kippenhahn plot 
                title        figure title
                t_eps        time eps at end for logrevtime
                plot_CO_ratio onn second axis True/False 
                '''
                if num_frame >= 0:
                    pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')[t0_model:]
		elif xax == 'model':
		    xaxisarray = self.get('model_number')[t0_model:]
		elif xax == 'logrevtime':
		    xaxisarray = self.get('star_age')
                    xaxisarray=np.log10(max(xaxisarray[t0_model:])+t_eps-xaxisarray[t0_model:])
		else:
		    print 't-surfabu error: invalid string for x-axis selction.'+ \
			  ' needs to be "time" or "model"'
	    
		star_mass         = self.get('star_mass')
		surface_c12       = self.get('surface_c12')
		surface_c13       = self.get('surface_c13')
		surface_n14       = self.get('surface_n14')
		surface_o16       = self.get('surface_o16')                
                
                target_n14 = -3.5

	
		COratio=(surface_c12*4.)/(surface_o16*3.)
                t0_mod=xaxisarray[t0_model]
                log10_c12=np.log10(surface_c12[t0_model:])

                symbs=['k:','-','--','-.','b:','-','--','k-.',':','-','--','-.']

		pyl.plot(xaxisarray,log10_c12,\
                             symbs[0],label='$^{12}\mathrm{C}$')
		pyl.plot(xaxisarray,np.log10(surface_c13[t0_model:]),\
                             symbs[1],label='$^{13}\mathrm{C}$')
		pyl.plot(xaxisarray,np.log10(surface_n14[t0_model:]),\
                             symbs[2],label='$^{14}\mathrm{N}$')
		pyl.plot(xaxisarray,np.log10(surface_o16[t0_model:]),\
                             symbs[3],label='$^{16}\mathrm{O}$')
#                pyl.plot([min(xaxisarray[t0_model:]-t0_mod),max(xaxisarray[t0_model:]-t0_mod)],[target_n14,target_n14])

		pyl.ylabel('mass fraction $\log X$')
		pyl.legend(loc=2)

		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
                elif xax == 'logrevtime':
		    pyl.xlabel('$\\log t-tfinal$')
                if plot_CO_ratio:
                    pyl.twinx()
                    pyl.plot(xaxisarray,COratio[t0_model:],'-k',label='CO ratio')
                    pyl.ylabel('C/O ratio')
                    pyl.legend(loc=4)
                pyl.title(title)
                if xax == 'logrevtime':
                    self.xlimrev()


# ... end t_surfabu

    def t_lumi(self,num_frame,xax):
		''' Luminosity evolution as a function of time or model
		
		num_frame    number of frame to plot this plot into
		xax          string that is either model or time to indicate what is 
			     to be used on the x-axis'''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"'
		
	    
		logLH   = self.get('log_LH')
		logLHe  = self.get('log_LHe')
	
		pyl.plot(xaxisarray,logLH,label='L_(H)')
		pyl.plot(xaxisarray,logLHe,label='L(He)')
		pyl.ylabel('log L')
		pyl.legend(loc=2)
	
	
		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')

    def t_surf_parameter(self,num_frame,xax):
		''' Surface parameter evolution as a function of time or model
		
		num_frame    number of frame to plot this plot into
		xax          string that is either model or time to indicate what is 
			     to be used on the x-axis'''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"'
		
	    
		logL    = self.get('log_L')
		logTeff    = self.get('log_Teff')
	
		pyl.plot(xaxisarray,logL,'-k',label='log L')
		pyl.plot(xaxisarray,logTeff,'-k',label='log Teff')
		pyl.ylabel('log L, log Teff')
		pyl.legend(loc=2)
	
	
		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')

    def kip_vline(self,modstart,modstop,sparse,outfile,xlims=[0.,0.],ylims=[0.,0.],ixaxis='log_time_left',mix_zones=5,burn_zones=50):
        '''This function creates a Kippenhahn plot with energy flux using
        vertical lines, better thermal pulse resolution.
        For a more comprehensive plot, your history.data or star.log file should contain columns
        called "mix_type_n","mix_qtop_n","burn_type_n" and "burn_qtop_n".
        The number of columns (i.e. the bbiggest value of n) is what goes in the
        arguments as mix_zones and burn_zones.
        DO NOT WORRY! if you do not have these columns, just leave the default
        values alone and the script should recognise that you do not have these
        columns and make the most detailed plot that is available to you.

        modstart:        model from which you want to plot (be careful if your history.data
                         or star.log output is sparse...)
        modstop:         model to which you wish to plot
        xlims,ylims:     plot limits, however these are somewhat obsolete now that we
                         have modstart and modstop. Leaving them as 0. is probably
                         no slower, and you can always zoom in afterwards in mpl
        outfile:         'filename + extension' where you want to save the figure
        ixaxis:          either 'log_time_left', 'age', or 'model_number'
        sparse:          x-axis sparsity
        mix_zones,
        burn_zones:      As described above, if you have more detailed output about
                         your convection and energy generation boundaries in columns
                         mix_type_n, mix_qtop_n, burn_type_n and burn_qtop_n, you need
                         to specify the total number of columns for mixing zones and
                         burning zones that you have. Can't work this out from your
                         history.data or star.log file? Check the history_columns.list that you used, it'll
                         be the number after "mixing regions" and "burning regions".
                         Can't see these columns? leave it and 2 conv zones and 2 burn
                         zones will be drawn using other data that you certainly should
                         have in your history.data or star.log file.'''


        xxyy=[self.get('star_age')[modstart:modstop],self.get('star_age')[modstart:modstop]]
        mup = max(float(self.get('star_mass')[0])*1.02,1.0)
        nmodels=len(self.get('model_number')[modstart:modstop])

        Msol=1.98892E+33

        engenstyle = 'full'

        dx = sparse
        x = np.arange(0, nmodels, dx)

        btypemax = 20
        btypemin = -20
        btypealpha=0.

	########################################################################
	#----------------------------------plot--------------------------------#
	fig = pl.figure()
#	fig.set_size_inches(16,9)
	fsize=15
        ax=pl.axes()

	if ixaxis == 'log_time_left':
	# log of time left until core collapse
	    gage= self.get('star_age')
	    lage=np.zeros(len(gage))
	    agemin = max(abs(gage[-1]-gage[-2])/5.,1.e-10)
	    for i in np.arange(len(gage)):
	        if gage[-1]-gage[i]>agemin:
	            lage[i]=np.log10(gage[-1]-gage[i]+agemin)
	        else :
	            lage[i]=np.log10(agemin)
	    xxx = lage[modstart:modstop]
	    print 'plot versus time left'
	    ax.set_xlabel('$\mathrm{log}_{10}(t^*) \, \mathrm{(yr)}$',fontsize=fsize)
	elif ixaxis =='model_number':
	    xxx= self.get('model_number')[modstart:modstop]
	    print 'plot versus model number'
	    ax.set_xlabel('Model number',fontsize=fsize)
	elif ixaxis =='age':
	    xxx= self.get('star_age')[modstart:modstop]/1.e6
	    print 'plot versus age'
	    ax.set_xlabel('Age [Myr]',fontsize=fsize)
        else:
            print 'ixaxis must be one of: log_time_left, age or model_number'
            sys.exit()

        if xlims == [0.,0.]:
            xlims[0] = xxx[0]
            xlims[1] = xxx[-1]
        if ylims == [0.,0.]:
            ylims[0] = 0.
            ylims[1] = mup


	print 'plotting patches'
	ax.plot(xxx[::dx],self.get('star_mass')[modstart:modstop][::dx],'k-')

	print 'plotting abund boundaries'
	ax.plot(xxx,self.get('h1_boundary_mass')[modstart:modstop],label='H boundary')
	ax.plot(xxx,self.get('he4_boundary_mass')[modstart:modstop],label='He boundary')
#	ax.plot(xxx,self.get('c12_boundary_mass')[modstart:modstop],label='C boundary')

        ax.axis([xlims[0],xlims[1],ylims[0],ylims[1]])

        ax.set_ylabel('Mass [M$_\odot$]')

        ########################################################################

        try:
            self.get('burn_qtop_1')
        except:
            engenstyle = 'twozone'
        if engenstyle == 'full':
            for i in range(len(x)):
                # writing reading status 
	        percent = int(i*100/len(x))
	        sys.stdout.flush()
	        sys.stdout.write("\rcreating color map1 " + "...%d%%" % percent)
	        for j in range(1,burn_zones+1):
	            ulimit=self.get('burn_qtop_'+str(j))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	            if j==1:
	                llimit=0.0
	            else:
	                llimit=self.get('burn_qtop_'+str(j-1))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	            btype=float(self.get('burn_type_'+str(j))[modstart:modstop][i*dx])
	            if llimit!=ulimit:
                        if btype>0.:
	                    #btypealpha = btype/btypemax
                            #ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=btypealpha)
                            pass
		        if btype<0.:
                            #btypealpha = (btype/btypemin)/5
		            #ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='r',alpha=btypealpha)
                            pass

        if engenstyle == 'twozone':
            for i in range(len(x)):
            # writing reading status 
                  percent = int(i*100/len(x))
                  sys.stdout.flush()
                  sys.stdout.write("\rcreating color map1 " + "...%d%%" % percent)
                  llimitl1=self.get('epsnuc_M_1')[modstart:modstop][i*dx]/Msol
                  ulimitl1=self.get('epsnuc_M_4')[modstart:modstop][i*dx]/Msol
                  llimitl2=self.get('epsnuc_M_5')[modstart:modstop][i*dx]/Msol
                  ulimitl2=self.get('epsnuc_M_8')[modstart:modstop][i*dx]/Msol
                  llimith1=self.get('epsnuc_M_2')[modstart:modstop][i*dx]/Msol
                  ulimith1=self.get('epsnuc_M_3')[modstart:modstop][i*dx]/Msol
                  llimith2=self.get('epsnuc_M_6')[modstart:modstop][i*dx]/Msol
                  ulimith2=self.get('epsnuc_M_7')[modstart:modstop][i*dx]/Msol
                  # lower thresh first, then upper thresh:
                  #if llimitl1!=ulimitl1:
                      #ax.axvline(xxx[i*dx],ymin=(llimitl1-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimitl1-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=1.)
                  #if llimitl2!=ulimitl2:
                      #ax.axvline(xxx[i*dx],ymin=(llimitl2-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimitl2-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=1.)
                  #if llimith1!=ulimith1:
                      #ax.axvline(xxx[i*dx],ymin=(llimith1-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimith1-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=4.)
                  #if llimith2!=ulimith2:
                      #ax.axvline(xxx[i*dx],ymin=(llimith2-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimith2-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=4.)

        mixstyle = 'full'
        try:
            self.get('mix_qtop_1')
        except:
            mixstyle = 'twozone'
        if mixstyle == 'full':
	        for i in range(len(x)):
	        # writing reading status 
	          percent = int(i*100/len(x))
	          sys.stdout.flush()
	          sys.stdout.write("\rcreating color map2 " + "...%d%%" % percent)
	          for j in range(1,mix_zones+1):
	            ulimit=self.get('mix_qtop_'+str(j))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	            if j==1:
	              llimit=0.0
	            else:
	              llimit=self.get('mix_qtop_'+str(j-1))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	            mtype=self.get('mix_type_'+str(j))[modstart:modstop][i*dx]
	            if llimit!=ulimit:
		        if mtype == 1:
                            ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='k',alpha=3., linewidth=.5)
        if mixstyle == 'twozone':
	        for i in range(len(x)):
	        # writing reading status 
	          percent = int(i*100/len(x))
	          sys.stdout.flush()
	          sys.stdout.write("\rcreating color map2 " + "...%d%%" % percent)
                  ulimit=self.get('conv_mx1_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          llimit=self.get('conv_mx1_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          if llimit!=ulimit:
                      ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='k',alpha=5.,linewidth=.5)
                  ulimit=self.get('conv_mx2_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          llimit=self.get('conv_mx2_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          if llimit!=ulimit:
	              ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='k',alpha=3.,linewidth=.5)

        print 'engenstyle was ', engenstyle
        print 'mixstyle was ', mixstyle
	print '\n finished preparing color map'

        #fig.savefig(outfile)
        pl.show()

    def kip_cont(self,ifig=110,modstart=0,modstop=-1,outfile='out.png',xlims=[0.,0.],ylims=[0.,0.],xres=1000,yres=1000,ixaxis='model_number',mix_zones=20,burn_zones=20,plot_radius=False,engenPlus=True,engenMinus=False,landscape_plot=False,rad_lines=False,profiles=[],showfig=True,outlines=True):
        '''This function creates a Kippenhahn plot with energy flux using
        contours.

        This plot uses mixing_regions and burning_regions written to
        your history.data or star.log. Set both variables in the
        log_columns.list file to 20 as a start. 


        The output log file should then contain columns called
        "mix_type_n", "mix_qtop_n", "burn_type_n" and "burn_qtop_n".
        The number of columns (i.e. the bbiggest value of n) is what
        goes in the arguments as mix_zones and burn_zones.  DO NOT
        WORRY! if you do not have these columns, just leave the
        default values alone and the script should recognise that you
        do not have these columns and make the most detailed plot that
        is available to you.

        Defaults are set to get some plot, that may not look great if
        you zoom in interactively. Play with xres and yres as well as
        setting the xlims to ylims to the region you are interested
        in.


        ifig             figure frame number, default 110
        modstart:        model from which you want to plot (be careful if your history.data
                         or star.log output is sparse...), =0 (default))start from beginning,
                         works even if log_cnt > 1
        modstop:         model to which you wish to plot, default -1 corresponds to end
                         [if log_cnt>1, devide modstart and modstop by log_cnt, this needs 
                         to be improved!]
        xlims[DEPPREC.], plot limits, however these are somewhat obsolete now that we
        ylims            have modstart and modstop. Leaving them as 0. is probably
                         no slower, and you can always zoom in afterwards in mpl.
                         ylims is important for well resolved thermal pulse etc plots;
                         it's best to get the upper and lower limits of he-intershell using
                         s.kippenhahn_CO(1,'model') first.
        outfile:         'filename + extension' where you want to save the figure
        ixaxis:          either 'log_time_left', 'age', or 'model_number'
        xres,yres:       plot resolution. Needless to say that increasing these
                         values will yield a nicer plot with some slow-down in
                         plotting time. You will most commonly change xres. For a
                         prelim plot, try xres~200, then bump it up to anywhere from
                         1000-10000 for real nicely resolved, publication quality
                         plots.
        mix_zones,
        burn_zones:      As described above, if you have more detailed output about
                         your convection and energy generation boundaries in columns
                         mix_type_n, mix_qtop_n, burn_type_n and burn_qtop_n, you need
                         to specify the total number of columns for mixing zones and
                         burning zones that you have. Can't work this out from your
                         history.data or star.log file? Check the history_columns.list 
                         that you used, it'll
                         be the number after "mixing regions" and "burning regions".
                         Can't see these columns? leave it and 2 conv zones and 2 burn
                         zones will be drawn using other data that you certainly should
                         have in your history.data or star.log file.
        plot_radius      Whether on a second y-axis you want to plot the radius of the surface
                         and the he-free core.
        engenPlus        Boolean whether or not to plot energy generation contours for eps_nuc>0.
        endgenMinus      Boolean whether or not to plot energy generation contours for eos_nuc<0.'''

        xxyy=[self.get('star_age')[modstart:modstop],self.get('star_age')[modstart:modstop]]
        mup = max(float(self.get('star_mass')[0])*1.02,1.0)
        nmodels=len(self.get('model_number')[modstart:modstop])

        if ylims == [0.,0.]:
            mup   = max(float(self.get('star_mass')[0])*1.02,1.0)
            mDOWN = 0.
        else:
            mup = ylims[1]
            mDOWN = ylims[0]

        # y-axis resolution
        ny=yres
        #dy=mup/float(ny)
        dy = (mup-mDOWN)/float(ny)

        # x-axis resolution
        maxpoints=xres
        dx=int(max(1,nmodels/maxpoints))

        #y = np.arange(0., mup, dy)
        y = np.arange(mDOWN, mup, dy)
        x = np.arange(0, nmodels, dx)
        Msol=1.98892E+33

        engenstyle = 'full'

	B1=np.zeros([len(y),len(x)],float)
	B2=np.zeros([len(y),len(x)],float)
        try:
            self.get('burn_qtop_1')
        except:
            engenstyle = 'twozone'
        if engenstyle == 'full' and (engenPlus == True or engenMinus == True):
            ulimit_array = np.array([self.get('burn_qtop_'+str(j))[modstart:modstop:dx]*self.get('star_mass')[modstart:modstop:dx] for j in range(1,burn_zones+1)])
            #ulimit_array = np.around(ulimit_array,decimals=len(str(dy))-2)
            llimit_array = np.delete(ulimit_array,-1,0)
            llimit_array = np.insert(ulimit_array,0,0.,0)
            #llimit_array = np.around(llimit_array,decimals=len(str(dy))-2)
            btype_array = np.array([self.get('burn_type_'+str(j))[modstart:modstop:dx] for j in range(1,burn_zones+1)])
            for i in range(len(x)):
                percent = int(i*100/len(x))
	        sys.stdout.flush()
	        sys.stdout.write("\rcreating color map burn " + "...%d%%" % percent)
                for j in range(burn_zones):
                    if btype_array[j,i] > 0. and abs(btype_array[j,i]) < 99.:
                        B1[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 10.0**(btype_array[j,i])
                    elif btype_array[j,i] < 0. and abs(btype_array[j,i]) < 99.:
                        B2[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 10.0**(abs(btype_array[j,i]))

        if engenstyle == 'twozone' and (engenPlus == True or engenMinus == True):
                V=np.zeros([len(y),len(x)],float)
                for i in range(len(x)):
                # writing reading status 
                  percent = int(i*100/len(x))
                  sys.stdout.flush()
                  sys.stdout.write("\rcreating color map1 " + "...%d%%" % percent)
                  llimitl1=self.get('epsnuc_M_1')[modstart:modstop][i*dx]/Msol
                  ulimitl1=self.get('epsnuc_M_4')[modstart:modstop][i*dx]/Msol
                  llimitl2=self.get('epsnuc_M_5')[modstart:modstop][i*dx]/Msol
                  ulimitl2=self.get('epsnuc_M_8')[modstart:modstop][i*dx]/Msol
                  llimith1=self.get('epsnuc_M_2')[modstart:modstop][i*dx]/Msol
                  ulimith1=self.get('epsnuc_M_3')[modstart:modstop][i*dx]/Msol
                  llimith2=self.get('epsnuc_M_6')[modstart:modstop][i*dx]/Msol
                  ulimith2=self.get('epsnuc_M_7')[modstart:modstop][i*dx]/Msol
                  # lower thresh first, then upper thresh:
                  if llimitl1!=ulimitl1:
                      for k in range(ny):
                          if llimitl1<=y[k] and ulimitl1>y[k]:
                              V[k,i]=10.
                  if llimitl2!=ulimitl2:
                      for k in range(ny):
                          if llimitl2<=y[k] and ulimitl2>y[k]:
                              V[k,i]=10.
                  if llimith1!=ulimith1:
                      for k in range(ny):
                          if llimith1<=y[k] and ulimith1>y[k]:
                              V[k,i]=30.
                  if llimith2!=ulimith2:
                      for k in range(ny):
                          if llimith2<=y[k] and ulimith2>y[k]:
                              V[k,i]=30.

        mixstyle = 'full'
        try:
            self.get('mix_qtop_1')
        except:
            mixstyle = 'twozone'
        if mixstyle == 'full':
	    Z=np.zeros([len(y),len(x)],float)
            ulimit_array = np.array([self.get('mix_qtop_'+str(j))[modstart:modstop:dx]*self.get('star_mass')[modstart:modstop:dx] for j in range(1,mix_zones+1)])
            llimit_array = np.delete(ulimit_array,-1,0)
            llimit_array = np.insert(ulimit_array,0,0.,0)
            mtype_array = np.array([self.get('mix_type_'+str(j))[modstart:modstop:dx] for j in range(1,mix_zones+1)])
            for i in range(len(x)):
                percent = int(i*100/len(x))
	        sys.stdout.flush()
	        sys.stdout.write("\rcreating color map mix " + "...%d%%" % percent)
                for j in range(mix_zones):
                    if mtype_array[j,i] == 1.:
                        Z[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 1.

        if mixstyle == 'twozone':
	        Z=np.zeros([len(y),len(x)],float)
	        for i in range(len(x)):
	        # writing reading status 
	          percent = int(i*100/len(x))
	          sys.stdout.flush()
	          sys.stdout.write("\rcreating color map2 " + "...%d%%" % percent)
                  ulimit=self.get('conv_mx1_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          llimit=self.get('conv_mx1_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          if llimit!=ulimit:
	              for k in range(ny):
		          if llimit<=y[k] and ulimit>y[k]:
		              Z[k,i]=1.
                  ulimit=self.get('conv_mx2_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          llimit=self.get('conv_mx2_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
	          if llimit!=ulimit:
	              for k in range(ny):
		          if llimit<=y[k] and ulimit>y[k]:
		              Z[k,i]=1.

	if rad_lines == True:
		masses = np.arange(0.1,1.5,0.1)
		rads=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
		modno=[]
		for i in range(len(profiles)):
			p=mesa_profile('./LOGS',profiles[i])
			modno.append(p.header_attr['model_number'])
			for j in range(len(masses)):
				idx=np.abs(p.get('mass')-masses[j]).argmin()
				rads[j].append(p.get('radius')[idx])

        print 'engenstyle was ', engenstyle
        print 'mixstyle was ', mixstyle
	print '\n finished preparing color map'

	########################################################################
	#----------------------------------plot--------------------------------#
	fig = pyl.figure(ifig)
	fsize=20
        if landscape_plot == True:
		fig.set_size_inches(9,4)
	        pl.gcf().subplots_adjust(bottom=0.15)
	        pl.gcf().subplots_adjust(right=0.85)
	params = {'axes.labelsize':  fsize,
	  'text.fontsize':   fsize,
	  'legend.fontsize': fsize,
	  'xtick.labelsize': fsize*0.8,
	  'ytick.labelsize': fsize*0.8,
	  'text.usetex': False}
	pyl.rcParams.update(params)

        #ax=pl.axes([0.1,0.1,0.9,0.8])

        #fig=pl.figure()
        ax=pl.axes()

	if ixaxis == 'log_time_left':
	# log of time left until core collapse
	    gage= self.get('star_age')
	    lage=np.zeros(len(gage))
	    agemin = max(abs(gage[-1]-gage[-2])/5.,1.e-10)
	    for i in np.arange(len(gage)):
	        if gage[-1]-gage[i]>agemin:
	            lage[i]=np.log10(gage[-1]-gage[i]+agemin)
	        else :
	            lage[i]=np.log10(agemin)
	    xxx = lage[modstart:modstop]
	    print 'plot versus time left'
	    ax.set_xlabel('$\mathrm{log}_{10}(t^*) \, \mathrm{(yr)}$',fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [xxx[0],xxx[-1]]
	elif ixaxis =='model_number':
	    xxx= self.get('model_number')[modstart:modstop]
	    print 'plot versus model number'
	    ax.set_xlabel('Model number',fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [self.get('model_number')[modstart],self.get('model_number')[modstop]]
	elif ixaxis =='age':
            if modstart != 0:
		    xxx= self.get('star_age')[modstart:modstop] - self.get('star_age')[modstart]
		    print 'plot versus age'
		    ax.set_xlabel('Age [yr] + '+str(self.get('star_age')[modstart]),fontsize=fsize)
            else:
                    xxx= self.get('star_age')[modstart:modstop]/1.e6
		    ax.set_xlabel('Age [Myr]',fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [xxx[0],xxx[-1]]

        ax.set_ylabel('$\mathrm{Mass }(M_\odot)$')

#	cmapMIX = matplotlib.colors.ListedColormap(['w','k'])
        cmapMIX = matplotlib.colors.ListedColormap(['w','#8B8386']) # rose grey
#        cmapMIX = matplotlib.colors.ListedColormap(['w','#23238E']) # navy blue
        cmapB1  = pyl.cm.get_cmap('Blues')
	cmapB2  = pl.cm.get_cmap('Reds')
	#cmapB1  = matplotlib.colors.ListedColormap(['w','b'])
	#cmapB2  = matplotlib.colors.ListedColormap(['r','w'])
        if ylims == [0.,0.]:
            ylims[0] = 0.
            ylims[1] = mup
        if ylims[0] != 0.:
            ax.set_ylabel('$\mathrm{Mass }(M_\odot)$ + '+str(ylims[0]))
            y = y - ylims[0]
            ylims[0] = y[0]
            ylims[1] = y[-1]

        print xxx[::dx]
        print y
        print Z

	print 'plotting contours'
	CMIX    = ax.contourf(xxx[::dx],y,Z, cmap=cmapMIX,alpha=0.6,levels=[0.5,1.5])
        #CMIX_outlines    = ax.contour(xxx[::dx],y,Z, cmap=cmapMIX, alpha=1.0,levels=[0.5,1.5])
        #CMIX    = ax.contourf(xxx[::dx],y,Z, cmap=cmapMIX, alpha=0.5)
        if outlines == True:
	        CMIX_outlines    = ax.contour(xxx[::dx],y,Z, cmap=cmapMIX)

        if engenstyle == 'full' and engenPlus == True:
                print B1
                print B2
        	CBURN1  = ax.contourf(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.5, locator=matplotlib.ticker.LogLocator())
                CB1_outlines  = ax.contour(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.7, locator=matplotlib.ticker.LogLocator())
        	#CBURN1  = ax.contourf(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.5)
                CBARBURN1 = pyl.colorbar(CBURN1)
                CBARBURN1.set_label('$|\epsilon_\mathrm{nuc}-\epsilon_{\\nu}| \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$',fontsize=fsize)
        if engenstyle == 'full' and engenMinus == True:
        	CBURN2  = ax.contourf(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.5, locator=matplotlib.ticker.LogLocator())
                CBURN2_outlines  = ax.contour(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.7, locator=matplotlib.ticker.LogLocator())
#        	CBURN2  = ax.contourf(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.5)
        	CBARBURN2 = pl.colorbar(CBURN2)
                if engenPlus == False:
        	    CBARBURN2.set_label('$|\epsilon_\mathrm{nuc}-\epsilon_{\\nu}| \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$',fontsize=fsize)
        if engenstyle == 'twozone' and (engenPlus == True or engenMinus == True):
                print V
                ax.contourf(xxx[::dx],y,V, cmap=cmapB1, alpha=0.5)

	print 'plotting patches'
	ax.plot(xxx[::dx],self.get('star_mass')[modstart:modstop][::dx],'k-')

	print 'plotting abund boundaries'
	ax.plot(xxx,self.get('h1_boundary_mass')[modstart:modstop],label='H boundary',linestyle='-')
	ax.plot(xxx,self.get('he4_boundary_mass')[modstart:modstop],label='He boundary',linestyle='--')
	ax.plot(xxx,self.get('c12_boundary_mass')[modstart:modstop],label='C boundary',linestyle='-.')

        ax.axis([xlims[0],xlims[1],ylims[0],ylims[1]])

        if plot_radius == True:
            ax2=pyl.twinx()
            ax2.plot(xxx,np.log10(self.get('he4_boundary_radius')[modstart:modstop]),label='He boundary radius',color='k',linewidth=1.,linestyle='-.')
            ax2.plot(xxx,self.get('log_R')[modstart:modstop],label='radius',color='k',linewidth=1.,linestyle='-.')
            ax2.set_ylabel('log(radius)')
	if rad_lines == True:
	    ax2=pyl.twinx()
	    for i in range(len(masses)):
		    ax2.plot(modno,np.log10(rads[i]),color='k')

        fig.savefig(outfile,dpi=300)
	if showfig == True:
	        pyl.show()
# we may or may not need this below
#        fig.clear()


class star_log(history_data):
        '''Class derived from history_data class (copy). Existing just (for compatibility 
        reasons) for older mesa python scripts.'''



# below are some utilities that the user typically never calls directly


def read_mesafile(filename,data_rows=0,only='all'):
    ''' private routine that is not directly called by the user
    '''
    f=open(filename,'r')
    vv=[]
    v=[]
    lines = [] 
    line  = ''
    for i in range(0,6):
        line = f.readline()
        lines.extend([line])
    
    hval  = lines[2].split()
    hlist = lines[1].split()
    header_attr = {}
    for a,b in zip(hlist,hval):
        header_attr[a] = float(b)  
    if only is 'header_attr':
        return header_attr

    cols    = {}
    colnum  = lines[4].split()
    colname = lines[5].split()
    for a,b in zip(colname,colnum):
        cols[a] = int(b)
            
    data = []
    for i in range(data_rows):
        line = f.readline()
        v=line.split()
        try: 
            vv=np.array(v,dtype='float64')
        except ValueError:
            for item in v:
                if item.__contains__('.') and not item.__contains__('E'):
                    v[v.index(item)]='0'
        data.append(vv)
    f.close()
    a=np.array(data) 
    data = []
    return header_attr, cols, a


def cleanstarlog(file_in):
    ''' cleaning history.data or star.log file, e.g. to take care of repetitive restarts
    
    private, should not be called by user directly

    file_in     typically the filename of the mesa output history.data or star.log file,
                creates a clean file called history.datasa or star.logsa

    (thanks to Raphael for providing this tool)            
    '''

    file_out=file_in+'sa'
    f = open(file_in)
    lignes = f.readlines()
    f.close()

    nb    = np.array([],dtype=int)   # model number
    nb    = np.concatenate((nb    ,[  int(lignes[len(lignes)-1].split()[ 0])])) 
    nbremove = np.array([],dtype=int)   # model number
    i=-1

    for i in np.arange(len(lignes)-1,0,-1):
        line = lignes[i-1]
        
        if i > 6 and line != "" :
            if int(line.split()[ 0])>=nb[-1]:
                nbremove = np.concatenate((nbremove,[i-1])) 
            else:
                nb = np.concatenate((nb    ,[  int(line.split()[ 0])])) 
    i=-1
    for j in nbremove:
        lignes.remove(lignes[j])
 
    fout=file(file_out,'w')
    for j in np.arange(len(lignes)):
        fout.write(lignes[j])
    fout.close()
 
