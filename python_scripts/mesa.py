#!/usr/bin/env python
import matplotlib
matplotlib.use('Qt4Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as cm

import numpy as np
import numpy.lib.recfunctions
from astropy import units as u
from astropy import constants as const
from multiprocessing import Pool,cpu_count


from scipy.interpolate import interp1d



def LoadOneProfile(filename):
	"""Load one profile.

	Use numpy's genfromtxt to read input file. Top 5 lines are skipped.

	Args:
	filename (str) -- name of the file to be loaded

	Returns:
	data stored as a numpy array
	"""
	# This function needs to be outside the class, otherwise it is not picklable and it does nto work with pool.async

	data_from_file=np.genfromtxt(filename, skip_header=5, names=True)

	return data_from_file

def InterpolateOneProfile(profile, NY, Yaxis, Ymin, Ymax, Variable):
	"""Interpolate along a stellar profile.

	Args:
	profile (ndarray) -- the data structure holding a stellar profile
	NY (int) -- the number of data points to interpolate
	Yaxis (str) -- the independent variable along which to interpolate.
		Possible values: mass, radius, q, log_mass, log_radius, log_q
	Ymin (flt) -- the minimum independent variable value
	Ymax (flt) -- the maximum independent variable value
	Variable (str) -- the variable to be interpolated
		Possible values: eps_nuc, velocity, entropy, total_energy, j_rot,
			eps_recombination, ionization_energy, energy, potential_plus_kinetic,
			extra_heat, v_div_vesc, v_div_csound, pressure, temperature, density,
			tau, opacity, gamma1, dq, L_div_Ledd, Lrad_div_Ledd. t_thermal, t_dynamical,
			t_dynamical_down, t_thermal_div_t_dynamical, omega_div_omega_crit
			super_ad, vconv, vconv_div_vesc, conv_vel_div_csound, total_energy_plus_vconv2


	Returns:
	data interpolated along profile as a 1-D numpy array of length NY
	"""
	# This function needs to be outside the class, otherwise it is not picklable and it does nto work with pool.async


	# Add the fields log_mass, log_q and log_radius, in case they are not stored in the profile files

	if (not "log_mass" in profile.dtype.names):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'log_mass',
							data = np.log10(profile['mass']), asrecarray=True)
		except Exception:
			raise ValueError("Column 'mass' is missing from the profile files")

	if (not "log_q" in profile.dtype.names):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'log_q',
							data = np.log10(profile['q']), asrecarray=True)
		except Exception:
			raise ValueError("Column 'q' is missing from the profile files")

	if (not "log_radius" in profile.dtype.names):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'log_radius',
							data = np.log10(profile['radius']), asrecarray=True)
		except Exception:
			try:
				profile = numpy.lib.recfunctions.append_fields(profile,'log_radius',
								data = profile['logR'], asrecarray=True)
			except Exception:
				raise ValueError("Column 'logR' or 'log_radius' or 'radius' is missing from the profile files")

	if (not "radius" in profile.dtype.names):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'radius',
							data = 10.**profile['log_radius'], asrecarray=True)
		except Exception:
			try:
				profile = numpy.lib.recfunctions.append_fields(profile,'radius',
								data = 10.**profile['logR'], asrecarray=True)
			except Exception:
				raise ValueError("Column 'logR' or 'log_radius' or 'radius' is missing from the profile files")

	if (not "j_rot" in profile.dtype.names and (Variable == 'log_j_rot' or Variable == 'j_rot')):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'j_rot',
							data = 10.**profile['log_j_rot'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_j_rot' is missing from the profile files")


	if (not "potential_plus_kinetic" in profile.dtype.names and (Variable == 'potential_plus_kinetic' )):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'potential_plus_kinetic',
							data = profile['total_energy'] - profile['energy'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'total_energy' and/or 'energy' is missing from the profile files")

	if (not "v_div_vesc" in profile.dtype.names and (Variable == 'v_div_vesc' )):
		G = const.G.to('cm3/(g s2)').value
		Msun = const.M_sun.to('g').value
		Rsun = const.R_sun.to('cm').value
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'v_div_vesc',
							data = profile['velocity']/np.sqrt(2.*G*(profile['mass']*Msun)/(profile['radius']*Rsun)), asrecarray=True)
		except Exception:
			raise ValueError("Column 'radius' and/or 'velocity' is missing from the profile files")

	if Variable == 'temperature' :
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'temperature',
							data = 10.**profile['logT'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'logT' is missing from the profile files")

	if Variable == 'pressure' :
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'pressure',
							data = 10.**profile['logP'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'logP' is missing from the profile files")

	if Variable == 'density' :
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'density',
							data = 10.**profile['logRho'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'logRho' is missing from the profile files")

	if Variable == 'gamma1' :
		try:
			profile['gamma1'] = 4./3.-profile['gamma1']
		except Exception:
			raise ValueError("Column 'gamma1' is missing from the profile files")

	if Variable == 'tau' :
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'tau',
							data = 10.**profile['logtau'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'logtau' is missing from the profile files")

	if (not "opacity" in profile.dtype.names and (Variable == 'opacity' )):
		raise ValueError("Column 'opacity' is missing from the profile files")

	if (not "v_div_csound" in profile.dtype.names and (Variable == 'v_div_csound' )):
		raise ValueError("Column 'v_div_csound' is missing from the profile files")

	if (not "dq" in profile.dtype.names and (Variable == 'dq' )):
		raise ValueError("Column 'dq' is missing from the profile files")

	if (Variable == 'L_div_Ledd'):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'L_div_Ledd',
							data = 10.**profile['log_L_div_Ledd'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_L_div_Ledd' is missing from the profile files")

	if (Variable == 'Lrad_div_Ledd'):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'Lrad_div_Ledd',
							data = 10.**profile['log_Lrad_div_Ledd'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_Lrad_div_Ledd' is missing from the profile files")

	if (Variable == 't_thermal'):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'t_thermal',
							data = 10.**profile['log_thermal_time_to_surface'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_thermal_time_to_surface' is missing from the profile files")

	if (Variable == 't_dynamical'):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'t_dynamical',
							data = 10.**profile['log_acoustic_depth'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_acoustic_depth' is missing from the profile files")

	if (Variable == 't_dynamical_down'):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'t_dynamical_down',
							data = 10.**profile['log_acoustic_radius'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_acoustic_radius' is missing from the profile files")

	if (Variable == 't_thermal_div_t_dynamical'):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'t_thermal_div_t_dynamical',
							data = 10.**profile['log_thermal_time_to_surface']/10.**profile['log_acoustic_depth'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'log_thermal_time_to_surface' and/or 'log_acoustic_depth' are missing from the profile files")

	if (not "omega_div_omega_crit" in profile.dtype.names and (Variable == 'omega_div_omega_crit' )):
		raise ValueError("Column 'omega_div_omega_crit' is missing from the profile files")

	if (not "super_ad" in profile.dtype.names and (Variable == 'super_ad' )):
		raise ValueError("Column 'super_ad' is missing from the profile files")

	if (not "vconv" in profile.dtype.names and (Variable == 'vconv' )):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'vconv',
							data = profile['conv_vel'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'conv_vel' is missing from the profile files")

	if (not "vconv_div_vesc" in profile.dtype.names and (Variable == 'vconv_div_vesc' )):
		G = const.G.to('cm3/(g s2)').value
		Msun = const.M_sun.to('g').value
		Rsun = const.R_sun.to('cm').value
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'vconv_div_vesc',
							data = profile['conv_vel']/np.sqrt(2.*G*(profile['mass']*Msun)/(profile['radius']*Rsun)), asrecarray=True)
		except Exception:
			raise ValueError("Column 'radius' and/or 'conv_vel' is missing from the profile files")


	if (not "conv_vel_div_csound" in profile.dtype.names and (Variable == 'conv_vel_div_csound' )):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'conv_vel_div_csound',
							data = profile['conv_vel']/profile['csound'], asrecarray=True)
		except Exception:
			raise ValueError("Column 'csound' and/or 'conv_vel' is missing from the profile files")


	if (not "total_energy_plus_vconv2" in profile.dtype.names and (Variable == 'total_energy_plus_vconv2' )):
		try:
			profile = numpy.lib.recfunctions.append_fields(profile,'total_energy_plus_vconv2',
							data = profile['total_energy'] + 0.5*np.power(profile['conv_vel'],2), asrecarray=True)
		except Exception:
			raise ValueError("Column 'total_energy' and/or 'conv_vel' is missing from the profile files")






	#Define the values that we want to itnerpolate along the  Y axis
	Y_to_interp = (np.arange(1,NY+1).astype(float))/float(NY+2) * (Ymax-Ymin) + Ymin

	# Find the elements of Y_to_iterp that are valid for the specific profile data file.
	valid_points  = np.where((Y_to_interp < np.max(profile[Yaxis])) &
							(Y_to_interp > np.min(profile[Yaxis])))

	invalid_points  = np.where((Y_to_interp > np.max(profile[Yaxis])) &
							(Y_to_interp < np.min(profile[Yaxis])))


	interp_func = interp1d(profile[Yaxis], profile[Variable])

	data1 = np.zeros(NY)
	data1[:] = float('nan')
	data1[valid_points] = interp_func(Y_to_interp[valid_points])

	return data1


class mesa(object):
	def __init__(self, **kwargs):
		"""An object containing the output data from a MESA run.

		Upon initialization, class loads data, interpolates along axes, and
		stores interpolated data in instance of class.

		Args:
		data_path (str) -- file path of data directory (default "./")
		NX (int) -- number of data points on the x-axis (default 1024)
		NY (int) -- number of data points on the y-axis (default 1024)
		Yaxis (str) -- y-axis plotting variable (default 'mass')
		Xaxis (str) -- x-axis plotting variable (default 'star_age')
		Variable (str) -- z-axis (color) plotting variable (default 'eps_nuc')
		cmap (str) -- color scheme (default 'coolwarm')
		cmap_dynamic_range (flt) -- range in decades for cmap (default 10)
		Xaxis_dynamic_range (flt) -- x-axis range (default float('Inf'))
		Yaxis_dynamic_range (flt) -- y-axis range (default 4)
		figure_format (str) -- figure output format (default "eps")
		font_small (int) -- figure small font size (default 16)
		font_large (int) -- figure large font size (default 20)
		file_out (str) -- output figure name (default 'figure')
		onscreen (bool) -- plot figure on screen (default False)
		parallel (bool) -- create object using multiple processors (default True)
		abundances (bool) -- include abundance data in object (default False)
		log_abundances (bool) -- abudance data stored as logs (default True)
		czones (bool) - include convective zone data (default False)
		signed_log_cmap (bool) - use absolute value for color map (default True)
		orbit (bool) -- include secondary's position, for profile plots (default False)
		tau10 (bool) -- include optical depth of 10 in plot (default True)
		tau100 (bool) -- include optical depth of 100 in plot (default False)
		Nprofiles_to_plot (int) -- number of profiles to plot (default 10)
		profiles_to_plot (list[int]) -- profile numbers to be plotted (default [])
		"""

		self._param = {'data_path':"./", 'NX':1024, 'NY':1024, 'Yaxis':'mass', 'Xaxis':'star_age',
					'Variable':'eps_nuc', 'cmap':'coolwarm', 'cmap_dynamic_range':10, 'Xaxis_dynamic_range':float('Inf'),
					'Yaxis_dynamic_range':4, 'figure_format':"eps", 'font_small':16, 'font_large':20, 'file_out':'figure',
					'onscreen':False, 'parallel':True, 'abundances':False, 'log_abundances':True, 'czones':False,
					'signed_log_cmap':True, 'orbit':False, 'tau10':True, 'tau100':False, 'Nprofiles_to_plot':10,
					'profiles_to_plot':[]}

		for key in kwargs:
			if (key in self._param):
				self._param[key] = kwargs[key]
			else:
				raise ValueError(key+" is not a valid parameter name")

		self.CheckParameters()
		self.LoadData()
		self.InterpolateData()


	@property
	def data_path(self):
	    return self._param['data_path']

	@property
	def NY(self):
	    return self._param['NY']

	@property
	def NX(self):
	    return self._param['NX']

	@property
	def Yaxis(self):
	    return self._param['Yaxis']

	@property
	def Xaxis(self):
	    return self._param['Xaxis']

	@property
	def Variable(self):
	    return self._param['Variable']

	@property
	def cmap(self):
	    return self._param['cmap']

	@property
	def cmap_dynamic_range(self):
	    return self._param['cmap_dynamic_range']

	@property
	def Xaxis_dynamic_range(self):
	    return self._param['Xaxis_dynamic_range']

	@property
	def Yaxis_dynamic_range(self):
	    return self._param['Yaxis_dynamic_range']

	@property
	def figure_format(self):
	    return self._param['figure_format']

	@property
	def font_small(self):
	    return self._param['font_small']

	@property
	def font_large(self):
	    return self._param['font_large']

	@property
	def file_out(self):
	    return self._param['file_out']

	@property
	def onscreen(self):
	    return self._param['onscreen']

	@property
	def parallel(self):
	    return self._param['parallel']

	@property
	def abundances(self):
	    return self._param['abundances']

	@property
	def abundances(self):
	    return self._param['log_abundances']

	@property
	def czones(self):
	    return self._param['czones']

	@property
	def signed_log_cmap(self):
	    return self._param['signed_log_cmap']

	@property
	def orbit(self):
	    return self._param['orbit']

	@property
	def tau10(self):
	    return self._param['tau10']

	@property
	def signed_log_cmap(self):
	    return self._param['tau100']

	@property
	def Nprofiles_to_plot(self):
	    return self._param['Nprofiles_to_plot']

	@property
	def profiles_to_plot(self):
	    return self._param['profiles_to_plot']





	def help(self):
	#TODO: add a list of all parameters, the default values and the possible option, add a list of functions, and an example
		pass





	def CheckParameters(self):
		"""Check parameters to make sure valid options have been assigned.

		Possible values:
		Xaxis -- model_number, star_age, inv_star_age, log_model_number,
			log_star_age, log_inv_star_age
		Yaxis -- mass, radius, q, log_mass, log_radius, log_q
		Variable -- eps_nuc, velocity, entropy, total_energy, j_rot,
			eps_recombination, ionization_energy, energy, potential_plus_kinetic,
			extra_heat, v_div_vesc, v_div_csound, pressure, temperature, density,
			tau, opacity, gamma1, dq,L_div_Ledd, Lrad_div_Ledd, t_thermal, t_dynamical,
			t_dynamical_down, t_thermal_div_t_dynamical, omega_div_omega_crit

		cmap -- colors allowed by colormap module in matplotlib
		"""
		cmaps=[m for m in cm.datad]
		if not (self._param['cmap'] in cmaps):
			raise ValueError(self._param['cmap']+"not a valid option for parameter cmap")
		if not (self._param['Yaxis'] in ['mass', 'radius', 'q', 'log_mass', 'log_radius', 'log_q']):
			raise ValueError(self._param['Yaxis']+"not a valid option for parameter Yaxis")
		if not (self._param['Xaxis'] in ['model_number', 'star_age', 'inv_star_age', 'log_model_number', 'log_star_age',
				'log_inv_star_age']):
			raise ValueError(self._param['Xaxis']+"not a valid option for parameter Xaxis")
		if not (self._param['Variable'] in ['eps_nuc', 'velocity', 'entropy', 'total_energy', 'j_rot', 'eps_recombination'
				, 'ionization_energy', 'energy', 'potential_plus_kinetic', 'extra_heat', 'v_div_vesc',
				'v_div_csound',	'pressure', 'temperature', 'density', 'tau', 'opacity', 'gamma1', 'dq',
				'L_div_Ledd', 'Lrad_div_Ledd', 't_thermal', 't_dynamical', 't_dynamical_down', 't_thermal_div_t_dynamical',
 				'omega_div_omega_crit','super_ad', 'vconv', 'vconv_div_vesc', 'conv_vel_div_csound', 'total_energy_plus_vconv2']):
			raise ValueError(self._param['Variable']+"not a valid option for parameter Variable")


		return

	def SetParameters(self,**kwargs):
		"""Change the value of a parameter held by instance of mesa class"""
		for key in kwargs:
			if (key in self._param):
				self._param[key] = kwargs[key]
			else:
				raise ValueError(key+" is not a valid parameter name")

		#Check if any of the parameters that changed require reloading and reinterrpolating the data
		for key in kwargs:
			if key in ['data_path','NX','NY','Xaxis','Yaxis','Variable','Xaxis_dynamic_range','Yaxis_dynamic_range']:
				self.InterpolateData()
				break


		self.CheckParameters()

		return



	def LoadData(self):
		"""Load data from mesa outputs.

		Data is loaded from both history.data and individual
		profile*.data files. 'profiles.index' is used.
		"""
		# Read history with numpy so that we keep the column names and then convert then convert to a record array
		self.history = np.genfromtxt(self._param['data_path']+"history.data", skip_header=5, names=True)

		# Read list of available profile files
		self._profile_index = np.genfromtxt(self._param['data_path']+"profiles.index",skip_header=1,usecols=(0,2),
								dtype=[('model_number',int),('file_index',int)])
		if not len(self._profile_index["file_index"]):
			raise(self._param['data_path']+"profiles.index"+" does not contain information about enough profiles files")


		#Create an array fromt eh history data that gives you the age of each model for wich you have output a profile
		model_age_from_history =  interp1d(self.history["model_number"], self.history["star_age"])
		self.profile_age = model_age_from_history(self._profile_index["model_number"])

		#Check that the age is increasing in consecutive profiles. If not, then MESA ma have done a back up in which
		#case we should remove these profiles
		idx_profiles_to_keep = np.where(self.profile_age[:-1] < self.profile_age[1:])[0]
		self._profile_index = self._profile_index[:][idx_profiles_to_keep]
		self.profile_age = self.profile_age[idx_profiles_to_keep]

		#Check that the age is increasing in lines of history.data. If not, then MESA ma have done a back up in which
		#case we should remove these profiles
		idx_history_lines_to_keep = np.where(self.history['star_age'][:-1] < self.history['star_age'][1:])[0]
		self.history = self.history[:][idx_history_lines_to_keep]



		self.Nprofile = len(self._profile_index["file_index"])


		self.profiles=[]
		#Load the profile files and interpolate along the Y axis
		if self._param['parallel'] :
		    # Creates jobserver with ncpus workers
		    pool = Pool(processes=cpu_count())
		    print "Process running in parallel on ", cpu_count(), " cores"
		    filenames = [self._param['data_path']+"profile"+str(self._profile_index["file_index"][i])+".data" for i in range(self.Nprofile)]
		    results = [pool.apply_async(LoadOneProfile, args = (filename,)) for filename in filenames]
		    Nresults=len(results)
		    for i in range(0,Nresults):
		        self.profiles.append(results[i].get())
			pool.close()

		else:
			print "Process running serially"
			for i in range(self.Nprofile):
				filename = self._param['data_path']+"profile"+str(self._profile_index["file_index"][i])+".data"
				self.profiles.append(LoadOneProfile(filename))




	def InterpolateData(self):
		"""Data is interpolated

		Interpolations are performed in 1D, across
		the defined x-axis.
		"""
		# Set the maximum and the minimum of the X axis
		if self._param['Xaxis'] == "model_number":
			self._Xmax = np.max(self._profile_index["model_number"])
			self._Xmin = np.min(self._profile_index["model_number"])
		elif self._param['Xaxis'] == "star_age":
			self._Xmax = np.max(self.profile_age)
			self._Xmin = np.min(self.profile_age)
		elif self._param['Xaxis'] == "inv_star_age":
			self._Xmax = np.min(self.profile_age[-1] - self.profile_age)
			self._Xmin = np.max(self.profile_age[-1] - self.profile_age)
		elif self._param['Xaxis'] == "log_model_number":
			self._Xmax = np.max(np.log10(self._profile_index["model_number"]))
			self._Xmin = max([np.min(np.log10(self._profile_index["model_number"])),
							self._Xmax-self._param['Xaxis_dynamic_range']])
		elif self._param['Xaxis'] == "log_star_age":
			self._Xmax = np.max(np.log10(self.profile_age))
			self._Xmin = max([np.min(np.log10(self.profile_age)), self._Xmax-self._param['Xaxis_dynamic_range']])
		elif self._param['Xaxis'] == "log_inv_star_age":
			self._Xmin = np.max(np.log10(2.*self.profile_age[-1] - self.profile_age[-2] - self.profile_age))
			self._Xmax = max([np.min(np.log10(2.*self.profile_age[-1] - self.profile_age[-2] - self.profile_age)),
							self._Xmin-self._param['Xaxis_dynamic_range']])
		else:
			raise(self._param['Xaxis']+" is not a valid option for Xaxis")


		# Set the maximum and the minimum of the Y axis
		if self._param['Yaxis'] == "mass":
			self._Ymax = np.max(self.history["star_mass"])
			self._Ymin = 0.
		elif self._param['Yaxis'] == "radius":
			self._Ymax = 10.**np.max(self.history["log_R"])
			self._Ymin = 0.
		elif self._param['Yaxis'] == "q":
			self._Ymax = 1.
			self._Ymin = 0.
		elif self._param['Yaxis'] == "log_mass":
			self._Ymax = np.max(np.log10(self.history["star_mass"]))
			self._Ymin = np.max(np.log10(self.history["star_mass"])) - self._param['Yaxis_dynamic_range']
		elif self._param['Yaxis'] == "log_radius":
			self._Ymax = np.max(self.history["log_R"])
			self._Ymin = np.max(self.history["log_R"]) - self._param['Yaxis_dynamic_range']
		elif self._param['Yaxis'] == "log_q":
			self._Ymax = 0.
			self._Ymin = -self._param['Yaxis_dynamic_range']
		else:
			raise(self._param['Yaxis']+" is not a valid option for Yaxis")


		#Define the values that we want to itnerpolate along the the X and Y axis
		X_to_interp = (np.arange(1,self._param['NX']+1).astype(float))/float(self._param['NX']+2) * (self._Xmax-self._Xmin) + self._Xmin

		self._data = np.zeros((self._param['NX'],self._param['NY']))
		data_all = np.zeros((self.Nprofile,self._param['NY']))
		data_all[:,:] = float('nan')



		#Load the profile files and interpolate along the Y axis
		if self._param['parallel'] :
		    # Creates jobserver with ncpus workers
		    pool = Pool(processes=cpu_count())
		    print "Process running in parallel on ", cpu_count(), " cores"
		    results = [pool.apply_async(InterpolateOneProfile, args = (profile, self._param['NY'], self._param['Yaxis'],
		    			self._Ymin, self._Ymax, self._param['Variable'],)) for profile in self.profiles]
		    Nresults=len(results)
		    for i in range(0,Nresults):
		        data_all[i,:] = results[i].get()
			pool.close()
		else:
			print "Process running serially"
			for i in range(self.Nprofile):
				data_all[i,:] = InterpolateOneProfile(self.profiles[i], self._param['NY'], self._param['Yaxis'], self._Ymin,
													self._Ymax, self._param['Variable'])






		for i in range(self._param['NY']):
			if self._param['Xaxis'] == "model_number":
				Xaxis_values = interp1d(self._profile_index["model_number"].astype(float), data_all[:,i])
			elif self._param['Xaxis'] == "star_age":
				Xaxis_values = interp1d(self.profile_age, data_all[:,i])
			elif self._param['Xaxis'] == "inv_star_age":
				Xaxis_values = interp1d(self.profile_age[-1]-self.profile_age, data_all[:,i])
			elif self._param['Xaxis'] == "log_model_number":
				Xaxis_values = interp1d(np.log10(self._profile_index["model_number"].astype(float)), data_all[:,i])
			elif self._param['Xaxis'] == "log_star_age":
				Xaxis_values = interp1d(np.log10(self.profile_age), data_all[:,i])
			elif self._param['Xaxis'] == "log_inv_star_age":
				Xaxis_values = interp1d(np.log10(2.*self.profile_age[-1]-self.profile_age[-2]-self.profile_age), data_all[:,i])
			self._data[:,i] = Xaxis_values(X_to_interp)



	def Kippenhahn(self):
		"""Generate a Kippenhahn diagram

		Resulting plot is created and outputed.
		"""
		######################################################################
		# Set the labels of the two axis
		if self._param['Yaxis'] == "mass":
			Ylabel = "Mass Coordinate [$M_{\odot}$]"
		elif self._param['Yaxis'] == "radius":
			Ylabel = "Radius Coordinate [$R_{\odot}$]"
		elif self._param['Yaxis'] == "q":
			Ylabel = "Dimentionless Mass Coordinate q"
		elif self._param['Yaxis'] == "log_mass":
			Ylabel = "log(Mass Coordinate [$M_{\odot}$])"
		elif self._param['Yaxis'] == "log_radius":
			Ylabel = "log(Radius Coordinate [$R_{\odot}$])"
		elif self._param['Yaxis'] == "log_q":
			Ylabel = "log(Dimentionless Mass Coordinate q)"

		if self._param['Xaxis'] == "model_number":
			Xlabel = "Model Number"
		elif self._param['Xaxis'] == "star_age":
			Xlabel = "Star Age [yr]"
		elif self._param['Xaxis'] == "inv_star_age":
			Xlabel = "Time since the end of evolution [yr]"
		elif self._param['Xaxis'] == "log_model_number":
			Xlabel = "log(Model Number)"
		elif self._param['Xaxis'] == "log_star_age":
			Xlabel = "log(Star Age [yr])"
		elif self._param['Xaxis'] == "log_inv_star_age":
			Xlabel = "log(Time since the end of evolution [yr])"

		if self._param['Variable'] == "eps_nuc":
			cmap_label = "log($\epsilon_{nuclear}$ [erg/s/gr])"
		elif self._param['Variable'] == "velocity":
			cmap_label = "log(radial velocity [cm/s])"
		elif self._param['Variable'] == "entropy":
			cmap_label = "log(specific entropy [erg/K/gr])"
		elif self._param['Variable'] == "total_energy":
			cmap_label = "log(specific total energy [erg/gr])"
		elif self._param['Variable'] == "j_rot":
			cmap_label = "log(specific angular momentum [cm$^2$/s])"
		elif self._param['Variable'] == "eps_recombination":
			cmap_label = "log($\epsilon_{recombination}$ [erg/s/gr])"
		elif self._param['Variable'] == "ionization_energy":
			cmap_label = "log(specific ionization energy [erg/gr])"
		elif self._param['Variable'] == "dq":
			cmap_label = "Cell Mass Fraction"
		elif self._param['Variable'] == "energy":
			cmap_label = "log(specific internal energy [erg/gr])"
		elif self._param['Variable'] == "potential_plus_kinetic":
			cmap_label = "log(specific potential+kinetic energy [erg/gr])"
		elif self._param['Variable'] == "extra_heat":
			cmap_label = "log(specific injected energy rate [erg/s/gr])"
		elif self._param['Variable'] == "v_div_vesc":
			cmap_label = "log($v/v_{esc}$)"
		elif self._param['Variable'] == "v_div_csound":
			cmap_label = "log($v/v_{sound}$)"
		elif self._param['Variable'] == "gamma1":
			cmap_label = "log($4/3-\Gamma_{1}$)"
		elif self._param['Variable'] == "opacity":
			cmap_label = "log(Opacity [cm$^2$/gr])"
		elif self._param['Variable'] == "pressure":
			cmap_label = "log(Pressure [gr/(cm s$^2$)])"
		elif self._param['Variable'] == "density":
			cmap_label = "log(Density [cm$^2$/gr])"
		elif self._param['Variable'] == "temperature":
			cmap_label = "log(Temperature [K])"
		elif self._param['Variable'] == "tau":
			cmap_label = "log(optical depth)"
		elif self._param['Variable'] == "dq":
			cmap_label = "log(dq)"
		elif self._param['Variable'] == "L_div_Ledd":
			cmap_label = "log(L/L$_{Edd}$)"
		elif self._param['Variable'] == "Lrad_div_Ledd":
			cmap_label = "log(L$_{rad}$/L$_{Edd}$)"
		elif self._param['Variable'] == "t_thermal":
			cmap_label = "log($\\tau_{thermal}$ [s])"
		elif self._param['Variable'] == "t_dynamical":
			cmap_label = "log($\\tau_{s.cr.,out}$ [s])"
		elif self._param['Variable'] == "t_dynamical_down":
			cmap_label = "log($\\tau_{s.cr.,in}$ [s])"
		elif self._param['Variable'] == "t_thermal_div_t_dynamical":
			cmap_label = "log($\\tau_{thermal}/\\tau_{s.cr.}$)"
		elif self._param['Variable'] == "omega_div_omega_crit":
			cmap_label = "log($\Omega/\Omega_{crit}$)"

		elif self._param['Variable'] == "super_ad":
			cmap_label = "log($\\nabla - \\nabla_{ad}$)"
		elif self._param['Variable'] == "vconv":
			cmap_label = "log($v_{conv}$ [cm/s])"
		elif self._param['Variable'] == "vconv_div_vesc":
			cmap_label = "log($v_{conv}/v_{esc}$)"
		elif self._param['Variable'] == "conv_vel_div_csound":
			cmap_label = "log($v_{conv}/v_{sound}$)"
		elif self._param['Variable'] == "total_energy_plus_vconv2":
			cmap_label = "log(specific total energy + $1/2v_{conv}^2$ [erg/gr])"






		if self._param['signed_log_cmap']:
			cmap_label = "sign x log(max(1,abs(" + cmap_label[4:]+"))"


		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)
		fig1.subplots_adjust(top=0.80, left=0.12, right=0.9, bottom=0.12)
		ax1.set_xlabel(Xlabel,fontsize=self._param['font_large'])
		ax1.set_ylabel(Ylabel,fontsize=self._param['font_large'])
		ax1.xaxis.set_tick_params(labelsize = self._param['font_small'])
		ax1.yaxis.set_tick_params(labelsize = self._param['font_small'])
		ax1.set_xlim([self._Xmin,self._Xmax])
		ax1.set_ylim([self._Ymin,self._Ymax])


		if self._param['signed_log_cmap']:
			data_to_plot = np.sign(np.transpose(self._data)) * np.log10(np.maximum(1.,np.abs(np.transpose(self._data))))
		else:
			data_to_plot = np.log10(np.transpose(self._data))

		# When using signed_log_cmap, ignore cmap_dynamic_range
		if self._param['signed_log_cmap']:
			self._param['cmap_dynamic_range'] = np.nanmax(data_to_plot) - np.nanmin(data_to_plot)
		Image1 = plt.imshow(data_to_plot, origin='lower', cmap=self._param['cmap'],
							extent=[self._Xmin, self._Xmax, self._Ymin,self._Ymax], vmax = np.nanmax(data_to_plot),
							vmin=np.nanmax(data_to_plot)-self._param['cmap_dynamic_range'])
		ax1.set_aspect('auto')

		cbaxes = fig1.add_axes([0.12, 0.86, 0.78, 0.05])
		cbar1 = fig1.colorbar(Image1, cax=cbaxes, orientation='horizontal')
		cbar1.ax.xaxis.set_ticks_position('top')
		cbar1.ax.xaxis.set_label_position('top')
		cbar1.ax.set_xlabel(cmap_label, fontsize=self._param['font_small'])
		cbar1.ax.xaxis.set_tick_params(labelsize = self._param['font_small'])


#				Xaxis_values = interp1d(self._profile_index["model_number"].astype(float), data_all[:,i])

		#Plotting the orbit of the companion star inside the common envelope
		if self._param['orbit']:
			if self._param['Xaxis'] == "model_number":
				X_axis_orbit = self.history['model_number']
			elif self._param['Xaxis'] == "star_age":
				X_axis_orbit = self.history['star_age']
			elif self._param['Xaxis'] == "inv_star_age":
				X_axis_orbit = self.history['star_age'][-1] - self.history['star_age']
			elif self._param['Xaxis'] == "log_model_number":
				X_axis_orbit = np.log10(self.history['model_number'])
			elif self._param['Xaxis'] == "log_star_age":
				X_axis_orbit = np.log10(self.history['star_age'])
			elif self._param['Xaxis'] == "log_inv_star_age":
				X_axis_orbit = np.log10(2.*self.history['star_age'][-1]-self.history['star_age'][-2]-self.history['star_age'])



			if self._param['Yaxis'] == "mass":
				Y_axis_orbit = self.history['CE_companion_position_m']
			elif self._param['Yaxis'] == "radius":
				Y_axis_orbit = self.history['CE_companion_position_r']
			elif self._param['Yaxis'] == "q":
				Y_axis_orbit = self.history['CE_companion_position_m']/self.history['star_mass']
			elif self._param['Yaxis'] == "log_mass":
				Y_axis_orbit = np.log10(self.history['CE_companion_position_m'])
			elif self._param['Yaxis'] == "log_radius":
				Y_axis_orbit = np.log10(self.history['CE_companion_position_r'])
			elif self._param['Yaxis'] == "log_q":
				Y_axis_orbit = np.log10(self.history['CE_companion_position_m']/self.history['star_mass'])

			ax1.plot(X_axis_orbit, Y_axis_orbit, linewidth=3, color='black')




		#Plotting the tau=10 surface
		if self._param['tau10']:
			if self._param['Xaxis'] == "model_number":
				X_axis_tau10 = self.history['model_number']
			elif self._param['Xaxis'] == "star_age":
				X_axis_tau10 = self.history['star_age']
			elif self._param['Xaxis'] == "inv_star_age":
				X_axis_tau10 = self.history['star_age'][-1] - self.history['star_age']
			elif self._param['Xaxis'] == "log_model_number":
				X_axis_tau10 = np.log10(self.history['model_number'])
			elif self._param['Xaxis'] == "log_star_age":
				X_axis_tau10 = np.log10(self.history['star_age'])
			elif self._param['Xaxis'] == "log_inv_star_age":
				X_axis_tau10 = np.log10(2.*self.history['star_age'][-1]-self.history['star_age'][-2]-self.history['star_age'])



			if self._param['Yaxis'] == "mass":
				Y_axis_tau10 = self.history['tau10_mass']
			elif self._param['Yaxis'] == "radius":
				Y_axis_tau10 = self.history['tau10_radius']
			elif self._param['Yaxis'] == "q":
				Y_axis_tau10 = self.history['tau10_mass']/self.history['star_mass']
			elif self._param['Yaxis'] == "log_mass":
				Y_axis_tau10 = np.log10(self.history['tau10_mass'])
			elif self._param['Yaxis'] == "log_radius":
				Y_axis_tau10 = np.log10(self.history['tau10_radius'])
			elif self._param['Yaxis'] == "log_q":
				Y_axis_tau10 = np.log10(self.history['tau10_mass']/self.history['star_mass'])

			ax1.plot(X_axis_tau10, Y_axis_tau10, "--",linewidth=2, color='lightgray')


		#Plotting the tau=100 surface
		if self._param['tau100']:
			if self._param['Xaxis'] == "model_number":
				X_axis_tau100 = self.history['model_number']
			elif self._param['Xaxis'] == "star_age":
				X_axis_tau100 = self.history['star_age']
			elif self._param['Xaxis'] == "inv_star_age":
				X_axis_tau100 = self.history['star_age'][-1] - self.history['star_age']
			elif self._param['Xaxis'] == "log_model_number":
				X_axis_tau100 = np.log10(self.history['model_number'])
			elif self._param['Xaxis'] == "log_star_age":
				X_axis_tau100 = np.log10(self.history['star_age'])
			elif self._param['Xaxis'] == "log_inv_star_age":
				X_axis_tau100 = np.log10(2.*self.history['star_age'][-1]- self.history['star_age'][-2] -self.history['star_age'])



			if self._param['Yaxis'] == "mass":
				Y_axis_tau100 = self.history['tau100_mass']
			elif self._param['Yaxis'] == "radius":
				Y_axis_tau100 = self.history['tau100_radius']
			elif self._param['Yaxis'] == "q":
				Y_axis_tau100 = self.history['tau100_mass']/self.history['star_mass']
			elif self._param['Yaxis'] == "log_mass":
				Y_axis_tau100 = np.log10(self.history['tau100_mass'])
			elif self._param['Yaxis'] == "log_radius":
				Y_axis_tau100 = np.log10(self.history['tau100_radius'])
			elif self._param['Yaxis'] == "log_q":
				Y_axis_tau100 = np.log10(self.history['tau100_mass']/self.history['star_mass'])

			ax1.plot(X_axis_tau100, Y_axis_tau100, "--",linewidth=2, color='darkgray')




		#Plot convecitve zones
#TODO: Add check that convective zones info are included in the history file
		if self._param['czones']:
			if self._param['Xaxis'] == "model_number":
				X_axis_czones = self.history['model_number']
			elif self._param['Xaxis'] == "star_age":
				X_axis_czones = self.history['star_age']
			elif self._param['Xaxis'] == "inv_star_age":
				X_axis_czones = self.history['star_age'][-1] - self.history['star_age']
			elif self._param['Xaxis'] == "log_model_number":
				X_axis_czones = np.log10(self.history['model_number'])
			elif self._param['Xaxis'] == "log_star_age":
				X_axis_czones = np.log10(self.history['star_age'])
			elif self._param['Xaxis'] == "log_inv_star_age":
				X_axis_czones = np.log10(2.*self.history['star_age'][-1]-self.history['star_age'][-2]-self.history['star_age'])


			if self._param['Yaxis'] == "mass":
				conv_mx1_top = interp1d(X_axis_czones, self.history['conv_mx1_top']*self.history['star_mass'])
				conv_mx1_bot = interp1d(X_axis_czones, self.history['conv_mx1_bot']*self.history['star_mass'])
				conv_mx2_top = interp1d(X_axis_czones, self.history['conv_mx2_top']*self.history['star_mass'])
				conv_mx2_bot = interp1d(X_axis_czones, self.history['conv_mx2_bot']*self.history['star_mass'])
			elif self._param['Yaxis'] == "radius":
				conv_mx1_top = interp1d(X_axis_czones, self.history['conv_mx1_top_r'])
				conv_mx1_bot = interp1d(X_axis_czones, self.history['conv_mx1_bot_r'])
				conv_mx2_top = interp1d(X_axis_czones, self.history['conv_mx2_top_r'])
				conv_mx2_bot = interp1d(X_axis_czones, self.history['conv_mx2_bot_r'])
			elif self._param['Yaxis'] == "q":
				conv_mx1_top = interp1d(X_axis_czones, self.history['conv_mx1_top'])
				conv_mx1_bot = interp1d(X_axis_czones, self.history['conv_mx1_bot'])
				conv_mx2_top = interp1d(X_axis_czones, self.history['conv_mx2_top'])
				conv_mx2_bot = interp1d(X_axis_czones, self.history['conv_mx2_bot'])
			elif self._param['Yaxis'] == "log_mass":
				conv_mx1_top = interp1d(X_axis_czones, np.log10(self.history['conv_mx1_top']*self.history['star_mass']))
				conv_mx1_bot = interp1d(X_axis_czones, np.log10(self.history['conv_mx1_bot']*self.history['star_mass']))
				conv_mx2_top = interp1d(X_axis_czones, np.log10(self.history['conv_mx2_top']*self.history['star_mass']))
				conv_mx2_bot = interp1d(X_axis_czones, np.log10(self.history['conv_mx2_bot']*self.history['star_mass']))
			elif self._param['Yaxis'] == "log_radius":
				conv_mx1_top = interp1d(X_axis_czones, np.log10(self.history['conv_mx1_top_r']))
				conv_mx1_bot = interp1d(X_axis_czones, np.log10(self.history['conv_mx1_bot_r']))
				conv_mx2_top = interp1d(X_axis_czones, np.log10(self.history['conv_mx2_top_r']))
				conv_mx2_bot = interp1d(X_axis_czones, np.log10(self.history['conv_mx2_bot_r']))
			elif self._param['Yaxis'] == "log_q":
				conv_mx1_top = interp1d(X_axis_czones, np.log10(self.history['conv_mx1_top']))
				conv_mx1_bot = interp1d(X_axis_czones, np.log10(self.history['conv_mx1_bot']))
				conv_mx2_top = interp1d(X_axis_czones, np.log10(self.history['conv_mx2_top']))
				conv_mx2_bot = interp1d(X_axis_czones, np.log10(self.history['conv_mx2_bot']))



			N_cz_lines = 200
			X_cz = np.arange(self._Xmin, self._Xmax, (self._Xmax-self._Xmin)/N_cz_lines)
			cz1_top = conv_mx1_top(X_cz)
			cz1_bot = conv_mx1_bot(X_cz)
			cz2_top = conv_mx2_top(X_cz)
			cz2_bot = conv_mx2_bot(X_cz)
			for i in range(N_cz_lines):
				ax1.plot([X_cz[i], X_cz[i]], [cz1_bot[i],cz1_top[i]], color='grey', alpha=0.5)
				ax1.plot([X_cz[i], X_cz[i]], [cz2_bot[i],cz2_top[i]], color='grey', alpha=0.5)




		#Plot surface and central abundances
		if self._param['abundances']:
			ax2 = ax1.twinx()
			if self._param['log_abundances']:
				ax2.set_ylim([1e-5,1.])
				ax2.set_yscale('log')
			else:
				ax2.set_ylim([0.,1.])
			ax2.set_ylabel('Mass Fraction', color='grey',fontsize=self._param['font_small'])
			for tl in ax2.get_yticklabels():
				tl.set_color('grey')
			ax2.set_aspect('auto')
			colors = ['black', 'blue', 'red', 'green', 'cyan', 'magenta', 'yellow', 'DarkBlue', 'Orange',
					'Brown', 'DarkCyan', 'DarkMagenta', 'DarkYellow']

			if self._param['Xaxis'] == "model_number":
				X_axis_abundances = self.history['model_number']
			elif self._param['Xaxis'] == "star_age":
				X_axis_abundances = self.history['star_age']
			elif self._param['Xaxis'] == "inv_star_age":
				X_axis_abundances = self.history['star_age'][-1]-self.history['star_age']
			elif self._param['Xaxis'] == "log_model_number":
				X_axis_abundances = np.log10(self.history['model_number'])
			elif self._param['Xaxis'] == "log_star_age":
				X_axis_abundances = np.log10(self.history['star_age'])
			elif self._param['Xaxis'] == "log_inv_star_age":
				X_axis_abundances = np.log10(2.*self.history['star_age'][-1]-self.history['star_age'][-1]-self.history['star_age'])



			if 'center_h1' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_h1'], ':', linewidth=2, color=colors[1])
				fig1.text(0.12, 0.81, "$H^1$", fontsize = self._param['font_small'], color=colors[1])
			if 'center_he4' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_he4'], ':', linewidth=2, color=colors[2])
				fig1.text(0.18, 0.81, "$He^4$", fontsize = self._param['font_small'], color=colors[2])
			if 'center_c12' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_c12'], ':', linewidth=2, color=colors[3])
				fig1.text(0.24, 0.81, "$C^{12}$", fontsize = self._param['font_small'], color=colors[3])
			if 'center_n14' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_n14'], ':', linewidth=2, color=colors[4])
				fig1.text(0.30, 0.81, "$N^{14}$", fontsize = self._param['font_small'], color=colors[4])
			if 'center_o16' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_o16'], ':', linewidth=2, color=colors[5])
				fig1.text(0.36, 0.81, "$O^{16}$", fontsize = self._param['font_small'], color=colors[5])
			if 'center_ne20' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_ne20'], ':', linewidth=2, color=colors[6])
				fig1.text(0.42, 0.81, "$Ne^{20}$", fontsize = self._param['font_small'], color=colors[6])
			if 'center_mg24' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_mg24'], ':', linewidth=2, color=colors[7])
				fig1.text(0.48, 0.81, "$Mg^{24}$", fontsize = self._param['font_small'], color=colors[7])
			if 'center_si28' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_si28'], ':', linewidth=2, color=colors[8])
				fig1.text(0.54, 0.81, "$Si^{28}$", fontsize = self._param['font_small'], color=colors[8])
			if 'center_fe56' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['center_fe56'], ':', linewidth=2, color=colors[9])
				fig1.text(0.60, 0.81, "$Fe^{56}$", fontsize = self._param['font_small'], color=colors[9])


			if 'surface_h1' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_h1'], '--', linewidth=1, color=colors[1])
				fig1.text(0.12, 0.81, "$H^1$", fontsize = self._param['font_small'], color=colors[1])
			if 'surface_he4' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_he4'], '--', linewidth=1, color=colors[2])
				fig1.text(0.18, 0.81, "$He^4$", fontsize = self._param['font_small'], color=colors[2])
			if 'surface_c12' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_c12'], '--', linewidth=1, color=colors[3])
				fig1.text(0.24, 0.81, "$C^{12}$", fontsize = self._param['font_small'], color=colors[3])
			if 'surface_n14' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_n14'], '--', linewidth=1, color=colors[4])
				fig1.text(0.30, 0.81, "$N^{14}$", fontsize = self._param['font_small'], color=colors[4])
			if 'surface_o16' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_o16'], '--', linewidth=1, color=colors[5])
				fig1.text(0.36, 0.81, "$O^{16}$", fontsize = self._param['font_small'], color=colors[5])
			if 'surface_ne20' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_ne20'], '--', linewidth=1, color=colors[6])
				fig1.text(0.42, 0.81, "$Ne^{20}$", fontsize = self._param['font_small'], color=colors[6])
			if 'surface_mg24' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_mg24'], '--', linewidth=1, color=colors[7])
				fig1.text(0.48, 0.81, "$Mg^{24}$", fontsize = self._param['font_small'], color=colors[7])
			if 'surface_si28' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_si28'], '--', linewidth=1, color=colors[8])
				fig1.text(0.54, 0.81, "$Si^{28}$", fontsize = self._param['font_small'], color=colors[8])
			if 'surface_fe56' in self.history.dtype.names:
				ax2.plot(X_axis_abundances, self.history['surface_fe56'], '--', linewidth=1, color=colors[9])
				fig1.text(0.60, 0.81, "$Fe^{56}$", fontsize = self._param['font_small'], color=colors[9])



		#Re-enforcing the calculated limits for X and Y axis. Without this there may be a white band on the right of the plot.
		ax1.set_xlim([self._Xmin,self._Xmax])
		ax1.set_ylim([self._Ymin,self._Ymax])
		if self._param['abundances']:
			ax2.set_xlim([self._Xmin,self._Xmax])
			if self._param['log_abundances']:
				ax2.set_ylim([1e-5,1.])
			else:
				ax2.set_ylim([0.,1.])


#		fig1.tight_layout()
		fig1.savefig(self._param['file_out']+"."+self._param['figure_format'], format=self._param['figure_format'])


		# fig2 = plt.figure()
		# ax2 = fig2.add_subplot(111)
		# fig2.subplots_adjust(top=0.99, left=0.12, right=0.99, bottom=0.12)
		# ax2.set_xlabel(Ylabel,fontsize=self._param['font_small'])
		# ax2.set_ylabel(cmap_label,fontsize=self._param['font_small'])
		# ax2.xaxis.set_tick_params(labelsize = self._param['font_small'])
		# ax2.yaxis.set_tick_params(labelsize = self._param['font_small'])
		# ax2.set_xlim([self._Ymin,self._Ymax])
		# ax2.set_ylim([np.nanmax(data_to_plot)-self._param['cmap_dynamic_range'],np.nanmax(data_to_plot)])
		#
		# ax2.plot(10.**self.profiles[10]['logR'], np.log10(self.profiles[10]['extra_heat']), linewidth=3, color='black')



		if self._param['onscreen']:
			plt.show()
			# fig1.canvas.manager.window.raise_()

		return


	def Demo(self):
#TODO: add function that has a script that test most of the functions of the class. This should serve as a self test for the developers
# to make sure that we are not brake between revisions
		pass
		return







if __name__ == "__main__":


#TODO Add linear cmap




	#Options for Xaxis: 'model_number', 'star_age', 'inv_star_age', 'log_model_number', 'log_star_age', 'log_inv_star_age'
	#Options for Yaxis: 'mass', 'radius', 'q', 'log_mass', 'log_radius', 'log_q'
	#Options for Variable: "eps_nuc", "velocity", "entropy", "total_energy", "j_rot", "eps_recombination", "ionization_energy",
	#						"energy", "potential_plus_kinetic", "extra_heat", "v_div_vesc", "v_div_csound"
	#						"pressure", "temperature", "density", "tau", "opacity", "gamma1", "dq"
	#						"super_ad", "vconv", "vconv_div_vesc", "conv_vel_div_csound", "total_energy_plus_vconv2"


	data_path = "/home/evol/fragkos/disk1/mesa_projects/CE/example_runs/heat_base_lim_vconv/LOGS/"
	a = mesa(data_path=data_path, parallel=True, abundances=False, log_abundances = True, Yaxis='mass', Xaxis="star_age",
		czones=True, Variable='vconv_div_vesc', orbit=True)
	a.SetParameters(onscreen=True, cmap = 'jet', cmap_dynamic_range=5, signed_log_cmap=False)

	a.Kippenhahn()
