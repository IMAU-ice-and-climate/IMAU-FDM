Definition of each variable in FDM. 

Model settings:
- nyears = simulation time [yr]
- nyearsSU = simulation time during the spin up period [yr]
- dtmodelExp = time step in model with explicit T-scheme [s] 
- dtmodelImp = time step in model with implicit T-scheme [s]
- dtmodel = time step in model with chosen scheme
- ImpExp = Impicit or Explicit scheme (1=Implicit/fast, 2= Explicit/slow)
- dtobs = time step in input data [s]
- dtSnow = duration of running average of variables used in snow parameterisation [s]
- dzmax = vertical model resolution [m]
- initdepth = initial depth of firn profile [m]
- th = theta (if theta=0.5 , it is a Crank Nicolson scheme) 
- startasice = indicates the initial rho-profile (1=linear, 2=ice)
- beginT = indicates the inital T-profile (0=winter, 1=summer, 2=linear)
- NoR = number of times the data series is repeated for the initial rho profile is constructed
- writeinspeed = frequency of writing speed components to file
- writeinprof = frequency of writing of firn profiles to file
- proflayers = number of output layers in the prof file
- writeindetail = frequency of writing to detailed firn profiles to file
- detlayers = number of output layers in the detailed prof file
- detthick = thickness of output layer in the detailed prof file
- lon_current = Lon; indicates the longitude gridpoint
- lat_current = Lat; indicates the latitude gridpoint
- Nlon = Number of longitude points in forcing
- Nlat = Number of latitude points in forcing
- Nt_forcing = Number of timesteps in forcing

- Nt_model_interpol = number of time steps in FDM per forcing time step
- Nt_model_tot = total number of time steps in FDM
- Nt_model_spinup = total number of time steps in spin-up period

- numOutputProf, numOutputSpeed, numOutputDetail = number of time steps between writing 2D data, speed components and detailed data to the output file []
- outputProf, outputSpeed, outputDetail = number of times that output is written to the files

- spinup_bound = maximum amount of spinups
- error_bound = boundary for elevation and FAC changes in spinup, when changes are smaller than the error_bound the spinup is done
- spinup_numb = number of spinups
- z_surf_error, fac_error = changes in elevation and FAC per spinup


Constants:
- pi = asin(1.) * 2 = 3.1415
- R = 8.3145 = gas constand [J mole-1 K-1]
- g = 9.81 = gravitational acceleration [m s-2]
- rhoi = 917. = density of ice [kg-3]
- Ec = 60000 = activation energy for creep [J mole-1]
- Eg = 42400 = activation energy for grain growth [J mole-1]
- Lh = 333500 = latent heat for fusion [J kg-1]


Model input (pay attention to units of variables, can be different!):
- LSM = ice mask
- Longitude = longitude as described in ice mask
- Latitude = latitude as described in ice mask
- ind_lon, ind_lat = index of gridpoint
- lon_grid, lat_grid = longitude and latitude of grid point
- lonfile, latfile = 
- AveMelt = average snow melt [kg m-2 s-1]
- AveAcc = average liquid and solid precipitation [kg m-2 s-1]
- AveWind = average 10m wind speed [m/s]
- AveTsurf = average surface temperature [K]
- AveSubl = average sublimation [kg m-2 s-1]
- AveSnowDrif = average snow drift [kg m-2 s-1]
- SnowMelt = timeseries of snowmelt [kg m-2 s-1]
- PreTot = timeseries of liquid and solid precipitation [kg m-2 s-1]
- PreSol = timeseries of snowfall [kg m-2 s-1]
- Sublim = timeseries of sublimation [kg m-2 s-1]
- TempSurf = timeseries of surface temperature [K]
- SnowDrif = timeseries of snowdrift [kg m-2 s-1]
- FF10m = timeseries of wind speed [m/s]
- remove_Melt, remove_Psol, remove_Ptot = amount of a variable that is removed as it is below a certain threshold
- PreLiq = liquid precipitation based on input of total and solid precipitation [kg m-2 s-1]

- acav = avarage accumulation at the specified grid point []
- ffav = average 10 m wind speed at the specified grid point [m/s]
- tsav = average surface temperatire at the specified grid point [K]
- TempFM, PSolFM, PLiqFM, SublFM, MeltFM, DriftFM, ff10FM = timeseries of variable after interpolation (some divided by Nt_model as it is a mass flux isntead of a rate)
- Rho0FM = fresh snow density calculated with interpolated values
- rho0, Ts, Psol, Pliq, Su, Me, Sd = variable at specified time step after interpolation
- numSnow =
- TempSnow = temperature of fersh snow [K]

Layer information in model
- Rho = layer density []
- M = layer mass [kg]
- T = layer temperature [K]
- Depth = depth middle of the layer [m]
- DZ = layer thickness [m]
- DenRho =
- Refreeze = 
- Year =
- dzmax = vertical model resolution []

Column information in model
- h_surf = surface height[]
- FirnAir = total firn air content in column []
- TotLwc = total liquid water content in column []
- IceMass = total ice mass in column []
- DenRho =

- vfc = vertical movement column due to firn compaction []
- msolin = mass added to upper layer []
- vacc = vertical movement column due to snowfall []
- vsub = vertical movement column due to sublimation []
- vsnd = vertical movement column due to snowdrift []
- vmelt = vertical movement column due to snowmelt []
- Mmelt = mass of melt calculated based on input of melt, precipitation, and temperature []
- vice =

Variables related to firn physics
- kice = thermal conductivity of ice []
- kice_ref = thermal conductivity of ice at reference 270.15 K []
- kcal = thermal conductivity of snow []
- kf = thermal conductivity of firn []
- kair = thermal conductivity of air []
- kair_ref = thermal conductivity of air at reference 270.15 K []
- ki = thermal conductivity []
- Diff = thermal diffusivity []
- MO = correction factor for firn densification which differs for densities above and below 550 kg m-3
- krate = 
- ci = heat capacity []
- kip = thermal conductivity of current layer []
- kis = thermal conductivity of layer below []
- kint = thermal conductivity at the interface []
- Gn = incoming flux from upper layer []
- Gs = outgoing flux to the lower layer []


Variables related to water physics
- cp0 = heat capacity at 273.15 K []
- cp = heat capacity []
- Efreeze = energy available for freezing in layer []
- Lh = latent heat of fusion []
- Mfreeze = mass that can be frozen with freezing energy in layer []
- poro = fraction of volume occupied by air
- maxpore = amount of liquid water the firn or snow can retain per mass unit
- MavailCol = mass of liquid water that can be retained in the layer [kg]
- MavailMax = mass of liquid water that can be retained in the layer (same as MavailCol, but different method) [kg]
- Mavail = mass of available pore space in the layer [kg]
- Madd = mass added to the layer [kg]
- Mrefreeze = total mass of refreozen water in column [kg]
- Refreeze = mass of refrozen water in layer [kg]
- Mlwc = mass liquid water content [kg]
- Mporespaceavail = mass of pore space availability in a layer after adding liquid water

Output variables 

Not done yet