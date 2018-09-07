class load_constants(object):
    """
    # Load a table of constants for thermodynamic calculations
    # Includes physical constants, default microphysical parameters
    #
    # Optional input argument chooses the set of constants to use
    #
    # 'default' : Exact thermodynamics under the assumption that isobaric specific heat capacities are constant
    # 'bolton'  : Use Bolton's formulas for saturation vapor pressure
    # 'teten'   : Use Teten's formulas for saturation vapor pressure
    # 'sam'     : Use constants consistent with SAM cloud-resolving model
    # 'fms'     : Use constants consistent with Simplified GCM based on the GFDL FMS model 
    #		(as used by O'Gorman & Schneider etc.)
    """
    def __init__(self,varargin='default',**kwargs):
        if type(varargin) == str:
            modeltype = varargin
        elif type(varargin) == list:
            modeltype = varargin[0]
        else:
            raise ValueError('Wrong input type! Please input either a str or list')

        self.modeltype = modeltype    
        ## Default microphysical parameters
        #self.gamma = 0 # Default precipitation fallout parameter (0 is reversible, 1 is pseudo-adiabatic)
        self.ice = 1              # include ice?
        self.deltaT = 40		# mixed-phase range (K)
        
        ## Physical self.nts (based on constants in CM1)
        # Isobaric heat capacities (J/kg/K)
        self.cp        = 1005.7    	# dry air
        self.cpv       = 1870.0		# water vapor
        self.cpl       = 4190.0		# liquid water
        self.cpi       = 2106.0		# solid water

        # Gas constants (J/kg/K)
        self.Rd        = 287.04		# dry air
        self.Rv        = 461.5

        # Latent heats (J/kg)
        self.Lv0       = 2501000.0	# liquid-gas
        self.Ls0       = 2834000.0	# solid-gas

        # gravitational acceleration (m/s^2)
        self.g         = 9.81

        # Liquid water density (kg/m^3)
        self.rhol = 1000

        ## Reference values
        self.T0        = 273.16 		# temperature (K)
        self.p00 	    = 100000		# Pressure (Pa)
        self.e0        = 611.2		# vapor pressure (Pa)

        ## Non thermodynamic constants
        self.sigma = 5.67e-8

        ## Altered constants for different thermodynamics
        if modeltype == 'default':
            next
        elif modeltype == 'bolton':
            next
        elif modeltype == 'teten':
            next
        elif modeltype == 'fms':
            # choose the following for consistency with the FMS code
            self.Rd                  = 287.04                 # J / kg / K
            self.eps                 = 0.622                  # gas self.nt water vapor [J/kg/K]
            self.Rv                  = self.Rd/self.eps
            self.kappa             = 2.0/7.0
            self.cp                = self.Rd/self.kappa     # specific heat at constant pressure
            self.cpv               = self.cp                 # specific heat water vapor [J/kg/K] - SAME AS FOR DRY AIR
            self.cpl               = self.cp                 # specific heat liquid water [J/kg/K] - SAME AS FOR DRY AIR
            self.cp_ocean          = 3989.24495292815         # heat capacity ocean water [J/kg/K]
            self.rhol              = 1000                     # density of liquid water [kg/m^3]
            self.Lv0               = 2.5e6                    # latent heat of evaporation [J/kg]
            self.ice         	  = 0
            self.g                 = 9.80665 

        elif modeltype == 'sam':
            self.cp = 1004.0      #                 ! Specific heat of air, J/kg/K
            self.Lv0 = 2.5104e+06 #                 ! Latent heat of condensation, J/kg
            self.Ls0 = 2.8440e+06 #                 ! Latent heat of sublimation, J/kg
            self.Rv = 461.        #                 ! Gas constant for water vapor, J/kg/K
            self.Rd = 287.        #                 ! Gas constant for dry air, J/kg/K
            self.deltaT = 20
            self.tbgmin = 253.16  #                 ! Minimum temperature for cloud water., K
            self.tbgmax = 273.16  #                 ! Maximum temperature for cloud ice, K
            self.tprmin = 268.16  #                 ! Minimum temperature for rain, K
            self.tprmax = 283.16  #                 ! Maximum temperature for snow+graupel, K
            self.tgrmin = 223.16  #                 ! Minimum temperature for snow, K
            self.tgrmax = 283.16  #                 ! Maximum temperature for graupel, K
            self.deltaT = 20			    # Consistent with tbgmin/tbgmax above
        else:
            raise ValueError('Cannot find the constants requested')

    #     ## Override constants if input
    #     if nargin >= 2 const.ice  = varargin{2} end 
    #     if nargin >= 3 const.deltaT  = varargin{3} end

        if type(varargin) == list:
            if len(varargin)>1:
                self.ice = varargin[1]
            if len(varargin)>2:
                self.deltaT = varargin[2]

        ## Derived parameters
        self.cv        = self.cp-self.Rd
        self.cvv       = self.cpv-self.Rv
        self.eps       = self.Rd/self.Rv
            
