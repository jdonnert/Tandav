
module TandavCodeParameters

export CodeParameters

type CodeParameters
	
	gamma	:: Float64	# adiabatic index
	xH		:: Float64	# Hydrogen Fraction
	boxsize	:: Float64	# Boxsize
	h		:: Float64	# hubble parameter
	omega0	:: Float64
	omegaL	:: Float64
	fRad	:: Bool		# Flag Radiative run
	fDouble	:: Bool		# Flag DOUBLEPRECISION
	fSfr	:: Bool		# Flag Star Formation
	fCool	:: Bool		# Flag Cooling
	fFeedB	:: Bool		# Flag Feedback
	fComov	:: Bool		# Flag COMOVING
	fPeriod	:: Bool		# Flag PERIODIC
	
	function CodeParameters(;gamma=5/3, xH=0.76, boxsize=0, h=1, omega0=1,
						 	omegaL=0.7, fRad=false, fDouble=false, fSfr=false,
						 	fCool=false, fFeedB=false, fComov=false, fPeriod=false)
	
		println("Setting Code Parameters:\n"*
		  		"    gamma     = $gamma\n"*
				"    xH        = $xH\n"*
				"    boxsize   = $boxsize\n"*
				"    omega0    = $omega0\n"*
				"    omegaL    = $omegaL\n"*
				"    fRad      = $fRad\n"* 
				"    fDouble   = $fDouble\n"*
				"    fSfr      = $fSfr\n"*   
				"    fCool     = $fCool\n"*  
				"    fFeedB    = $fFeedB\n"* 
				"    fComov    = $fComov\n"* 
				"    fPeriod   = $fPeriod\n")

		new(gamma, xH, boxsize, h, omega0, omegaL, fRad, fDouble, fSfr,
	  		fCool, fFeedB, fComov, fPeriod)
	end
end

end # module
