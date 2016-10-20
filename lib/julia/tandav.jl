# Tandav Code module.
# Define a super type "Tandav" that contains a full self-contained description 
# of a simulation

module Tandav

importall CGSUnits

importall TandavCodeUnits
importall TandavCodeParameters
importall TandavCodeSnapshotIO

export TandavCodeObject
export Density, NumberDensity, U2T, T2U, ThermalEnergyDensity, SoundSpeed 

type TandavCodeObject # super type, glues all submodules
	
	z		:: Float64 			# Current Redshift
	
	Par 	:: CodeParameters
	Unit	:: CodeUnits
	Const	:: CodeConstants

	function TandavCodeObject(z=0; 	
					gamma=5/3, xH=0.76, boxsize=0, h=1, omega0=1,
					omegaL=0.7, fRad=false, fDouble=false, fSfr=false,
					fCool=false, fFeedB=false, fComov=false, fPeriod=false,
					unitName="kpc, 1e10 Msol, km/s", length=kpc2cm,
					mass=1e10*Msol, vel=1e5) 
	
		println("\nSetting CodeObject with \n"*
		  		"    z=$z\n")

		Par = CodeParameters(;gamma=gamma, xH=xH, boxsize=boxsize, h=h, 
					   		omega0=omega0, omegaL=omegaL, fRad=fRad, 
					   		fDouble=fDouble, fSfr=fSfr, fCool=fCool, 
							fFeedB=fFeedB, fComov=fComov, fPeriod=fPeriod)

		Unit = CodeUnits(unitName, length, mass, vel; show=true)
		
		Const = CodeConstants(Unit, show=true)

		new(z, Par, Unit, Const)
	end

end

# here start the function wrappers for all subtyp methods

function Density(t::TandavCodeObject, rho)

	return Density(t.Unit, rho; h=t.Par.h, z=t.z)
end

function NumberDensity(t::TandavCodeObject, rho)

	return NumberDensity(t.Unit, rho; h=t.Par.h, z=t.z, xH=t.Par.xH)
end

function Pressure(t::TandavCodeObject, rho, u)

	return Pressure(t.Unit, rho, u; h=t.Par.h, z=t.z, xH=t.Par.xH, gam=t.Par.Gamma)
end

function U2T(t::TandavCodeObject, u)

	return U2T(t.Unit, u; xH=t.Par.xH, gam=t.Par.gamma, rad=t.Par.fRad )
end

function T2U(t::TandavCodeObject, temp)
  
	return U2T(t.Unit, temp; xH=t.Par.xH, gam=t.Par.gamma, rad=t.Par.fRad)
end

function ThermalEnergyDensity(t::TandavCodeObject, rho, u)

	return ThermalEnergyDensity(t.Unit, rho, u; h=t.Par.h, z=t.h, xH=t.Par.xH,
							 	gam=t.Par.gamma, radiative=t.Par.fRad)
end

function SoundSpeed(t::TandavCodeObject, u)

	return SoundSpeed(t.Unit, u; gam=t.Par.gamma)
end

end # module
