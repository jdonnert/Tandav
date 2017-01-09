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
export ReadHead, ReadSnap

type TandavCodeObject # super type, glues all submodules
	
	z		:: Float64 			# Current Redshift
	
	Par 	:: CodeParameters
	Unit	:: CodeUnits
	Const	:: CodeConstants

	# inner constructor
	function TandavCodeObject(z=0 ; 	
					gamma=5/3, xH=0.76, boxsize=0, h=1, omega0=1,
					omegaL=0.7, fDouble=false, fSfr=false,
					fCool=false, fFeedB=false, fComov=false, fPeriod=false,
					unitName="kpc, 1e10 Msol, km/s", length=kpc2cm,
					mass=1e10*Msol, vel=1e5) 
	
		println("\nSetting CodeObject with \n"*
		  		  "    z=$z\n")

		Par = CodeParameters(;gamma=gamma, xH=xH, boxsize=boxsize, h=h, 
					   		omega0=omega0, omegaL=omegaL,  
					   		fDouble=fDouble, fSfr=fSfr, fCool=fCool, 
							fFeedB=fFeedB, fComov=fComov, fPeriod=fPeriod)

		Unit = CodeUnits(unitName, length, mass, vel; show=true)
		
		Const = CodeConstants(Unit, show=true)

		new(z, Par, Unit, Const)
	end

end

# here start the function wrappers for all subtype methods

function Density(t::TandavCodeObject, rho::Array)

	return Density(t.Unit, rho; h=t.Par.h, z=t.z)
end


function NumberDensity(t::TandavCodeObject, rho::Array)

	return NumberDensity(t.Unit, rho; h=t.Par.h, z=t.z, xH=t.Par.xH)
end

function Pressure(t::TandavCodeObject, rho::Array, u::Array)

	return Pressure(t.Unit, rho, u; h=t.Par.h, z=t.z, xH=t.Par.xH, gam=t.Par.Gamma)
end

function U2T(t::TandavCodeObject, u::Array)

	return U2T(t.Unit, u; xH=t.Par.xH, gam=t.Par.gamma, rad=t.Par.fCool )
end

function T2U(t::TandavCodeObject, temp::Array)
  
	return T2U(t.Unit, temp; xH=t.Par.xH, gam=t.Par.gamma, rad=t.Par.fCool)
end

function ThermalEnergyDensity(t::TandavCodeObject, rho::Array, u::Array)

	return ThermalEnergyDensity(t.Unit, rho, u; h=t.Par.h, z=t.h, xH=t.Par.xH,
							 	gam=t.Par.gamma, radiative=t.Par.fCool)
end

function SoundSpeed(t::TandavCodeObject, u::Array)

	return SoundSpeed(t.Unit, u; gam=t.Par.gamma)
end


function ReadSnap(t::TandavCodeObject, fname::AbstractString, label::String; 
				  pType=0x07, debug=false)
	
	head = ReadHead(t, fname; debug=debug)

	data = ReadSnap(fname, label; pType=pType, debug=debug)

	return data
end

function ReadHead(t::TandavCodeObject, fname::AbstractString; debug=false)

	head = ReadHead(fname; debug=debug)

	if debug == true

		println("Updating TandavCodeObject from HEAD of Snapshot ")

		Update!(t.z, head.redshift, "z")	
		Update!(t.Par.boxsize, head.boxsize, "boxsize")
		Update!(t.Par.omega0, head.omega0, "omega0")
		Update!(t.Par.omegaL, head.omegaL, "omegaL")
		Update!(t.Par.h, head.hbpar, "h")

		Update!(t.Par.fDouble, head.flag_double, "fDouble")
		Update!(t.Par.fSfr, head.flag_sfr, "fSfr")
		Update!(t.Par.fCool, head.flag_cooling, "fCool")
		Update!(t.Par.fFeedB, head.flag_feedback, "fFeedB")
		Update!(t.Par.fComov, head.flag_comoving, "fComov")
		Update!(t.Par.fPeriod, head.flag_periodic, "fPeriod")
	
		println()
	end

	return head
end

function Update!(a, b, name::AbstractString)
	
	if a != b
		
		a = b

		println("   HEAD -> $name = $a")
	end
end

end # module
