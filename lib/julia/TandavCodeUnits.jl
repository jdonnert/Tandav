# provide code unit definition and conversion

module TandavCodeUnits

importall CGSUnits

export CodeUnits, CodeConstants
export Density, NumberDensity, Pressure, U2T, T2U, ThermalEnergyDensity, 
	   SoundSpeed

type CodeUnits

	unitName::AbstractString
	length::Float64
	mass::Float64
	vel::Float64
	time::Float64
	energy::Float64

	function CodeUnits(unitName, length, mass, vel; show=false) # inner constructor
		
		time = length/vel
		energy = mass*vel^2

		if show == true
			println("Setting Units: $unitName")
			println("    length = $length cm")
			println("    mass   = $mass g")
			println("    vel    = $vel cm/s")
			println("    time   = $time s")
			println("    energy = $energy gcm²/s\n")
		end

		new(unitName, length, mass, vel, time, energy)
	end

end # CodeUnits

type CodeConstants
	
	G :: Float64	# Newtons constans

	function CodeConstants(unit::CodeUnits; show=false)
	
		G = grav / unit.length^3 * unit.mass * unit.time^2
		
		if show == true
			
			println("Constants in Code units: ")
			println("    G      = $G cm³/g/s²")
		end

		new(G)
	end

end

# conversions Code -> physical

function Density(unit::CodeUnits, rho; h=1, z=0)

	return float(rho) * (1+z)^3 * unit.mass/unit.length^3 * h^2
end

function NumberDensity(unit::CodeUnits, rho; h=1, z=0, xH=0.76)

	n2ne = (xH+0.5*(1-xH))/(2*xH+0.75*(1-xH))
	umu = 4./(5.*xH+3.)

	return Density(unit, rho; h=h, z=z) * n2ne/(umu*mp)
end

function Pressure(unit::CodeUnits, rho, U; h=1, z=0, xH=0.76, gam=5/3)

	rho_cgs = Density(unit, rho, h=h, z=z)

	return rho_cgs .* U * (gam-1) * unit.vel^2

end

function U2T(unit::CodeUnits, U; xH=0.76, gam=5/3, rad=false)

	yhelium = (1-xH)/(4*xH)

	mean_mol_weight = (1+4*yhelium)/(1+3*yhelium+1)

	return U * (gam-1) * unit.vel^2 * mp * mean_mol_weight / k_boltz
end

function T2U(unit::CodeUnits, T; xH=0.76, gam=5/3, rad=false)

	return T / U2T(unit, 1, xH=xH, gam=gam, rad=rad)
end

function ThermalEnergyDensity(unit::CodeUnits, rho, u; h=1, z=0, xH=0.76, gam=5/3, 
							  radiative=false)

	rho_cgs = Density(unit, rho, h=h, z=z)
	
	T = U2T(unit, u, xH-xH, gam=gam, radiative=radiative)

    umol = 4./(5.*xH+3.)
	
	return rho_cgs .* T *k_boltz/(mp*umol) 
end

function SoundSpeed(unit::CodeUnits, U; gam=5/3)

	return sqrt(U * gam * (gam-1))
end

end # module
