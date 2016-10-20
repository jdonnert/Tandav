module Cosmology

importall CGSUnits

export Params
export Show, Hubble

type Params

	Model 		:: String	# Name of model
	H0			:: Float64		# Hubble constant at z=0
	Omega_b		:: Float64		# Baryon density parameter
	Omega_M		:: Float64  	# Matter density
	Omega_L		:: Float64		# Dark energy density parameter
	Omega_r		:: Float64		# Photon energer density parameter
	sigma_8		:: Float64		# Fluctuation amplitude at 8/h Mpc
	#Delta_r2	:: Float64		# Curvature fluctuation amplitude
	#n_s		:: Float64 	 	# scalar spectral index
	#tau		:: Float64 	 	# reionisation optical depth
	#w			:: Float64 	 	# equation of state parameter
	#z_star		:: Float64 	 	# redshift of decoupling
	#t_star		:: Float64 	 	# Age of decoupling
	#z_reion	:: Float64 	 	# Redshift of reionisation

	# Derived Parameters	

	Omega_tot	:: Float64  	# total density
	hbpar		:: Float64  	# h parameter
	H100		:: Float64  	# hubble constant with 100 km/s/Mpc
	#t0 			:: Float64		# Age of universe

	function Params(i=1; show=false) # inner constructor
	
		H100 = 100 * 1e5 / mpc2cm

		if i == 0
	
			model = "NONE"
			H0 = H100
			Omega_L = 1
			Omega_b = 1
			Omega_M = 1
			sigma_8 = 1
			Omega_r = 1

		elseif i == 1
	
			model = "Concordance"
			H0 = 70 * 1e5 / mpc2cm
			Omega_L = 0.7
			Omega_b = 0.03
			Omega_M = 0.3
			sigma_8 = 0.8
			Omega_r = 0

		else
			throw(TypeError())
		end
		
		hbpar = H0 / H100
		Omega_tot = Omega_M+Omega_L+Omega_r

		cosmo = new(model, H0, Omega_b, Omega_L, Omega_r, sigma_8, Omega_tot, hbpar, H100)

		if show == true
			Show(cosmo)
		end

		return cosmo

	end # constructor

end # type

function Show(cosmo::Params)

	println("$(cosmo.Model) Cosmology : ")
	println("   H0       = $(cosmo.H0/(1e5 / mpc2cm)) km/s/Mpc")
	println("   Omega_L  = $(cosmo.Omega_L)")
	println("   Omega_M  = $(cosmo.Omega_M)")
	println("   Omega_b  = $(cosmo.Omega_b)")
	println("   Omega_r  = $(cosmo.Omega_r)")
	println("   sigma_8  = $(cosmo.sigma_8)")
	#println("   Horizon  = $(d_hubble/(1e3*kpc2cm))")

	println( "\nMethods (all cgs) :")
 	println( "  scale factor at age t :           a(cosmo,t)")
	println( "  Hubble paramter:                  H(cosmo,z)")
	println( "  Comoving Distance:                d_comov(cosmo,z)")
	println( "  Transverse Comoving Distance :    d_transcomov(vz)")
	println( "  Angular Diameter Distance :       d_ang(cosmo,z)")
	println( "  Luminosity Distance :             d_lum(cosmo,z)")
	println( "  Critical Density :                rho_crit(cosmo,z)")
	println( "  Overdensity Parameter :           Delta(cosmo,z)")
	println( "  Redshift to time :                z2t(cosmo,z)")
	println( "                                    t2z(cosmo,t)")
	println( "  On the Sky :                      arcmin2kpc(cosmo,arcmin, z)")
	println( "                                    kpc2arcmin(cosmo,kpc, z)")
	println( "  Luminosity to flux :              lum2flux(cosmo,Lum, z)")
	println( "                                    flux2lum(cosmo,flux, z)")
end


function T2A(t)

	a = 1
end

function Hubble(z)

	H_a = 1
end

function DComoving(z) # comoving distance

	d_c = 1
end

function DTComoving(z) #transverse comoving distance

	d_tc = 1
end

function DAngular(z) # angular diameter distance

	d_a = 1
end

function DLuminosity(z) # luminosity distance

	d_lum = 1
end

function Z2T(z) # redshift -> time

	t = 1 
end

function T2Z(t) # time -> redshift

	z = 1
end

function Arcmin2Kpc(x)

	kpc = 1
end

function Kpc2Arcmin(x)

	arcmin = 1
end

function RhoCrit(z)

	rho = 1
end

function Delta(z)

	delta = 1
end

function Lum2Flux(L)

	F = 1
end

function Flux2Lum(F)

	L = 1
end

function DHubble(H0)

	d_hubble = 1
end

function THubble(H0)

	t_hubble = 1
end

end # module Cosmology
