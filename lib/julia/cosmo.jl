module Cosmology

import CGS_Units
import Base.show
import Base.convert

type Params
	Model 		:: "Concordance" # Name of model
	t0 			:: Float64	# Age of universe
	H0			:: Float64	# Hubble constant at z=0
	Omega_b		:: Float64	# Baryon density parameter
	Omega_c		:: Float64	# Dark matter density parameter
	Omega_L		:: Float64	# Dark energy density parameter
	Omega_r		:: Float64	# Photon energer density parameter
	sigma_8		:: Float64	# Fluctuation amplitude at 8/h Mpc
	Delta_r2	:: Float64	# Curvature fluctuation amplitude
	n_s			:: Float64  # scalar spectral index
	tau			:: Float64  # reionisation optical depth
	w			:: Float64  # equation of state parameter
	z_star		:: Float64  # redshift of decoupling
	t_star		:: Float64  # Age of decoupling
	z_reion		:: Float64  # Redshift of reionisation

	# Derived Parameters	

	Omega_M		:: Float64  # Matter density
	Omega_k		:: Float64  # curvature density
	Omega_tot	:: Float64  # total density
	hbpar		:: Float64  # h parameter
	H100		:: Float64  # hubble 100 constant
	d_hubble	:: Float64  # hubble distance
	t_hubble	:: Float64  # hubble time 

end # type cosmo

function show!(cosmo::Params)

	println("Cosmology $(cosmo.Model) : ")
	println("   H0       = $(cosmo.H0/3.2407765e-20) km/s/Mpc")
	println("   Omega_L  = $(cosmo.Omega_L)")
	println("   Omega_M  = $(cosmo.Omega_M)")
	println("   Omega_b  = $(cosmo.Omega_b)")
	println("   Omega_r  = $(cosmo.Omega_r)")
	println("   sigma_8  = $(cosmo.sigma_8)")
	println("   Horizon  = $(cosmo.d_hubble/(1e3*kpc2cm))")

	println( "Methods (all cgs) :")
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

	return, cosmo
end

function convert(i::Type{Int}, x::Params) # custom constructor

	if i = 0
	
	elseif i = 1

	elseif i = 2
	
	else

	end

end


function t2a(cosmo::Params, Real::t)

	a = 1
end

function hbpar(cosmo::Params, Real::z)

	H_a = 1
end

function d_cmv(cosmo::Params, Real::z) # comoving distance

	d_c = 1
end

function d_tcmv(cosmo::Params, Real::z)

	d_tc = 1
end

function d_ang(cosmo::Params, Real::z)

	d_a = 1
end

function d_lum(cosmo::Params, Real::z)

	d_lum = 1
end

function z2t(cosmo::Params, Real::z)

	t = 1 
end

function t2z(cosmo::Params, Real::t)

	z = 1
end

function arcmin2kpc(cosmo::Params, Real::x)

	kpc = 1
end

function kpc2arcmin(cosmo::Params, Real::x)

	arcmin = 1
end

function rho_crit(cosmo::Params, Real::z)

	rho = 1
end

function delta(cosmo::Params, Real::z)

	delta = 1
end

function lum2flux(cosmo::Params, Real::L)

	F = 1
end

function flux2lum(cosmo::Params, Real::F)

	L = 1
end

end # module cosmology
