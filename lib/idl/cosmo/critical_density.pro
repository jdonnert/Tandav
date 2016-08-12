; critical density,  Mo, v.d.Bosch, White (3.63)
function critical_density, z, Omega_tot, H100, hbpar, Omega_M, Omega_r, Omega_L

	@set_cgs

	H0_cgs = H100 * hbpar

	H_z = hubble_parameter(z, H0_cgs, Omega_M, Omega_r, Omega_L)

	rho_crit = Omega_tot * 3 * H_z^2 / (8*!DPI*G)

	return, rho_crit
end
