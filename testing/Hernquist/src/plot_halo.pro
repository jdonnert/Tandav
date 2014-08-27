pro density_profile


	return
end

pro velocity_distribution_function
	
	r = sqrt(pos[0,*]^2 + pos[1,*]^2 + pos[2,*]^2)
	pot = hernquist_potential(r, a_hernq, mass)
	E = -1D * pot
	fE = hernquist_distribution_function(E, a_hernq, mass)

	plot, E, fE, psym=3, /xlog, /ylog

	v = sqrt(vel[0,*]^2 + vel[1,*]^2 + vel[2,*]^2)

	Etot = 0.5 * v^2 + pot
	one = Etot * 0 + 1

	tmp = bin_arr(one, pos=-Etot, bin_pos=bin_pos, nbins=100, cnt=cnt, /log)

	oplot, bin_pos, cnt/mass, psym=10
	
	return
end

function hernquist_density_profile, r, a, mass

	common  globals, tandav, cosmo

	return, Mass/(2*!pi) * a/r / (r+a)^3
end

function hernquist_potential, r, a, mass

	common  globals, tandav, cosmo

	return, -tandav.Grav * mass / (r+a)
end

function hernquist_distribution_function, E, a, mass

	common  globals, tandav, cosmo
	
	prefac = 1D / (sqrt(2) * (2*!pi)^3 * (tandav.Grav*mass*a)^(1.5) )

	q2 = a * E / (tandav.Grav * mass)

	fE = prefac * mass * sqrt(q2) / (1-q2)^2 $
		*( (1 - 2*q2) * (8*q2*q2 - 8*q2 - 3) $
		+ (3*asin(sqrt(q2)))/sqrt(q2*(1-q2)) )

	return, fE
end
