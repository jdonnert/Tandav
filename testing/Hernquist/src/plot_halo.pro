pro density_profile


	return
end

pro velocity_distribution_function

	common globals, tandav, cosmo

	; model

	mass = 1d15 
	a_hernq = 924D 		

	Emax = tandav.grav*mass/a_hernq
	Emin = 1d-3 * Emax

	N =  1000L
	di = alog10(Emax/Emin) / (N-1)
	pot = -Emin * 10D^(lindgen(N)*di)
	E = -1D * pot
	fE = hernquist_distribution_function(E, a_hernq, mass)

	plot, E/tandav.grav*a_hernq/Mass, fE*(tandav.grav*mass*a_hernq)^(1.5),  $
		/ylog, xstyle=1

	; snap

	fname = 'IC_Hernquist_Halo'

	pos = tandav.readsnap(fname, 'POS')
	vel = tandav.readsnap(fname, 'VEL')

	r = sqrt(pos[0,*]^2 + pos[1,*]^2 + pos[2,*]^2)
	
	pot = hernquist_potential(r, a_hernq, mass)

	v = sqrt(vel[0,*]^2 + vel[1,*]^2 + vel[2,*]^2)

	Etot = 0.5 * v^2 + pot
	one = Etot * 0 + 1

	tmp = bin_arr(one, pos=-Etot , bin_pos=bin_pos, nbins=500, cnt=cnt)

	oplot, bin_pos,  cnt/mass, psym=10
	stop
	return
end

function hernquist_density_profile, r, a, mass

	return, Mass/(2*!pi) * a/r / (r+a)^3
end

function hernquist_potential, r, a, mass

	common  globals, tandav, cosmo

	return, -tandav.Grav * mass / (r+a)
end

function hernquist_distribution_function, E, a, mass

	common  globals, tandav, cosmo
	
	prefac = 1D / (sqrt(2D) * (2D*!PI)^3 * (tandav.Grav*mass*a)^(1.5D) )

	q2 = a * E / (tandav.Grav * mass)

	fE = prefac * mass * sqrt(q2) / (1D - q2)^2 $
		*( (1 - 2*q2) * (8D*q2*q2 - 8D*q2 - 3D) $
		+ (3D*asin(sqrt(q2)))/sqrt(q2*(1D - q2)) )

	return, fE
end
