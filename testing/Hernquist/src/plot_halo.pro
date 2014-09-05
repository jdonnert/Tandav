pro density_profile, snap=snap

	common globals, tandav, cosmo

	mass = 1d5 
	a_hernq = 924D 		

	rmin = 1d
	rmax = 1d6

	N = 1000L

	di = alog10(rmax/rmin) / (N-1)

	r = rmin * 10D^(lindgen(N)*di)

	rho_analytic = hernquist_density_profile(r, a_hernq, mass)

	plot, r, rho_analytic, /ylog, /xlog, xrange=[rmin, rmax]

	if not keyword_set(snap) then $
		fname = 'IC_Hernquist_Halo' $
	else $
		fname = 'snap_'+strn(snap, len=3, padc='0')

	pos = tandav.readsnap(fname, 'POS', head=head)

	pos[0,*] -= median(pos[0,*])
	pos[1,*] -= median(pos[1,*])
	pos[2,*] -= median(pos[2,*])

	r = sqrt(pos[0,*]^2 + pos[1,*]^2 + pos[2,*]^2)

	one = r * 0 + 1

	tmp = bin_arr(one, pos=r , bin_pos=bin_pos, nbins=100, cnt=cnt, /log)

	Vshell = make_array(100, /double, val=0)
	Vshell[0] = 4D*!pi/3D * bin_pos[0]^3

	for i = 0L, 99 do $
		Vshell[i] =  4D*!pi/3D * (bin_pos[i]^3 - bin_pos[i-1]^3)

	rho = cnt * head.massarr[1] / Vshell

	oplot, bin_pos, rho, psym=10

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

	plot, E/tandav.grav*a_hernq/Mass, fE, /ylog, xstyle=1, xrange=[0,1.3]

	; snap

	fname = 'IC_Hernquist_Halo'

	pos = tandav.readsnap(fname, 'POS', head=head)
	vel = tandav.readsnap(fname, 'VEL')

	r = sqrt(pos[0,*]^2 + pos[1,*]^2 + pos[2,*]^2)
	
	pot = hernquist_potential(r, a_hernq, mass)

	v = sqrt(vel[0,*]^2 + vel[1,*]^2 + vel[2,*]^2)

	Etot = 0.5 *  v^2 + pot
	one = Etot * 0 + 1

	tmp = bin_arr(one, pos=abs(Etot) , bin_pos=bin_pos, nbins=100, cnt=cnt, /log)

	oplot, bin_pos/tandav.grav*a_hernq/Mass,  (bin_pos)^2*cnt/mass/2d37/head.npart[1], psym=10
	
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
