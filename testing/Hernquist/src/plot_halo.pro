
;for i =0,25 do begin & density_profile, snap=i & wait, 0.2 & end
;for i =0,25 do begin & potential_profile, snap=i & wait, 0.2 & end

pro show_halo, snap=snap, periodic=periodic

	common globals, tandav, cosmo

	if not keyword_set(snap) then $
		snap = 0
	
	fname = 'snap_'+strn(snap, len=3, padc='0')

	pos = tandav.readsnap(fname, 'POS', head=head)

	plot, pos[0,*], pos[1,*], psym=3, /iso, xrange=minmax(pos[0,*]), $
		yrange=minmax(pos[1,*]), xstyle=1, ystyle=1

	return
end

pro density_profile, nsnap=nsnap, periodic=periodic
	
	; we expect a little evolution here, because f(E) in the ICs 
	; does not account for softening (Barnes 2012)

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

	for snap = 0, nsnap-1 do begin

		fname = 'snap_'+strn(snap, len=3, padc='0')

		pos = tandav.readsnap(fname, 'POS', head=head)

		print, median(pos[0,*]), median(pos[1,*]), median(pos[2,*])

		pos[0,*] -= median(pos[0,*])
		pos[1,*] -= median(pos[1,*])
		pos[2,*] -= median(pos[2,*])

		r = sqrt(pos[0,*]^2 + pos[1,*]^2 + pos[2,*]^2)

		one = r * 0 + 1

		tmp = bin_arr(one, pos=r , bin_pos=bin_pos, nbins=100, cnt=cnt, /log)

		Vshell = make_array(100, /double, val=0) ; volume
		Vshell[0] = 4D*!pi/3D * bin_pos[0]^3

		for i = 0L, 99 do $
			Vshell[i] =  4D*!pi/3D * (bin_pos[i]^3 - bin_pos[i-1]^3)

		rho = cnt * head.massarr[1] / Vshell

		oplot, bin_pos, rho, psym=10

	end

	return
end

pro potential_profile, snap=snap

	common globals, tandav, cosmo

	mass = 1d5 
	a_hernq = 924D 		

	rmin = 1d
	rmax = 1d6

	N = 1000L

	di = alog10(rmax/rmin) / (N-1)

	r = rmin * 10D^(lindgen(N)*di)

	pot_analytic = hernquist_potential(r, a_hernq, mass)

	plot, r, -pot_analytic, /ylog, /xlog, xrange=[rmin, rmax]

	if not keyword_set(snap) then $
		snap = 0

	fname = 'snap_'+strn(snap, len=3, padc='0')

	pos = tandav.readsnap(fname, 'POS', head=head)

	pos[0,*] -= median(pos[0,*])
	pos[1,*] -= median(pos[1,*])
	pos[2,*] -= median(pos[2,*])

	r = sqrt(pos[0,*]^2 + pos[1,*]^2 + pos[2,*]^2)

	gpot = tandav.readsnap(fname, "GPOT")

	oplot, r, -gpot, psym=3, color=7839259

	return 
end

pro conservation

	common globals, tandav, cosmo

	nsnap = 19

	X = make_array(nsnap, val=0D)
	Y = make_array(nsnap, val=0D)
	Z = make_array(nsnap, val=0D)
	E = make_array(nsnap, val=0D)
	P = make_array(nsnap, val=0D)

	AngP = make_array(3, nsnap, val=0D)
	time = make_array(nsnap, val=0D)

	for snap = 0, nsnap-1 do begin
	
		fname = 'snap_'+strn(snap, len=3, padc='0')

		pos = tandav.readsnap(fname, 'POS', head=head)
		vel = tandav.readsnap(fname, 'VEL', head=head)
		pot = tandav.readsnap(fname, 'GPOT')

		v = length(vel)

		npart = head.npart[1]
		mpart = head.massarr[1]

		X[snap] = total(pos[0,*])/npart
		Y[snap] = total(pos[1,*])/npart
		Z[snap] = total(pos[2,*])/npart

		E[snap] = total(0.5 * mpart * v^2 + pot)
		
		P[snap] = mpart * total(v)

		AngP[0, snap] = mpart * total(pos[1,*]*vel[2,*] - pos[2,*]*vel[1,*])
		AngP[1, snap] = mpart * total(pos[2,*]*vel[0,*] - pos[0,*]*vel[2,*])
		AngP[2, snap] = mpart * total(pos[0,*]*vel[1,*] - pos[1,*]*vel[0,*])

		time[snap] = head.time
		
	end

	plot, time, (X-X[0])/X[0], yrange=[-0.1, 0.1]
	oplot, time, (Y-Y[0])/Y[0]
	oplot, time, (Z-Z[0])/Z[0]
		
	oplot, time, (E-E[0])/E[0], col=11759733

	oplot, time, (P-P[0])/P[0], color=7839259

	oplot, time, (AngP[0,*]-AngP[0,0])/AngP[0,0], linestyle=0, color=155609
	oplot, time, (AngP[1,*]-AngP[1,0])/AngP[1,0], linestyle=1, color=155609
	oplot, time, (AngP[2,*]-AngP[2,0])/AngP[2,0], linestyle=2, color=155609

	xyouts, 1, -0.05, "R"
	xyouts, 2, -0.05, "E", col=11759733
	xyouts, 3, -0.05, "P", col=7839259
	xyouts, 4, -0.05, "AngP", col=155609

	stop

	plot, time, R
	plot, time, E
	plot, time, P
	plot, time, abs(AngP[1,*])
	oplot, time, abs(AngP[0,*])
	oplot, time, abs(AngP[2,*])


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

	oplot, bin_pos/tandav.grav*a_hernq/Mass,  $
		(bin_pos)^2*cnt/mass/2d37/head.npart[1], psym=10
	
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




