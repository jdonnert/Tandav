pro plot_orbit

	@set_cgs

	common globals, gadget, cosmo

	M = [1d15, 1d10] *Msol / Gadget.mass

	e = 0.9 ; eccentricity

	x0 = 200D

	eta = 0 ; eccentric anomaly
	psi0 = -2 * atan( sqrt( (1+e) / (1-e) ) * tan(0.5*eta)  ) ; eq. 3.29

	pos = [ [0, 0, 0], [x0, 0, 0] ]

	r = sqrt( pos[0,1]^2 + pos[1,1]^2 + pos[2,1]^2 )
	
	v = sqrt( gadget.grav * M[0] * (1 + e) / r ) ; Keplers eq.
	vel = [ [0,0,0], [0,v,0] ] ; L = r x v in y direction

	L = [ pos[1,1]*vel[2,1] - vel[1,1]*pos[2,1], $
		  pos[2,1]*vel[0,1] - vel[2,1]*pos[0,1], $
		  pos[0,1]*vel[1,1] - vel[0,1]*pos[1,1]  ]
	L2 = L[0]^2  + L[1]^2 + L[2]^2 ; square of angular momentum L = r x vel
		
	C  = e/L2 * gadget.grav*M[0]
	a = L2 / (gadget.grav*M[0]*(1-e^2))

	Torbit = 2*!pi * sqrt(A^3 / (gadget.grav*M[0]))

	print, 'Torbit=', Torbit

	En = - gadget.grav * M[0] /2 /a * M[1]

	nstep = 1000

	x = make_array(nstep, /double)
	y = make_array(nstep, /double)

	for i=0,nstep-1 do begin

		psi = i*2*!pi / nstep
	
		r = a*(1-e^2) / (1 + e*cos(psi-psi0))

		x[i] = r * cos(psi)
		y[i] = r * sin(psi)
	end

	!p.multi[1] = 2

	plot, x, y, /iso, xrange=[-4000, 4000], yrange=[-4000,4000]

	oplot, x[0]*[1,1], y[0]*[1,1], psym=4

	oplot, [0,0], [0,0], psym=1

	nsnap = 933

	E = make_array(nsnap, /double)
	x = make_array(nsnap, /double)
	y = make_array(nsnap, /double)
	z = make_array(nsnap, /double)
	
	vx = make_array(nsnap, /double)
	vy = make_array(nsnap, /double)
	vz = make_array(nsnap, /double)

	for i = 0, nsnap-1 do begin
	
		fname = 'snap_'+strn(i, len=3, padc='0')

		pos = gadget.readsnap(fname, 'POS', head=head)
		vel = gadget.readsnap(fname, 'VEL')
	
		vx[i] = vel[0,1]
		vy[i] = vel[1,1]
		vz[i] = vel[2,1]

		x[i] = pos[0,1]
		y[i] = pos[1,1]
		z[i] = pos[2,1]

		;oplot, [1,1]*x[i], [1,1]*y[i], psym=3, col=color(0)

		E[i] =  0.5 * M[1] * (vx[i]^2+vy[i]^2+vz[i]^2) $
			- gadget.grav* M[1]*M[0]/ sqrt(x[i]^2+y[i]^2+z[i]^2)

	end

	oplot, x, y, col=color(1), psym=3
	
	; energy

	t = findgen(nsnap)/(nsnap-1) *head.time
	
	plot, t/Torbit, abs( (E-E[0]) / E[0] ), /ylog,  yrange=[1e-3,1]

	stop



	return 

end
