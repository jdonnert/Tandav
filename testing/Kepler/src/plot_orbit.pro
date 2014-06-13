pro plot_orbit

	common globals, gadget, cosmo

	M = [1d11, 1d0] / gadget.mass 

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

	readcol, "test", cnt, x,y,z,vx,vy,vz,ax,ay,az, t
	
	oplot, x, y, psym=3, col=color(0)

	; energy

	E = 0.5 * M[1] * (vx^2+vy^2+vz^2) - gadget.grav*M[0]/ sqrt(x^2+y^2+z^2)
	
	plot, abs(E), /ylog, psym=3

	stop



	return 

end
