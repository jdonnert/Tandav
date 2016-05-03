pro softening

	common globals, tandav, cosmo

	N = 1024
	rmin = 1e-4
	rmax = 1d4

	i = lindgen(N)
	
	dr = alog10(rmax/rmin) / (N-1)
	r = rmin * 10D^(i * dr)

	eps = 1.85D

	h = eps * (-105D)/32D

	kernel = make_array(N, /double, val=0)
	force = make_array(N, /double, val=0)
	poten = make_array(N, /double, val=0)
	
	for i = 0,N-1 do begin

		kernel[i] = K1(r[i], eps)

		force[i] = force_K1(1,1,r[i], eps)
		
		poten[i] = potential_K1(1, r[i], eps)
	end

	!p.multi[1] = 3

	plot, r, kernel, xrange=[rmin, 3], xstyle=1, xtitle='r [cm]', $
		ytitle=textoidl('\rho [g/cm^3]')
	plot, r, force, xtitle='r [cm]', ytitle='F [g cm/s^2]', $
		xrange=[rmin, 3]
	plot, r, -poten,  xtitle='r [cm]', xrange=[rmin, 3], $
		ytitle=textoidl('\phi [g cm^2/s^2]')

	return 
end

function K1, r, h ; Dehnen 2001, fig 1 & 2

	u = r/h

	if u gt 1 then $
		return,  0

	return, 105D/64/!pi/h^3 * (5 - 9*u^2) * (1 - u^2)

end

function force_K1, M0, M1, r, h
	
	grav_const = 6.6720000e-08

	u = r/h

	if u gt 1 then $
		return,  grav_const * M0 * M1 / r^2

	rinv2 = u * (135*u^4 - 294*u^2 + 175D)/16D/h^2 

	return,  grav_const * M0 * M1 * rinv2
	
end

function potential_K1, M, r, h

	grav_const = 6.6720000e-08

	u = r/h

	if r gt h then $
		return, -grav_const * M  / r

	rinv = ( 45D*u^6  - 147D*u^4 + 175D * u^2 - 105D) /32D/h

	return, grav_const * M * rinv
end

function WC2, r, h

	u = r/h

	if u gt 1 then $
		return, 0

	return, 21D/(2D*!pi) /h^3 * (1-u)^4 * (1+4*u)
end

function force_WC2, M0, M1, r, eps

	grav_const = 6.6720000e-08

	h = 3D/eps; 

	u = double(r/h)

	if r gt h then $
		return,  grav_const * M0 * M1 / r^2

	rinv2 = (14*u - 84*u^3 + 140*u^4 - 90*u^5 + 21*u^6) / h^2

	return, grav_const * M0 * M1  * rinv2
end

function potential_WC2, M, r ,eps

	grav_const = 6.6720000e-08

	h = 3/eps 

	if r ge h then $
		return, -grav_const * M  / r

	u = double(r/h)

	rinv = (7*u^2 - 21*u^4 + 28*u^5 - 15*u^6 + 3*u^7 - 3D) / h

	print, u, r, h,rinv 

	return,  grav_const * M * rinv
end
