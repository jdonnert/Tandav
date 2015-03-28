; test the kick & drift factor numerical integration against the
; analytic solution for a EdS universe.
; also test spline interpolation

pro EdS_KDFactor ; in dt  integration

	Mpc2cm = 3.085d24
	yr2sec = 31556926d
	km2cm = 1d5
	t_unit = Mpc2cm/1e3 / 1d5

	H0 = 100D ; km/s/Mpc
	H0_cgs = H0 * km2cm / Mpc2cm
	H0_code = H0_cgs  * t_unit

	z0 = 140D
	z1 = 0D

	a0 = 1/(1+z0)
	a1 = 1/(1+z1)

	N = 512L
	i = lindgen(N)

	da = alog(a1/a0)/(N-1)
	a = a0 * exp(da * i)

	d_a = drift_factor_EdS_tandav( a0, a, H0_code)
	k_a = kick_factor_EdS_tandav( a0, a, H0_code)

	readcol, 'test_128', a_low, da_low, ka_low
	readcol, 'test_16384', a_high, da_high, ka_high

	plot, a_high, da_high, xtitle='a', ytitle='Drift & Kick factor [s]', $
		/ylog, /xlog

	oplot, a_high, ka_high

	plot, a_low, da_low, xtitle='a', ytitle='Drift & Kick factor [s]', $
		/ylog, /xlog

	oplot, a_low, ka_low

	err_da = abs(da_low - da_high)/da_high
	err_ka = abs(ka_low - ka_high)/ka_high

	plot, a_high, err_da, /ylog, psym=10
	oplot, a_high, err_ka, color=color(1)

stop
	return 
end

function drift_factor_EdS_tandav, a0, a1, H0

	return, -2D/H0 * (1/sqrt(a1) - 1/sqrt(a0))
end

function kick_factor_EdS_tandav, a0, a1, H0

	return, 2D/H0 * (sqrt(a1) - sqrt(a0))
end
