; Pierpaoli 2001, Boehringer+ 2012
function overdensity_parameter, z, Omega_L, Omega_M

	cij = [ [546.67, -137.82, 94.083, -204.68,  111.51], 	$
		    [-1745.6, 627.22,  -1175.2, 2445.7,  -1341.7], 	$
		    [3928.8, -1519.3, 4015.8,  -8415.3, 4642.1],	$
			[-4384.8, 1748.7,  -5362.1, 11257.,  -6218.2], 	$
			[1842.3, -765.53, 2507.7, -5210.7, 2867.5] ]; Pierpaoli+ 01 Table 1

	x = Omega_M - 0.2
	y = Omega_L

	Delta = 0;

	for i = 0, 5 do begin
		for j = 0, 5 do begin

			Delta += cij[i][j] * x^i * y^j
		
		end
	end

	return, Omega_M * Delta
end
