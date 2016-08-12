;  Mo, v.d.Bosch, White (2.62, 3.75), Boehringer 2012 (5)
function E, z, Omega_M, Omega_r, Omega_L

	Omega_tot = Omega_M+Omega_r+Omega_L

	return, sqrt(Omega_L + Omega_M*(1+z)^3 + Omega_r*(1+z)^4 $
						 + (1D - Omega_tot)*(1+z)^2) 
end

;  Mo, v.d.Bosch, White (3.74)
function hubble_parameter, z, H0, Omega_M, Omega_r, Omega_L
    
    return, H0 / E(z, Omega_M, Omega_r, Omega_L)   
end


