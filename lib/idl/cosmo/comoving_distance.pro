function comoving_distance_numerical, z, H0, Omega_M, Omega_L

    @set_cgs

    print, "ERROR not implemented yet!"
    stop

    result =  1 

    return, c/H0 * result
end

; Wickramasinghe+ 2010
function comoving_distance_approx, z, H0, Omega_L 
    
    @set_cgs

    alpha = 1 + 2*Omega_L/ (1 - Omega_L) / (1+z)^3
    x = alog(alpha + sqrt(alpha^2 - 1))
    Sigma_z = 3 * x^(1./3.) * 2^(2./3.) * (1 - x^2/252. + x^4/21060)

    alpha =  1 + 2*Omega_L/ (1 - Omega_L)
    x = alog(alpha + sqrt(alpha^2 - 1))
    Sigma_0 = 3 * x^(1./3.) * 2^(2./3.) * (1 - x^2/252. + x^4/21060)

    prefac = c/(3*H0)/(Omega_L^(1./6.)*(1-Omega_L)^(1./3.))

    return, prefac * (Sigma_0 - Sigma_z)
end

function comoving_distance, z, H0, Omega_M, Omega_L
    
    Omega_k = curvature_density(Omega_M, Omega_L)

    if Omega_k eq 0 and Omega_L ne 0 then $
        return, comoving_distance_approx(z, H0, Omega_L) 

    return, comoving_distance_numerical(z, H0, Omega_M, Omega_L)
end
