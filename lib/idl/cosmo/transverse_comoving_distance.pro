function transverse_comoving_distance, z, H0, Omega_M, Omega_L
    
    @set_cgs

    Omega_k = curvature_density(Omega_M, Omega_L)
    sqrt_Ok = sqrt(abs(Omega_k))
    d_hubble = c/double(H0)

    d_comov = comoving_distance(z, H0, Omega_M, Omega_L)

    if Omega_k eq 0 then $  ; 90% enough
        return, d_comov

    if Omega_k gt 0 then $
        return, d_hubble/sqrt_Ok * sinh(sqrt_Ok * d_comov/d_hubble )


    return, d_hubble/sqrt_Ok * sin(sqrt_Ok * d_comov/d_hubble )
end


