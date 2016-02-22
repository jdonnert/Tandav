function luminosity_distance, z, H0, Omega_M, Omega_L
    return, (1+z)*transverse_comoving_distance(z, H0, Omega_M, Omega_L)
end
