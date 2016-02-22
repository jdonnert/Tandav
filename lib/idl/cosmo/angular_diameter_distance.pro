function angular_diameter_distance, z, H0, Omega_M, Omega_L

    d_transcomov = transverse_comoving_distance(z, H0, Omega_M, Omega_L)

    return, d_transcomov/(1D0 +z)
end
