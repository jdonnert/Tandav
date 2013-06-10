function gadget_pressure, rho, u, gamma=gamma, xH=xH, z=z, h=h, uVel=uVel, uMass=uMass, uLength=uLength

    if not (exist(rho) or exist(U)) then begin
        print, "Pressure( rho, u, gamma=gamma, xH=xH, z=z, h=h, uVel=uVel, uMass=uMass, uLength=uLength )"
        return, -1
    end

    if not exist(gamma) then $
        gamma = 5D/3D

    rho_cgs = gadget_density(rho, z=z, h=h, xH=xH, uMass=uMass, uLength=uLength)

    return,  rho_cgs * u * (gamma-1) * uVel^2
end
