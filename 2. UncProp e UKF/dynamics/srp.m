function acc = srp(r_sc, t)

%Compute solar radiation pressure perturbating accelerations

%Recall Eros ID and Eros Position in EclipJ2000
id_eros = cspice_spkobj('kernels\erosephem_1999004_2002181.bsp', 1);
r_SCI = cspice_spkgeo(id_eros, t, 'ECLIPJ2000', 10);
r_SCI = r_SCI(1:3);

%Rotate to eros fixed frame
R = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t);

%Compute s/c sun vector
r_body = R*r_SCI + r_sc;
dist = norm(r_body)/astroConstants(2);

%Compute acceleration
F = 1361*(1/dist)^2;
reflectivity = 0.5;
ballistic = 0.01*10^-6; %literature
acc = F/astroConstants(5) * (1+reflectivity) * ballistic * r_body/norm(r_body);

end
