function dy = dynamicsEllipsoid( t, y, mass, omega, C20, C22)

%Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);
G = astroConstants(1);

%Compute model spherical armonic perturbation and solar radiation pressure
acc_shp = Recursive_Spherical_armonics_perturbation_ellipsoid(y, mass, C20, C22);

% Set the derivatives of the state
dy = [ v
 (-G*mass/rnorm^3)*r - 2*cross(omega, v) - cross(omega, cross(omega, r)) + acc_shp];

end
