function dy = dynamicsTrue( t, y, mass, omega)

%Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);
G = astroConstants(1);

% Add spherical armonics and solar radiation pressure
acc_shp = rhs_true(y, mass, 15);

% Set the derivatives of the state
dy = [ v
 (-G*mass/rnorm^3)*r - 2*cross(omega, v) - cross(omega, cross(omega, r)) + acc_shp];

end
