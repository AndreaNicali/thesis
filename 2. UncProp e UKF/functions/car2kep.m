function [kep] = car2kep(r, v, mu)
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
%   [kep] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in 
% radians.
%
% INPUT:
% r            [3x1]  Position vector           [km]
% v            [3x1]  Velocity vector           [km/s] 
%
% OUTPUT:
% kep [6x1],which contains:
% a            [1x1]  Semi-major axis           [km]
% e            [1x1]  Eccentricity              [-]
% i            [1x1]  Inclination               [rad]
% OM           [1x1]  RAAN                      [rad]
% om           [1x1]  Periapsis anomaly         [rad]
% th           [1x1]  True anomaly              [rad]
% mu           [1x1]  Gravitational parameter   [km^3/s^2]
%
% CONTRIBUTORS:
%   Francesco Mazzetto
%   Tommaso Bustreo
%   Andrea Nicali
%
% VERSIONS:
%   2021-11-20: First version
%   2022-12-1: Added singularities

r_norm = norm(r); 
v_norm = norm(v);

h = cross(r, v);
h_norm = norm(h);

i = acos(h(3) / h_norm);

e = 1/mu*((v_norm^2 - mu/r_norm)*r - dot(r, v)*v);
e_norm = norm(e);

eps = 1/2*v_norm^2 - mu/r_norm;
a = -mu/(2*eps);

N = cross([0 0 1], h);
N_norm = norm(N);

if norm(N) == 0
    N = [1; 0; 0];
else
    N = N/norm(N);
end

if(N(2) >= 0) 
    OM = acos(N(1));
else
    OM = 2*pi - acos(N(1));
end

if(e(3) >= 0) 
    om = acos(dot(N, e)/(e_norm));
else
    om = 2*pi - acos(dot(N, e)/(e_norm));
end

v_r = dot(r, v)/r_norm;

if(v_r >= 0) 
    th = acos(dot(e, r)/(e_norm*r_norm));
else
    th = 2*pi - acos(dot(e, r)/(e_norm*r_norm));
end 

%singularities

if e_norm == 0
    om = 0;
end
if i == 0 % e = N = [1 0 0]
        OM = 0;

end
e=e_norm;
kep = [a e i OM om th];