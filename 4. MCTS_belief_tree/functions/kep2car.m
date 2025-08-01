function [r,v] = kep2car(kep, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
% 
% PROTOTYPE:
% [r, v] = kep2car(kep)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
% kep [6x1],which contains:
% a            [1x1]  Semi-major axis           [km]
% e            [1x1]  Eccentricity              [-]
% i            [1x1]  Inclination               [rad]
% OM           [1x1]  RAAN                      [rad]
% om           [1x1]  Periapsis anomaly         [rad]
% th           [1x1]  True anomaly              [rad]
% mu           [1x1]  Gravitational parameter   [km^3/s^2]
%
% OUTPUT:
% r            [3x1]  Position vector           [km]
% v            [3x1]  Velocity vector           [km/s]   
%
% CONTRIBUTORS:
%   Francesco Mazzetto
%   Tommaso Bustreo
%   Andrea Nicali
%
% VERSIONS:
%   2021-11-20: First version
%   2022-12-1: Added singularities

a = kep(1); e = kep(2); i = kep(3);OM = kep(4); om = kep(5); th = kep(6);
% Straight semi-side
if e==1
    p=2*a;
else
    p=a*(1-e^2);
end
    % Radius from the ellipse fire
    r=p/(1+e*cos(th));

    % Radius and velocity in the perifocal reference system
    rpf=r*[cos(th); sin(th); 0];
    vpf=(sqrt(mu/p)).*[-sin(th); e+cos(th); 0];

    % Rotational matrix of OM around the versor k
    R_3_OM=[cos(OM)  sin(OM)  0;
        -sin(OM) cos(OM)  0;
        0        0        1];

    % Rotational matrix of i around the versor i'
    R_1_i=[1    0        0;
        0    cos(i)   sin(i);
        0    -sin(i)  cos(i)];

    % Rotational matrix of om around the versor k''
    R_3_om=[cos(om)   sin(om)  0;
        -sin(om)  cos(om)  0;
        0         0        1];

    % Complete matrix from perifocal system to earth centered intertial (ECI)
    T_pf_eci=(R_3_om*R_1_i*R_3_OM)';

    % Radius and velocity in ECI
    reci=T_pf_eci*rpf;
    veci=T_pf_eci*vpf;

    % Output
    r=reci;
    v=veci;
