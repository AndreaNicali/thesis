function [X,Y,Z] = plotOrbit(kepE1,mu,deltaTh)

% plotOrbit - Plots the arc of length deltaTh of the orbit described by kepE1.
% 
% CALLED FUNCTIONS:
%       kep2car.m
%       timeOfFlight.m
%
% PROTOTYPE:
%  [X,Y,Z] = plotOrbit(kepE1,mu,deltaTh)
%
% DESCRIPTION:
%  Plots the arc of the orbit described by a set of orbital
%  elements for a specific arc length.
%  Outputs the initial position in Cartesian elements.
%
% INPUT:
%   kepE1      [1x6]    Orbital elements                [km,rad]
%   mu         [1x1]    Gravitational parameter         [km^3/s^2]
%   deltaTh    [1x1]    Arc length                      [rad]
%
% OUTPUT:
%   X          [1xn]    X position                      [km]
%   Y          [1xn]    Y position                      [km]
%   Z          [1xn]    Z position                      [km]
%
% CONTRIBUTORS:
%   Francesco Mazzetto
%   Tommaso Bustreo
%   Andrea Nicali

% Standard condition
if deltaTh == 0
    [r, ~] = kep2car(kepE1, mu);
    X = r(1);
    Y = r(2);
    Z = r(3);
else
delta_t = timeOfFlight(kepE1(1),kepE1(2),kepE1(6),kepE1(6)+deltaTh,mu);

[r, v] = kep2car(kepE1, mu);
y0 = [r, v];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

[ ~, state ] = ode113( @(t,y) ode_2bp(t,y,mu), 0:3600:delta_t , y0, options );
X = state(:, 1);
Y = state(:, 2);
Z = state(:, 3);
end

