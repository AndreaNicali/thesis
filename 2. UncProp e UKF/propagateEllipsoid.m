function [tt, xx, PHIf]  = propagateEllipsoid(x0, et_vec, mass, omega, C20, C22)
%propagate state and STM in the Earth Moon rotating planar adimensional frame from
%time t0 to time tf

% Initialize State Transition Matrix at t0
Phi0 = eye(6);

% Append to initial conditions the conditions for the STM
x0Phi0 = [x0; Phi0(:)];

% Perform integration
options = odeset('reltol', 1e-11, 'abstol', 1e-11);
[tt , xx] = ode78(@(t,x) dynamicsEllipsoidSTM(t, x, mass, omega, C20, C22), et_vec, x0Phi0, options);

% Extract state vector and State Transition Matrix
xf = xx(end,1:6)';
PHIf = reshape(xx(end,7:end),6,6);

xx = xx(:, 1:6);
end