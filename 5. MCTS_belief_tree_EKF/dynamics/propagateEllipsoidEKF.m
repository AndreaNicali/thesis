function [tt, xx, PHIf]  = propagateEllipsoidEKF(x0, et_vec, mass, omega, C20, C22, tao)
%propagate state and STM in the Earth Moon rotating planar adimensional frame from
%time t0 to time tf

% Initialize State Transition Matrix at t0
Phi0 = eye(9);

% Append to initial conditions the conditions for the STM
x0Phi0 = [x0; Phi0(:)];

% Perform integration
options = odeset('reltol', 1e-11, 'abstol', 1e-11);
[tt , xx] = ode78(@(t,x) dynamicsEllipsoidEKF(t, x, mass, omega, C20, C22, tao), et_vec, x0Phi0, options);

% Extract state vector and State Transition Matrix

PHIf = reshape(xx(end,10:end),9,9);

xx = xx(:, 1:9);
end