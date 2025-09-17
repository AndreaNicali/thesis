function [TT, XX] = compute_trajectory(x0, t_vec, u_vec, dyn, options)

% COMPUTE_TRAJECTORY Propagates a trajectory with impulsive maneuvers.
%
% INPUT:
%   x0        - initial state (6x1: position + velocity)
%   t_vec     - vector of impulse application times (length n+1)
%   u_vec     - 3xn matrix of impulses (∆v) to be applied (in m/s)
%   model_dyn - handle to the dynamics function (e.g., @(t,x) dynamics(t,x))
%   options   - options for the ODE solver (e.g., odeset)
%
% OUTPUT:
%   TT - concatenated time vector of the entire trajectory
%   XX - concatenated states (each row is a state)

    n = size(u_vec, 2);  % numero di segmenti di traiettoria
    x_curr = x0;
    TT = {};
    XX = {};
    % TT(1) = t_vec(1);
    % XX(1, :) = x0 + [zeros(3,1); u_vec(:, 1)];
    
    for i = 1:n
        x_start = x_curr + [zeros(3,1); u_vec(:, i)];
        t_span = t_vec(i):100:t_vec(i+1);
        [t_seg, x_seg] = ode78(@(t,x) dyn(t,x), t_span, x_start, options);
    
        TT{i} = t_seg;
        XX{i} = x_seg;

        x_curr = x_seg(end, :)';
    end

end

