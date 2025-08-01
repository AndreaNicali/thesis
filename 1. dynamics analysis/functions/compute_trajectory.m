function [TT, XX] = compute_trajectory(x0, t_vec, u_vec, model_dyn, options)

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
    
    t_step = 200;
    TT = zeros(t_step*n, 1);
    XX = zeros(t_step*n, 6);
    idx = 1;
    
    for i = 1:n
        x_curr = x_curr + [zeros(3,1); u_vec(:, i)];
        t_span = linspace(t_vec(i), t_vec(i+1), t_step);
        [t_seg, x_seg] = ode78(@(t,x) model_dyn(t,x), t_span, x_curr, options);
    
        range = idx:idx + length(t_seg) - 1;
        TT(range) = t_seg(:);
        XX(range, :) = x_seg;
        idx = idx + length(t_seg);
    
        x_curr = x_seg(end, :)';
    end
    
    % Se t_seg ha meno di t_step elementi (es. per interruzioni), tronca
    TT = TT(1:idx-1);
    XX = XX(1:idx-1, :);
