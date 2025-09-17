function [xx,tt] = integrate_ode_reachability(xx0,t0,tf,tstep, dynamics)
% This function integrate a specific set of ODEs
% INPUT 
%           -xx0: initial condition [Nx1]
%           - t0: initial time
%           - tf: final time
%           - tstep_max:max step size allowed
%           - spacecraft_data: spacecraft_data class
%           - rhs_model: integer indicating the right hand side 
% OUTPUT
%           - xx: integrated flow [Nxnt]
%           - tt: time vector [ntx1]
% Author: Antonio Rizza - 22/11/2024


tt = t0:tstep:tf;
if tt(end) ~= tf
    tt = [tt tf];
end

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
[tt , xx] = ode78(@(t,x) dynamics(t, x), tt, xx0, options);

end
