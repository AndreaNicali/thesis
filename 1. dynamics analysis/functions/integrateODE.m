function [xx,tt] = integrateODE(xx0, tt)

T_rotation_eros = 5.27025547*3600; 
omega_body = 2*pi/T_rotation_eros*[0; 0; 1];
mass_eros = 6.687e15;
C20 = -0.052461839309700;
C22 = 0.082399387985800;

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
[~ , xx] = ode78(@(t,x) model_rotating_dynamics(t, x, mass_eros, omega_body, C20, C22), tt, xx0, options);

end
