function [u_pert, P_man] = pertThrust(u_nom, sigma_magn, sigma_align, err_magn, err_align, err_theta)

uu_mag = norm(u_nom);
R_t2b = thrust2body(u_nom);

if sigma_magn == 0
    m_err = mvnrnd(0, sigma_magn^2);
    alpha_err = mvnrnd(0, sigma_align^2);
    theta_err = -pi + 2*pi*rand(1);
else
    m_err = err_magn;
    alpha_err = err_align;
    theta_err = err_theta;
end

u_pert = R_t2b*( (1+m_err)*uu_mag*[cos(alpha_err); sin(alpha_err)*cos(theta_err); sin(alpha_err)*sin(theta_err)]);

N = 0.25*(1+sigma_magn^2)*uu_mag^2;
P = exp(-sigma_align^2);

P_t = diag([2*N*(1+P^2)-P*uu_mag^2, N*(1-P^2), N*(1-P^2)]);
P_man = R_t2b*P_t*R_t2b';
end

function R = thrust2body(u)

x_t = u/norm(u);
ref = [0;0;1];
if abs(dot(x_t,ref)) > 0.9
    ref = [0;1;0]; 
end
y_t = cross(x_t, ref);
y_t = y_t/norm(y_t);
z_t = cross(x_t, y_t);

R = [x_t, y_t, z_t];
end