%This script propagate the trajectory around 433 Eros both in a truth 
%and a model environment. Truth environment comprehends spherical armonics
%perturbation up to order 15 and solar radiation pressure. In model
%environment, the asteroid is modeled as an ellipsoid hence it only has C20 and C22
%components of spherical armonic togheter with SRP.
%Moreover, in this script an analysis of the order of magnitude varying
%distance is performed. 
%To validate the models, propagation is confronted with gravity harmonics
%matlab function

clc
clear
close all

addpath(genpath('mice\'))
addpath(genpath('sgp4\'))
addpath(genpath('kernels\'))
addpath(genpath('functions\'))
addpath(genpath('3dModel\'))
addpath(genpath('dynamics\'))

cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\erosephem_1999004_2002181.bsp')
cspice_furnsh('kernels\erosatt_1998329_2001157_v01.bpc');

%% Dynamics around Eros in fixed body frame 

et_i = cspice_str2et('2000-08-01 T07:00:00');
et_f = cspice_str2et('2000-08-02 T10:00:00');

%Data for body frame propagation
T_rotation_eros = 5.27025547*3600; 
omega_body = 2*pi/T_rotation_eros*[0; 0; 1];
mass_eros = 6.687e15;
R_in2body_i = cspice_pxform('ECLIPJ2000', 'IAU_EROS', et_i); 

%Initial Conditions in EclipJ2000
r0 = [0.71; -45; 0];                    %[km]
v0 = [2.24; -0.035; 2.21]*10^(-3);      %[km/s]

%Rotation in body fixed frame
r0 = R_in2body_i*r0;
v0 = R_in2body_i*v0 - cross(omega_body, r0);

%%
%Propagation of the truth dynamics (SH up to 15th order)
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
[t, xx_rot_truth] = ode78(@(t,x) dynamics15(t, x, mass_eros, omega_body), [et_i, et_f], [r0; v0], options);

[t, xx_rot_matlab] = ode78(@(t,x) matlab_model_rotating_dynamics(t, x, mass_eros, omega_body, 15), [et_i, et_f], [r0; v0], options);

%Plot
figure(1)
plotEros
hold on
plot3(xx_rot_truth(:, 1), xx_rot_truth(:, 2), xx_rot_truth(:, 3), 'b', 'LineWidth', 1);
start_body = plot3(xx_rot_truth(1,1), xx_rot_truth(1, 2), xx_rot_truth(1, 3), 'r.', 'MarkerSize', 20);
ending_body = plot3(xx_rot_truth(end ,1), xx_rot_truth(end, 2), xx_rot_truth(end, 3), 'g.', 'MarkerSize', 20);
plot3(xx_rot_matlab(:, 1), xx_rot_matlab(:, 2), xx_rot_matlab(:, 3), 'r--', 'LineWidth', 1.5);

title('Body Frame Trajectory')
xlabel('X [Km]')
ylabel('Y [Km]')
zlabel('Z [Km]')
legend([start_body, ending_body], 'starting position', 'final position')

%% EROS ELLIPSOIDAL MODEL

% Semiaxis (in km) from Science
a = 20.591;
b = 5.711;
c = 5.332;

%Surface face refinement
n1 = 70;
n2 = 30;

%Create a triangular mesh of the ellipsoid
[F, V, N, C20, C22, A] = Ellipsoid(a, b, c, n1, n2);


%Propagation of the model trajectory
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
[t_model, xx_rot_mod] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), [et_i, et_f], [r0; v0], options);

[t, xx_rot_matlab] = ode78(@(t,x) matlab_model_rotating_dynamics(t, x, mass_eros, omega_body, 2), [et_i, et_f], [r0; v0], options);


%Plot
figure(3)
model_body = plot3(xx_rot_mod(:, 1), xx_rot_mod(:, 2), xx_rot_mod(:, 3), 'b', 'LineWidth', 1);
hold on
plot3(xx_rot_matlab(:, 1), xx_rot_matlab(:, 2), xx_rot_matlab(:, 3), 'r--', 'LineWidth', 1.5);
patch('Faces', F, 'Vertices', V);
plot3(xx_rot_mod(1,1), xx_rot_mod(1, 2), xx_rot_mod(1, 3), 'r.', 'MarkerSize', 20)
plot3(xx_rot_mod(end ,1), xx_rot_mod(end, 2), xx_rot_mod(end, 3), 'g.', 'MarkerSize', 20)
grid on
grid minor
axis equal


%% Perturbation magnitude of the truth model

distances = linspace(16, 1000, 1000);
acc_sph = zeros(length(distances), 1);
acc_srp = zeros(length(distances), 1);
acc_nominal = zeros(length(distances), 1);
acc_sph_upto2 = zeros(length(distances), 1);

for i = 1:length(distances)
    r = [0.5; 0.5; sqrt(2)/2]*distances(i);
    r = [1; 0; 0]*distances(i);

    acc_sph_vec = Spherical_armonics_perturbation(r, mass_eros, 15);
    acc_sph(i) = norm(acc_sph_vec);

    acc_sph_upto2_vec = Recursive_Spherical_armonics_perturbation_ellipsoid(r, mass_eros, C20, C22);
    acc_sph_upto2(i) = norm(acc_sph_upto2_vec);

    acc_srp_vec = srp(r, et_i);
    acc_srp(i) = norm(acc_srp_vec);

    acc_nominal(i) = norm(astroConstants(1)*mass_eros/norm(r)^3*r);
end

figure(5)
solradpr = loglog(distances, acc_srp*10^3, 'r', 'LineWidth', 1.5);
hold on
sh = loglog(distances, acc_sph*10^3, 'b', 'LineWidth', 1.5);
central = loglog(distances, acc_nominal*10^3, 'g', 'LineWidth', 1.5);
sh2 = loglog(distances, acc_sph_upto2*10^3, 'c', 'LineWidth', 1.5);

xlabel('Distance [km]')
ylabel('Perturbation magnitude [m/s^2]')
grid on
legend([central, sh, sh2, solradpr ], {'Central body', 'Spherical Armonics up to 15', 'Spherical Armonics Ellipsoid', 'Solar Rad Pressure'});
xlim([16, 1000])

