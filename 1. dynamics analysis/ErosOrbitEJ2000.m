clc
clear
close all

clc
clear
close all

addpath(genpath('kernels\'))
addpath(genpath('functions\'))
addpath(genpath('Model3d\'))
addpath(genpath('dynamics\'))

cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\erosephem_1999004_2002181.bsp')
cspice_furnsh('kernels\erosatt_1998329_2001157_v01.bpc');

%% Orbit of 433 Eros around the Sun, visualization of body frame
%cacacaaacacacac
id_eros = cspice_spkobj('kernels\erosephem_1999004_2002181.bsp', 1);

%Define a time vector for the propagations
et_i = cspice_str2et('2000-08-01 T00:00:00');
et_f = cspice_str2et('2001-08-01 T00:00:00');
n_step = 3600;
step = (et_f-et_i) / n_step;
et_vec = et_i:step:et_f;

eros_eph = zeros(n_step, 6);

%Get EROS position in EclipJ2000
for i = 1:n_step
    eros_eph(i, :) = cspice_spkgeo(id_eros, et_vec(i), 'ECLIPJ2000', 10);
end

%Plot
figure(1)
plot_eros_eph = plot3(eros_eph(:, 1), eros_eph(:, 2), eros_eph(:, 3), 'r', 'LineWidth', 1.5);
hold on
plot_sun = plot3(0, 0, 0, 'r.', 'MarkerSize', 30);
grid on
grid minor
title('433 Eros Orbit around the Sun')

%Plot Reference frames
rotmat = cspice_pxform('ECLIPJ2000', 'IAU_EROS', et_i);
X_in = [1; 0; 0]*10^8;
Y_in = [0; 1; 0]*10^8;
Z_in = [0; 0; 1]*10^8;

X_body = rotmat'*X_in;
Y_body = rotmat'*Y_in;
Z_body = rotmat'*Z_in;

Inertial = plot3([0, X_in(1)], [0, X_in(2)], [0, X_in(3)], 'k', 'LineWidth', 2);
plot3([0, Y_in(1)], [0, Y_in(2)], [0, Y_in(3)], 'k', 'LineWidth', 2);
plot3([0, Z_in(1)], [0, Z_in(2)], [0, Z_in(3)], 'k', 'LineWidth', 2);

body_fixed = plot3([0, X_body(1)], [0, X_body(2)], [0, X_body(3)], 'r', 'LineWidth', 2);
plot3([0, Y_body(1)], [0, Y_body(2)], [0, Y_body(3)], 'r', 'LineWidth', 2);
plot3([0, Z_body(1)], [0, Z_body(2)], [0, Z_body(3)], 'r', 'LineWidth', 2);

legend([Inertial, body_fixed], {'Inertial Frame', 'Body-Fixed Frame'});