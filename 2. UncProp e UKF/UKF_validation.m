%In this script the validation of UKF is performed. Measurements are taken
%from the true dynamics (15th order SH) as the line of sight towards 100
%features randomly distributed around the asteroid. The line of sight is
%computed as coaltitude and azimuth angle in the feature topocentric frame.
%A measurement can be used if the feature is visible from the spacecraft
%(inside the camera field of view and with a proper emission angle) and
%illuminated (with a proper incidence angle).

clc
clear
close all

addpath(genpath('kernels\'))
addpath(genpath('functions\'))
addpath(genpath('3dModel\'))
addpath(genpath('dynamics\'))
addpath(genpath('Uncertainties\'))

cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\erosephem_1999004_2002181.bsp')
cspice_furnsh('kernels\erosatt_1998329_2001157_v01.bpc');

%% UT propagation 

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

x_mean = [21.095697268161683; -31.527866008753950; -23.491639672861801;  -0.011418159898777;  -0.009247442529600;  0.001997382918943];

%Initial Covariance
sigma_pos = 0.5e-1;
sigma_vel = 0.3*10^(-4);
P_initial = [eye(3)*sigma_pos^2, zeros(3); zeros(3),  eye(3)*sigma_vel^2];
n_sample = 1;
x_0_vec = zeros(6, n_sample);

for i = 1:size(x_0_vec, 2)
    x_0_vec(:, i) = mvnrnd(x_mean, P_initial);
end

%Epochs
et_i = cspice_str2et('2000-08-01 T07:00:00');
et_f = cspice_str2et('2000-08-01 T11:00:00');
et_step = 60;
et_vec = et_i:60:et_f;

% Semiaxis (in km) of the asteroid model 
a = 20.591;
b = 5.711;
c = 5.332;
n1 = 40;
n2 = 20;

% Asteroid model generation
[F, V, N, C20, C22, A] = Ellipsoid(a, b, c, n1, n2);

%Data for body frame propagation
T_rotation_eros = 5.27025547*3600; 
omega_body = 2*pi/T_rotation_eros*[0; 0; 1];
mass_eros = 6.687e15;

%Set some assumed to be known features on the asteroid surface, the index
%corresponds to the face
features.nav_index = [   7,   12,   17,   24,   53,   62,   79,   86,   88,  116,  120, ...
  123,  162,  166,  170,  174,  188,  201,  217,  223,  236,  246, ...
  247,  270,  307,  326,  338,  340,  362,  370,  397,  405,  409, ...
  413,  414,  429,  465,  513,  536,  576,  588,  593,  622,  642, ...
  657,  674,  707,  708,  719,  749,  761,  763,  781,  786,  817, ...
  819,  836,  847,  849,  850,  867,  881,  901,  907,  911,  912, ...
  934,  953,  999, 1004, 1009, 1025, 1040, 1041, 1052, 1066, 1090, ...
 1102, 1105, 1106, 1127, 1205, 1223, 1226, 1227, 1264, 1266, 1282, ...
 1338, 1343, 1348, 1358, 1385, 1402, 1410, 1415, 1431, 1447, 1456, ...
 1466];

known_map = ones(size(F, 1), 1);

% FOV of the Navigation camera
fov1 = 20*pi/180;
fov2 = 20*pi/180;

% Sintonize spacecraft_data struct
spacecraft_data.data_guidance.fov1 = fov1;
spacecraft_data.data_guidance.fov2 = fov2;
spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.features = features;
spacecraft_data.data_asteroids.known_map = known_map;
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;
spacecraft_data.data_asteroids.omega = omega_body;
spacecraft_data.data_asteroids.mass = mass_eros;

%Iniatilize measurements struct
measurements = struct();

%Measurements noise
sigma1 = 100; %37.2
sigma2 = 100; %90
sigma_angular = ( (sigma1/3600*pi/180)^2 + (sigma2/3600*pi/180)^2 ) ;
sigma_range = ( (1/1000)^2 + (10/1000)^2 ) ;

%Integration of truth and model dynamics in the given time interval
[~ , xx_model] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), et_vec, x_mean, options);

%Process noise
tao = et_step;
sigma_acc = 5e-9;

q11 = (1/4)*tao^4;
q12 = (1/2)*tao^3;
q22 = tao^2;

Qk = sigma_acc^2 * [ ...
    q11, 0,   0,   q12, 0,   0;
    0,   q11, 0,   0,   q12, 0;
    0,   0,   q11, 0,   0,   q12;
    q12, 0,   0,   q22, 0,   0;
    0,   q12, 0,   0,   q22, 0;
    0,   0,   q12, 0,   0,   q22];

%Measurement noise
sigma_meas = sigma_angular;

%Initialize storing vectors
sigma_r = zeros(size(et_vec));
sigma_v = zeros(size(et_vec));
sigma_r(1) = sqrt(trace(P_initial(1:3, 1:3)));
sigma_v(1) = sqrt(trace(P_initial(4:6, 4:6)));

filter_error_vec = zeros(size(x_0_vec, 2), length(et_vec));
filter_error_x_vec = zeros(size(x_0_vec, 2), length(et_vec));
filter_error_y_vec = zeros(size(x_0_vec, 2), length(et_vec));
filter_error_z_vec = zeros(size(x_0_vec, 2), length(et_vec));
sigma_r_vec = zeros(size(x_0_vec, 2), length(et_vec));
sigma_x_vec = zeros(size(x_0_vec, 2), length(et_vec));
sigma_y_vec = zeros(size(x_0_vec, 2), length(et_vec));
sigma_z_vec = zeros(size(x_0_vec, 2), length(et_vec));

sigma_r_vec(:, 1) = sqrt(trace(P_initial(1:3, 1:3)));
sigma_x_vec(:, 1) = sqrt(P_initial(1, 1));
sigma_y_vec(:, 1) = sqrt(P_initial(2, 2));
sigma_z_vec(:, 1) = sqrt(P_initial(3, 3));

P_tutte = zeros(6,6, length(et_vec));
P_tutte(:, :, 1) = P_initial;

h = waitbar(0, 'Propagating Uncertainties...'); 


for j = 1:size(x_0_vec, 2)
    
    x_0 = x_0_vec(:, j);

    r_start = xx_model(1, 1:3)';
    v_start = xx_model(1, 4:6)';
    xx(1, :) = xx_model(1, :);
    P0 = P_initial;

    [~ , xx_true] = ode78(@(t,x) dynamics15(t, x, mass_eros, omega_body), et_vec, x_0, options);

    %Obtain measurements from true trajectory
    [measurements] = feature_measurements(xx_true, et_vec, spacecraft_data);
    
    %Perturb measurements 
    for i = 1:length(measurements.coaltitude)
        measurements.coaltitude.val(i) = mvnrnd(measurements.coaltitude.val(i), sigma_angular);
        measurements.azimuth.val(i) = mvnrnd(measurements.azimuth.val(i), sigma_angular);
    end

    for i = 1:(length(et_vec)-1)
    
        if any(et_vec == et_vec(i+1))
            [xx_filtered, P_filtered] = UKF([r_start; v_start], [et_vec(i), et_vec(i+1)], measurements, P0, spacecraft_data, sigma_meas, Qk);
            
            r_start = xx_filtered(end, 1:3)';
            v_start = xx_filtered(end, 4:6)';
            P0 = P_filtered(:, :, end);
    
            xx(i+1, :) = xx_filtered(end, :);
            sigma_r(i+1) = sqrt(trace(P0(1:3, 1:3)));
            sigma_v(i+1) = sqrt(trace(P0(4:6, 4:6)));
            P_tutte(:, :, i+1) = P0;
            sigma_x_vec(j, i+1) = 3*sqrt(P0(1,1));
            sigma_y_vec(j, i+1) = 3*sqrt(P0(2,2));
            sigma_z_vec(j, i+1) = 3*sqrt(P0(3,3));
            sigma_r_vec(j, i+1) = 3*sqrt(trace(P0(1:3, 1:3)));
    
        else
            [xx_prop, P_prop] = UT(r_start, v_start, P0, et_vec(i), et_vec(i+1), C20, C22);
            
            r_start = xx_prop(1:3);
            v_start = xx_prop(4:6);
            P0 = P_prop;
            P_tutte(:, :, i+1) = P0;
    
            xx(i+1, :) = xx_filtered(end, :);
            sigma_r(i+1) = sqrt(trace(P0(1:3, 1:3)));
            sigma_v(i+1) = sqrt(trace(P0(4:6, 4:6)));
            P_tutte(:, :, i+1) = P0;
            sigma_x_vec(j, i) = 3*sqrt(P0(1,1));
            sigma_y_vec(j, i) = 3*sqrt(P0(2,2));
            sigma_z_vec(j, i) = 3*sqrt(P0(3,3));
            sigma_r_vec(j, i) = 3*sqrt(trace(P0(1:3, 1:3)));

        end
            
    
    end
    
    filter_error_pos = xx(:, 1:3) - xx_true(:, 1:3);
    
    filter_error_x_vec(j, :) = filter_error_pos(:, 1);
    filter_error_y_vec(j, :) = filter_error_pos(:, 2);
    filter_error_z_vec(j, :) = filter_error_pos(:, 3);

    filter_error_vec(j, :) = sqrt(sum(filter_error_pos.^2, 2));
    waitbar(j /(size(x_0_vec, 2)), h);
end
close(h);

%%

figure(1);

% --- Subplot 1 ---
subplot(2,2,1);
for i = 1:size(filter_error_vec, 1)
    err = plot(et_vec, abs(filter_error_vec(i, :)), 'r', 'LineWidth', 1.5);
    hold on;
    sig = plot(et_vec, sigma_r_vec(i, :), 'g', 'LineWidth', 1.5);
end
grid on; grid minor;
title('Position error');
xlabel('Time [s]'); ylabel('Error [km]');

% --- Subplot 2 ---
subplot(2,2,2);
for i = 1:size(filter_error_x_vec, 1)
    err = plot(et_vec, filter_error_x_vec(i, :), 'r', 'LineWidth', 1.5);
    hold on;
    sig = plot(et_vec, sigma_x_vec(i, :), 'g', 'LineWidth', 1.5);
    plot(et_vec, -sigma_x_vec(i, :), 'g', 'LineWidth',1.5);

end
grid on; grid minor;
title('rx error');
xlabel('Time [s]'); ylabel('Error [km]');
legend([sig, err], '3*sqrt(P_x)', 'norm of error on position component')

% --- Subplot 3 ---
subplot(2,2,3);
for i = 1:size(filter_error_y_vec, 1)
    err = plot(et_vec, filter_error_y_vec(i, :), 'r', 'LineWidth', 1.5);
    hold on;
    sig = plot(et_vec, sigma_y_vec(i, :), 'g', 'LineWidth', 1.5);
    plot(et_vec, -sigma_y_vec(i, :), 'g', 'LineWidth',1.5);
end
legend([sig, err], '3*sqrt(P_y)', 'norm of error on position component')
grid on; grid minor;
title('ry error');
xlabel('Time [s]'); ylabel('Error [km]');

% --- Subplot 4 ---
subplot(2,2,4);
for i = 1:size(filter_error_z_vec, 1)
    err = plot(et_vec, filter_error_z_vec(i, :), 'r', 'LineWidth', 1.5);
    hold on;
    sig = plot(et_vec, sigma_z_vec(i, :), 'g', 'LineWidth', 1.5);
    plot(et_vec, -sigma_z_vec(i, :), 'g', 'LineWidth',1.5);

end
legend([sig, err], '3*sqrt(P_z)', 'norm of error on position component')

grid on; grid minor;
title('rz error');
xlabel('Time [s]'); ylabel('Error [km]');
