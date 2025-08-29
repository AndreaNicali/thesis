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

%% Initial conditions 

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

x_mean = [21.095697268161683; -31.527866008753950; -23.491639672861801;  -0.011418159898777;  -0.009247442529600;  0.001997382918943];

%Initial Covariance
sigma_pos = 0.5e-1;
sigma_vel = 0.3*10^(-4);
P_initial = [eye(3)*sigma_pos^2, zeros(3); zeros(3),  eye(3)*sigma_vel^2];
n_sample = 30;
x_0_vec = zeros(6, n_sample);

for i = 1:size(x_0_vec, 2)
    x_0_vec(:, i) = mvnrnd(x_mean, P_initial) ;
end

%Epochs
et_i = cspice_str2et('2000-08-01 T07:00:00');
et_f = cspice_str2et('2000-08-01 T19:00:00');
et_step = 60;
et_vec = et_i:et_step:et_f;

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

known_map = ones(size(F, 1), 1);

features.nav_index = round(linspace(1, length(known_map), 300));

% FOV of the Navigation camera
fov1 = 15*pi/180;
fov2 = 15*pi/180;

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
spacecraft_data.data_asteroids.features.known_map_features = known_map;
spacecraft_data.data_asteroids.navigation = features.nav_index;

%Iniatilize measurements struct
measurements = struct();

%Measurements noise
sigma1 = 100; %37.2
sigma2 = 100; %90
sigma_angular = sqrt( (sigma1/3600*pi/180)^2 + (sigma2/3600*pi/180)^2 ) ;

%Process noise
sigma_acc = 5e-9;

%Integration of truth and model dynamics in the given time interval
[~ , xx_model] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), et_vec, x_mean, options);

%Measurement noise
sigma_meas = sigma_angular;

%Initialize storing vectors

P_tutte = zeros(6,6, length(et_vec));
P_tutte(:, :, 1) = P_initial;

error_r = zeros(size(x_0_vec, 2), length(et_vec)); 
error_v = zeros(size(x_0_vec, 2), length(et_vec)); 
error_x = zeros(size(x_0_vec, 2), length(et_vec)); 
error_y = zeros(size(x_0_vec, 2), length(et_vec)); 
error_z = zeros(size(x_0_vec, 2), length(et_vec)); 

sigma_r_vec = zeros(size(x_0_vec, 2), length(et_vec)); 
sigma_v_vec = zeros(size(x_0_vec, 2), length(et_vec)); 
sigma_x_vec = zeros(size(x_0_vec, 2), length(et_vec)); 
sigma_y_vec = zeros(size(x_0_vec, 2), length(et_vec)); 
sigma_z_vec = zeros(size(x_0_vec, 2), length(et_vec)); 


h = waitbar(0, 'Propagating Uncertainties...'); 

for j = 1:size(x_0_vec, 2)
    
    x_0 = x_0_vec(:, j);

    r_start = x_0_vec(1:3, j)';
    v_start = x_0_vec(4:6, j)';
    xx(1, :) = x_0_vec(1:6, j);
    P0 = P_initial;

     [eval_times , xx_true] = ode78(@(t,x) dynamics15(t, x, mass_eros, omega_body), et_vec, x_0, options);
     %[~ , xx_true] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), et_vec, x_0, options);
     %sigma_acc = 0;

    %Obtain measurements from true trajectory
    [measurements] = measurementsFun(xx_true, eval_times, spacecraft_data);
    pert_meas = zeros(size(measurements.val));
    
    %Perturb measurements 
    for i = 1:size(measurements.val, 2)
        vals = measurements.val(:, i);
        elevation = mvnrnd(measurements.val(1,i), sigma_meas^2);
        azimuth = mvnrnd(measurements.val(2,i), sigma_meas^2);
        pert_meas(:, i) = [elevation;azimuth];
    end
    measurements.val = pert_meas;

    [xx_filtered, P_filtered, flag] = EKF([r_start, v_start], eval_times, measurements, P0, spacecraft_data, sigma_meas, sigma_acc, xx_true);
    
    error_r(j, :) = vecnorm(xx_filtered(:, 1:3)-xx_true(:, 1:3), 2, 2);
    error_v(j, :) = vecnorm(xx_filtered(:, 4:6)-xx_true(:, 4:6), 2, 2);
    error_x(j, :) = vecnorm(xx_filtered(:, 1)-xx_true(:, 1), 2, 2);
    error_y(j, :) = vecnorm(xx_filtered(:, 2)-xx_true(:, 2), 2, 2);
    error_z(j, :) = vecnorm(xx_filtered(:, 3)-xx_true(:, 3), 2, 2);

    for i = 1:size(P_filtered, 3)
        sigma_r_vec(j, i) = 3*sqrt(trace(P_filtered(1:3, 1:3, i)));
        sigma_v_vec(j, i) = 3*sqrt(trace(P_filtered(4:6, 4:6, i)));
        sigma_x_vec(j, i) = 3*sqrt(trace(P_filtered(1, 1, i)));
        sigma_y_vec(j, i) = 3*sqrt(trace(P_filtered(2, 2, i)));
        sigma_z_vec(j, i) = 3*sqrt(trace(P_filtered(3, 3, i)));
    end

    waitbar(j /(size(x_0_vec, 2)), h);
 
end
close(h);

figure(1)
hold on
grid on
grid minor
for i = 1:size(x_0_vec, 2)
    plot(et_vec, sigma_r_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, -sigma_r_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, error_r(i, :), 'b', 'LineWidth', 1)
end
    
figure(2)
hold on
grid on
grid minor
for i = 1:size(x_0_vec, 2)
    plot(et_vec, sigma_v_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, -sigma_v_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, error_v(i, :), 'b', 'LineWidth', 1)
end

figure(3)
hold on
grid on
grid minor
for i = 1:size(x_0_vec, 2)
    plot(et_vec, sigma_x_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, -sigma_x_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, error_x(i, :), 'b', 'LineWidth', 1)
end

figure(4)
hold on
grid on
grid minor
for i = 1:size(x_0_vec, 2)
    plot(et_vec, sigma_y_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, -sigma_y_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, error_y(i, :), 'b', 'LineWidth', 1)
end

figure(5)
hold on
grid on
grid minor
for i = 1:size(x_0_vec, 2)
    plot(et_vec, sigma_z_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, -sigma_z_vec(i, :), 'g', 'LineWidth', 1.5)
    plot(et_vec, error_z(i, :), 'b', 'LineWidth', 1)
end