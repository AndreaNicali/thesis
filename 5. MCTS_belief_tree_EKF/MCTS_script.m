clc
clear
close all

addpath(genpath('kernels\'))
addpath(genpath('functions\'))
addpath(genpath('3dModel\'))
addpath(genpath('dynamics\'))
addpath(genpath('Uncertainties\'))
addpath(genpath('Reachability\'))
addpath(genpath('MCTS functions\'))

cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\erosephem_1999004_2002181.bsp')
cspice_furnsh('kernels\erosatt_1998329_2001157_v01.bpc');

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

%% Initial data
%Initial Conditions
t0 = cspice_str2et('2000-08-01 T01:00:00');

r0 = [8.5855;   44.6428;   -4.4817]; %km
v0 = [0.0164;   -0.0032;    0.0018]; %km/s
eta0 = zeros(1, 3);

sigma_pos = 0.01; %km
sigma_vel = 0.01*10^-3; %km/s
sigma_eta = 10^(-12);
P0 = zeros(9);
P0(1:3, 1:3) = eye(3)*sigma_pos^2;
P0(4:6, 4:6) = eye(3)*sigma_vel^2;
P0(7:9, 7:9) = eye(3)*sigma_eta^2;

%Data for body frame propagation
T_rotation_eros = 5.27025547*3600; %s
omega_body = 2*pi/T_rotation_eros*[0; 0; 1]; %rad/s
mass_eros = 6.687e15; %kg

% Semiaxis (in km) of the asteroid model 
a = 20.591;
b = 5.711;
c = 5.332;
M = 1000;

% Features positioning and values
centres = [120, 965];
vals = [1, 2];

% Asteroid model generation
[F, V, N, C20, C22, A, score] = EllipsoidGenerationFib(a, b, c, M, centres, vals);

%Set initial asteroid knowledge
% face_centroids = (V(F(:,1), :) + V(F(:,2), :) + V(F(:,3), :)) / 3;
% known_map = double( ...
%     face_centroids(:,1) > 0 & ...
%     face_centroids(:,2) > 0 & ...
%     face_centroids(:,3) > 0 );

known_map = zeros(size(F, 1), 1);
%Features for navigation
nav_index = round(linspace(1, length(known_map), 300));

% figure()
% plotEllipsoidWithKnownRegion(F,V,score,ones(size(F, 1), 1));
% figure()
% plotEllipsoidWithFeatures(F, V, ones(size(F, 1), 1), nav_index);

alpha = [1/3; 1/3; 1/3]; %[mapping; exploitation; navigation]
%Set up data for reachability
spacecraft_data = struct();
spacecraft_data.data_guidance.ReachabilityScoreComputation = 1;
spacecraft_data.data_guidance.ReachabilityExplorationScheme = 2;
spacecraft_data.data_guidance.DeltaV_max = 2e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 2*3600;
spacecraft_data.data_guidance.DeltaT_after_man = 0.3*3600;
spacecraft_data.data_guidance.r_impact = 24;
spacecraft_data.data_guidance.r_escape = 100;
spacecraft_data.data_guidance.scientificFov1 = 8*pi/180;
spacecraft_data.data_guidance.scientificFov2 = 8*pi/180;
spacecraft_data.data_guidance.navigationFov1 = 15*pi/180;
spacecraft_data.data_guidance.navigationFov2 = 15*pi/180;
spacecraft_data.data_guidance.alpha = alpha;
spacecraft_data.data_guidance.eta0 = eta0;
spacecraft_data.data_guidance.sigma_magn = 0.05/3;
spacecraft_data.data_guidance.sigma_align = 2/3*pi/180;

spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.mass = 6.687e15;
spacecraft_data.data_asteroids.omega = 2*pi/(5.27025547*3600)*[0; 0; 1];
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;
spacecraft_data.data_asteroids.mapping.known_map = known_map; %known zones
spacecraft_data.data_asteroids.mapping.incidence = [0, 85]*pi/180; 
spacecraft_data.data_asteroids.mapping.emission = [0, 85]*pi/180;
spacecraft_data.data_asteroids.features.score = score;
spacecraft_data.data_asteroids.features.known_map_features = known_map;
spacecraft_data.data_asteroids.navigation_features = nav_index;

model_dyn = @(t, x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22);
truth_dyn = @(t, x) dynamicsTrue(t, x, mass_eros, omega_body);

%% RUN MCTS

iterations = 50; %Number of iterations per tree (As a reference, for 200 iterations 45 min/1 h are required)
n_trees = 4; %Number of trees

[all_trees, real_trajectory, filter_trajectory, P_all, tt_all, ...
          total_mapping_score, total_exploiting_score, total_nav_score, ...
          action_times, all_flag, spacecraft_data_out] = ...
    runMCTSBatch(spacecraft_data, r0, v0, t0, P0, iterations, n_trees, truth_dyn, options);

%% Plots
%Plot errori posizione
final_scores = spacecraft_data_out.data_asteroids.features.score;
final_known_map = spacecraft_data_out.data_asteroids.mapping.known_map;

figure(1)
tresig = [];
error = vecnorm(real_trajectory(:, 1:3)-filter_trajectory(:, 1:3), 2, 2);
error_plot = plot(( tt_all(1:end)-tt_all(1) )/3600, error, 'b', 'LineWidth',1.5);
hold on
for i = 1:size(P_all, 3)
    tresig = [tresig, 3*sqrt(trace(P_all(1:3, 1:3, i)))];
    % if all_flag(i) == 1
    %     xline(( tt_all(i)-tt_all(1) )/3600, 'Color', [207 233 255]/255, 'LineWidth', 0.5);
    % end
end
for i = 1:length(action_times)
    act = xline( (action_times(i)-tt_all(1) )/3600, 'k--', 'LineWidth', 1);
end
tresig_plot = plot(( tt_all(1:end)-tt_all(1) )/3600, tresig, 'g', 'LineWidth',1.5);
plot(( tt_all(1:end)-tt_all(1) )/3600, -tresig, 'g', 'LineWidth',1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('[km]')
legend([tresig_plot, error_plot, act], '3sigma r', 'error r', 'action')

%Plot errori velocità
figure(2)
tresigv = [];
errorv = vecnorm(real_trajectory(:, 4:6)-filter_trajectory(:, 4:6), 2, 2);
error_plot = plot(( tt_all(1:end)-tt_all(1) )/3600, errorv*10^3, 'b', 'LineWidth',1.5);
hold on
for i = 1:size(P_all, 3)
    tresigv = [tresigv, 3*sqrt(trace(P_all(4:6, 4:6, i)))];
    % if all_flag(i) == 1
    %     xline(( tt_all(i)-tt_all(1) )/3600, 'Color', [207 233 255]/255, 'LineWidth', 0.5);
    % end
end
for i = 1:length(action_times)
    act = xline( (action_times(i)-tt_all(1) )/3600, 'k--', 'LineWidth', 1);
end
tresig_plot = plot(( tt_all(1:end)-tt_all(1) )/3600, tresigv*10^3, 'g', 'LineWidth',1.5);
plot(( tt_all(1:end)-tt_all(1) )/3600, -tresigv*10^3, 'g', 'LineWidth',1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('[m/s]')
legend([tresig_plot, error_plot, act], '3sigma v', 'error v', 'action')

%Plot Score
figure(3)
m = plot(( tt_all(1:end)-tt_all(1) )/3600, cumsum(total_mapping_score*alpha(1)), 'g', 'LineWidth', 1.5);
hold on
e = plot(( tt_all(1:end)-tt_all(1) )/3600, cumsum(total_exploiting_score*alpha(2)), 'b', 'LineWidth', 1.5);
n = plot(( tt_all(1:end)-tt_all(1) )/3600, cumsum(total_nav_score*alpha(3)), 'r', 'LineWidth', 1.5);
for i = 1:length(action_times)
    act = xline( (action_times(i)-tt_all(1) )/3600, 'k--', 'LineWidth', 1);
end
grid on
grid minor
xlabel('Time [h]')
ylabel('Score [-]')
legend([m, e, n, act], 'Mapping Score*alpha1', 'Exploit Score*alpha2', 'Navigation Score*alpha3', 'Actions')

%Plot traiettoria in body frame
figure(4)
plot3(real_trajectory(:, 1), real_trajectory(:, 2), real_trajectory(:, 3), 'b', 'LineWidth', 1.5)
hold on
for i = 1:length(action_times)
    pos = find(tt_all == action_times(i));
    man = plot3(real_trajectory(pos, 1), real_trajectory(pos, 2), real_trajectory(pos, 3), 'b.', 'MarkerSize', 15);
end
plotEllipsoidWithKnownRegion(F, V, final_scores, final_known_map)
grid on
grid minor
xlabel('X [km]')
zlabel('Z [km]')
ylabel('Y [km]')
legend(man, 'Manoeuvring Point')

%Plot traiettoria in inertial
figure(5)
trajectory_in = zeros(size(real_trajectory(:, 1:3)));
for i = 1:length(tt_all)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_all(i)); 
    trajectory_in(i, :) = (R_body2in*real_trajectory(i, 1:3)')';
end
plotEllipsoidWithKnownRegion(F, V, score, final_known_map);
plot3(trajectory_in(:, 1), trajectory_in(:, 2), trajectory_in(:, 3), 'b', 'LineWidth', 1.5)
hold on
for i = 1:length(action_times)
    pos = find(tt_all == action_times(i));
    man = plot3(trajectory_in(pos, 1), trajectory_in(pos, 2), trajectory_in(pos, 3), 'b.', 'MarkerSize', 15);
end
%plotSunPointingVector(tt_all, 'ECLIPJ2000')
title('Trajectory in ECLIPJ2000')
grid on
grid minor
xlabel('X [km]')
zlabel('Z [km]')
ylabel('Y [km]')
legend(man, 'Manoeuvring Point')

