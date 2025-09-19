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
centres = [310, 1700];
vals = [1, 2];

% Asteroid model generation
[F, V, N, C20, C22, A, score] = EllipsoidGenerationFib(a, b, c, M, centres, vals);

%Set initial asteroid knowledge
face_centroids = (V(F(:,1), :) + V(F(:,2), :) + V(F(:,3), :)) / 3;
known_map = double( ...
    face_centroids(:,1) > 0 & ...
    face_centroids(:,2) > 0 & ...
    face_centroids(:,3) < 0 );

known_map = zeros(size(F, 1), 1);
%Features for navigation
nav_index = round(linspace(1, length(known_map), 300));

% figure()
plotEllipsoidWithKnownRegion(F,V,score,known_map);
% figure()
%plotEllipsoidWithFeatures(F, V, ones(size(F, 1), 1), nav_index);

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
spacecraft_data.data_guidance.planningTime = 5*60;

%Set up data for the score function
spacecraft_data.data_guidance.scientificFov1 = 8*pi/180;
spacecraft_data.data_guidance.scientificFov2 = 8*pi/180;
spacecraft_data.data_guidance.navigationFov1 = 15*pi/180;
spacecraft_data.data_guidance.navigationFov2 = 15*pi/180;
spacecraft_data.data_guidance.alpha = alpha;
spacecraft_data.data_guidance.detPRef = -42;
spacecraft_data.data_guidance.DU = 40;
spacecraft_data.data_guidance.minLandmToScore = 5;

%Set up data for navigation and uncertainties
spacecraft_data.data_guidance.eta0 = eta0;
spacecraft_data.data_guidance.tao = 24*3600;

%Set up data for MCTS
spacecraft_data.data_guidance.ka = 3;
spacecraft_data.data_guidance.alpha_a = 0.1;
spacecraft_data.data_guidance.ko = 3;
spacecraft_data.data_guidance.alpha_o = 0.1;
spacecraft_data.data_guidance.gamma = 1;
spacecraft_data.data_guidance.expl_const = sqrt(2)/1.5;

%Set up data for asteroid state
spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.mass = 6.687e15;
spacecraft_data.data_asteroids.omega = 2*pi/(5.27025547*3600)*[0; 0; 1];
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;
spacecraft_data.data_asteroids.mapping.known_map = known_map; %known zones
spacecraft_data.data_asteroids.mapping.incidence = 85*pi/180; 
spacecraft_data.data_asteroids.mapping.emission = 85*pi/180;
spacecraft_data.data_asteroids.features.score = score;
spacecraft_data.data_asteroids.features.known_map_features = known_map;
spacecraft_data.data_asteroids.navigation_features = nav_index;

%SETUP 1: DETERMINISTIC
spacecraft_data.data_guidance.modelDynamics = @(t, x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22);
spacecraft_data.data_guidance.trueDynamics = @(t, x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22);
spacecraft_data.data_guidance.sigma_magn = 0;
spacecraft_data.data_guidance.sigma_align = 0;
spacecraft_data.data_guidance.process_noise = 0;
spacecraft_data.data_guidance.measurement_noise = 0;
sigma_pos = 0; %km
sigma_vel = 0; %km/s
sigma_eta = 0;
P0 = zeros(6);
P0(1:3, 1:3) = eye(3)*sigma_pos^2;
P0(4:6, 4:6) = eye(3)*sigma_vel^2;
%P0(7:9, 7:9) = eye(3)*sigma_eta^2;

%SETUP 2: NAVIGATION UNCERTAINTIES+CONTROL UNCERTAINTIES
spacecraft_data.data_guidance.modelDynamics = @(t, x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22);
spacecraft_data.data_guidance.trueDynamics = @(t, x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22);
spacecraft_data.data_guidance.sigma_magn =  0.05/3;
spacecraft_data.data_guidance.sigma_align = 2/3*pi/180;
spacecraft_data.data_guidance.process_noise = 0;
spacecraft_data.data_guidance.measurement_noise = sqrt( (100/3600*pi/180)^2 + (100/3600*pi/180)^2 ) ;
sigma_pos = 0.01; %km
sigma_vel = 0.01*10^-3; %km/s
sigma_eta = 10^(-12);
P0 = zeros(6);
P0(1:3, 1:3) = eye(3)*sigma_pos^2;
P0(4:6, 4:6) = eye(3)*sigma_vel^2;
%P0(7:9, 7:9) = eye(3)*sigma_eta^2;

%SETUP 3: DYN + CONTROL UNCERTAINTIES
spacecraft_data.data_guidance.modelDynamics = @(t, x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22);
spacecraft_data.data_guidance.trueDynamics = @(t, x) dynamicsTrue(t, x, mass_eros, omega_body);
spacecraft_data.data_guidance.sigma_magn = 0.05/3;
spacecraft_data.data_guidance.sigma_align = 2/3*pi/180;
spacecraft_data.data_guidance.process_noise = 5*10^-9;
spacecraft_data.data_guidance.measurement_noise = sqrt( (100/3600*pi/180)^2 + (100/3600*pi/180)^2 ) ;
sigma_pos = 0.01; %km
sigma_vel = 0.01*10^-3; %km/s
sigma_eta = 10^(-12);
P0 = zeros(6);
P0(1:3, 1:3) = eye(3)*sigma_pos^2;
P0(4:6, 4:6) = eye(3)*sigma_vel^2;
% %P0(7:9, 7:9) = eye(3)*sigma_eta^2;

%% RUN MCTS

load('sigma_magn_precomputed.mat');
load('sigma_align_precomputed.mat');
load('sigma_meas_precomputed.mat');
load('theta_precomputed.mat');
n_run = 100; %Max 100

for i = 1:n_run
    spacecraft_data.MC.pert_magn = sigma_magn_pre(i, :);
    spacecraft_data.MC.pert_align = sigma_align_pre(i, :);
    spacecraft_data.MC.pert_theta = theta_pre(i, :);
    spacecraft_data.MC.iter_man = 1;
    spacecraft_data.MC.pert_meas = sigma_meas_pre(i, :);
    spacecraft_data.MC.iter_meas = 1;

end
iterations = 250; %Number of iterations per tree (As a reference, for 200 iterations 45 min/1 h are required)
n_trees = 4; %Number of trees

profile clear
profile on
[all_trees, real_trajectory_MCTS, filter_trajectory_MCTS, ref_trajectory_MCTS, P_all_MCTS, real_times_MCTS, ...
 mapping_scores_MCTS, exploiting_scores_MCTS, nav_scores_MCTS, action_times_MCTS, all_flag_MCTS, planned_actions_MCTS, ...
 real_trajectory_total_MCTS, ref_trajectory_total_MCTS, filter_trajectory_total_MCTS, tt_all_MCTS, ...
 mapping_score_total_MCTS, exploiting_score_total_MCTS, nav_score_total_MCTS, P_total_MCTS, spacecraft_data_out_MCTS, correction1, correction2, actual_maneuvres]= ...
    runMCTSBatch(spacecraft_data, r0, v0, t0, P0, iterations, n_trees, options);
profile off
profile viewer

%%
% [real_trajectory_MCTS, filter_trajectory_MCTS, P_all_MCTS, tt_all_MCTS, ...
%           total_mapping_score_MCTS, total_exploiting_score_MCTS, total_nav_score_MCTS, ...
%           action_times_MCTS, all_flag_MCTS, action_list_MCTS, spacecraft_data_out_MCTS] = ...
%     recoverMCTSBatch(spacecraft_data, r0, v0, t0, P0, all_trees, n_trees, options);

%% GREEDY
n_set = size(planned_actions_MCTS, 1);
n_actions = zeros(n_set, 1);
for i = 1:size(planned_actions_MCTS, 1)
    for j = 1:size(planned_actions_MCTS, 2)
        if ~isempty(planned_actions_MCTS{i, j})
            n_actions(i) = n_actions(i)+1;
        end
    end
end

[real_trajectory_GREEDY, filter_trajectory_GREEDY, ref_trajectory_GREEDY, P_all_GREEDY, real_times_GREEDY, ...
 mapping_scores_GREEDY, exploiting_scores_GREEDY, nav_scores_GREEDY, action_times_GREEDY, all_flag_GREEDY, planned_actions_GREEDY, ...
 real_trajectory_total_GREEDY, ref_trajectory_total_GREEDY, filter_trajectory_total_GREEDY, tt_all_GREEDY, ...
 mapping_score_total_GREEDY, exploiting_score_total_GREEDY, nav_score_total_GREEDY, P_total_GREEDY, spacecraft_data_out_GREEDY] = ...
    greedyApproach(spacecraft_data, r0, v0, t0, P0, n_actions, n_set, options);

%% Plots
%Plot errori posizione
final_scores_MCTS = spacecraft_data_out_MCTS.data_asteroids.features.score;
final_known_map_MCTS = spacecraft_data_out_MCTS.data_asteroids.mapping.known_map;
final_scores_GREEDY = spacecraft_data_out_GREEDY.data_asteroids.features.score;
final_known_map_GREEDY = spacecraft_data_out_GREEDY.data_asteroids.mapping.known_map;

action_times_list_MCTS = [];
action_times_list_GREEDY = [];
for i = 1:size(action_times_GREEDY, 1)
    action_times_list_MCTS = [action_times_list_MCTS, cell2mat(action_times_MCTS(i, :))];
    action_times_list_GREEDY = [action_times_list_GREEDY, cell2mat(action_times_GREEDY(i, :))];
end
action_times_list_GREEDY = action_times_list_GREEDY(action_times_list_GREEDY>0);
action_times_list_MCTS = action_times_list_MCTS(action_times_list_MCTS>0);

figure(1)
subplot(1,2,1)
tresig = [];
error = vecnorm(real_trajectory_total_MCTS(:, 1:3)-filter_trajectory_total_MCTS(:, 1:3), 2, 2);
error_plot = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, error, 'b', 'LineWidth',1.5);
hold on
for i = 1:size(P_total_MCTS, 3)
    tresig = [tresig, 3*sqrt(trace(P_total_MCTS(1:3, 1:3, i)))];
end
for i = 1:length(action_times_list_MCTS)
    act = xline( (action_times_list_MCTS(i)-tt_all_MCTS(1) )/3600, 'k--', 'LineWidth', 1);
end
tresig_plot = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, tresig, 'g', 'LineWidth',1.5);
plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, -tresig, 'g', 'LineWidth',1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('[km]')
legend([tresig_plot, error_plot], '3sigma r', 'error r')
title('MCTS')

subplot(1,2,2)
tresig = [];
error = vecnorm(real_trajectory_total_GREEDY(:, 1:3)-filter_trajectory_total_GREEDY(:, 1:3), 2, 2);
error_plot = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, error, 'b', 'LineWidth',1.5);
hold on
for i = 1:size(P_total_GREEDY, 3)
    tresig = [tresig, 3*sqrt(trace(P_total_GREEDY(1:3, 1:3, i)))];
end
for i = 1:length(action_times_list_GREEDY)
    act = xline( (action_times_list_GREEDY(i)-tt_all_GREEDY(1) )/3600, 'k--', 'LineWidth', 1);
end
tresig_plot = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, tresig, 'g', 'LineWidth',1.5);
plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, -tresig, 'g', 'LineWidth',1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('[km]')
legend([tresig_plot, error_plot], '3sigma r', 'error r')
title('GREEDY')


%Plot errori velocità
figure(2)
subplot(1,2,1)
tresigv = [];
errorv = vecnorm(real_trajectory_total_MCTS(:, 4:6)-filter_trajectory_total_MCTS(:, 4:6), 2, 2);
error_plot = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, errorv*10^3, 'b', 'LineWidth',1.5);
hold on
for i = 1:size(P_total_MCTS, 3)
    tresigv = [tresigv, 3*sqrt(trace(P_total_MCTS(4:6, 4:6, i)))];
end
for i = 1:length(action_times_list_MCTS)
    act = xline( (action_times_list_MCTS(i)-tt_all_MCTS(1) )/3600, 'k--', 'LineWidth', 1);
end
tresig_plot = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, tresigv*10^3, 'g', 'LineWidth',1.5);
plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, -tresigv*10^3, 'g', 'LineWidth',1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('[m/s]')
legend([tresig_plot, error_plot], '3sigma v', 'error v')
title('MCTS')


subplot(1,2,2)
tresigv = [];
errorv = vecnorm(real_trajectory_total_GREEDY(:, 4:6)-filter_trajectory_total_GREEDY(:, 4:6), 2, 2);
error_plot = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, errorv*10^3, 'b', 'LineWidth',1.5);
hold on
for i = 1:size(P_total_GREEDY, 3)
    tresigv = [tresigv, 3*sqrt(trace(P_total_GREEDY(4:6, 4:6, i)))];
end
for i = 1:length(action_times_list_GREEDY)
    act = xline( (action_times_list_GREEDY(i)-tt_all_GREEDY(1) )/3600, 'k--', 'LineWidth', 1);
end
tresig_plot = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, tresigv*10^3, 'g', 'LineWidth',1.5);
plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, -tresigv*10^3, 'g', 'LineWidth',1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('[m/s]')
legend([tresig_plot, error_plot], '3sigma v', 'error v')
title('GREEDY')


%Plot Score
m_MCTS = cumsum(mapping_score_total_MCTS*alpha(1));
e_MCTS = cumsum(exploiting_score_total_MCTS*alpha(2));
n_MCTS = cumsum(nav_score_total_MCTS*alpha(3));
m_GREEDY = cumsum(mapping_score_total_GREEDY*alpha(1));
e_GREEDY = cumsum(exploiting_score_total_GREEDY*alpha(2));
n_GREEDY = cumsum(nav_score_total_GREEDY*alpha(3));
ymin_all = min([m_MCTS(:); e_MCTS(:); n_MCTS(:); m_GREEDY(:); e_GREEDY(:); n_GREEDY(:)]);
ymax_all = max([m_MCTS(:); e_MCTS(:); n_MCTS(:); m_GREEDY(:); e_GREEDY(:); n_GREEDY(:)]);

figure(3)
subplot(1,2,1)
m = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, m_MCTS, 'g', 'LineWidth', 1.5);
hold on
e = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, e_MCTS, 'b', 'LineWidth', 1.5);
n = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, n_MCTS, 'r', 'LineWidth', 1.5);
for i = 1:length(action_times_list_MCTS)
    act = xline( (action_times_list_MCTS(i)-tt_all_MCTS(1) )/3600, 'k--', 'LineWidth', 1);
end
grid on
grid minor
xlabel('Time [h]')
ylabel('Score [-]')
legend([m, e, n], 'Mapping Score*alpha1', 'Exploit Score*alpha2', 'Navigation Score*alpha3')
title('MCTS')

subplot(1,2,2)
m = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, m_GREEDY, 'g', 'LineWidth', 1.5);
hold on
e = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, e_GREEDY, 'b', 'LineWidth', 1.5);
n = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, n_GREEDY, 'r', 'LineWidth', 1.5);
for i = 1:length(action_times_list_GREEDY)
    act = xline( (action_times_list_GREEDY(i)-tt_all_GREEDY(1) )/3600, 'k--', 'LineWidth', 1);
end
grid on
grid minor
xlabel('Time [h]')
ylabel('Score [-]')
legend([m, e, n], 'Mapping Score*alpha1', 'Exploit Score*alpha2', 'Navigation Score*alpha3')
title('GREEDY')

figure(12)
subplot(1,2,1)
m = plot(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, m_MCTS+e_MCTS+n_MCTS, 'g', 'LineWidth', 1.5);
hold on
for i = 1:length(action_times_list_MCTS)
    act = xline( (action_times_list_MCTS(i)-tt_all_MCTS(1) )/3600, 'k--', 'LineWidth', 1);
end
grid on
grid minor
xlabel('Time [h]')
ylabel('Score [-]')
title('MCTS')

subplot(1,2,2)
m = plot(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, m_GREEDY+e_GREEDY+n_GREEDY, 'g', 'LineWidth', 1.5);
hold on
for i = 1:length(action_times_list_GREEDY)
    act = xline( (action_times_list_GREEDY(i)-tt_all_GREEDY(1) )/3600, 'k--', 'LineWidth', 1);
end
grid on
grid minor
xlabel('Time [h]')
ylabel('Score [-]')
title('GREEDY')

%Plot traiettoria in body frame
figure(4)
subplot(1,2,1)
plot3(real_trajectory_total_MCTS(:, 1), real_trajectory_total_MCTS(:, 2), real_trajectory_total_MCTS(:, 3), 'b', 'LineWidth', 1.5)
hold on
for i = 1:length(action_times_list_MCTS)
    pos = find(abs(tt_all_MCTS - action_times_list_MCTS(i)) < 1e-9, 1, 'first');
    man_MCTS = plot3(real_trajectory_total_MCTS(pos, 1), real_trajectory_total_MCTS(pos, 2), real_trajectory_total_MCTS(pos, 3), 'b.', 'MarkerSize', 15);
end
plotEllipsoidWithKnownRegion(F, V, final_scores_MCTS, final_known_map_MCTS)
grid on
grid minor
xlabel('X [km]')
zlabel('Z [km]')
ylabel('Y [km]')
legend(man_MCTS, 'Manoeuvring Point')
title('MCTS')

subplot(1,2,2)
plot3(real_trajectory_total_GREEDY(:, 1), real_trajectory_total_GREEDY(:, 2), real_trajectory_total_GREEDY(:, 3), 'b', 'LineWidth', 1.5)
hold on
for i = 1:length(action_times_list_GREEDY)
    pos = find(abs(tt_all_GREEDY - action_times_list_GREEDY(i)) < 1e-9, 1, 'first');
    man_GREEDY = plot3(real_trajectory_total_GREEDY(pos, 1), real_trajectory_total_GREEDY(pos, 2), real_trajectory_total_GREEDY(pos, 3), 'b.', 'MarkerSize', 15);
end
plotEllipsoidWithKnownRegion(F, V, final_scores_GREEDY, final_known_map_GREEDY)
grid on
grid minor
xlabel('X [km]')
zlabel('Z [km]')
ylabel('Y [km]')
legend(man_GREEDY, 'Manoeuvring Point')
title('GREEDY')

%Plot traiettoria in inertial
figure(5)
trajectory_in_MCTS = zeros(size(real_trajectory_total_MCTS(:, 1:3)));
for i = 1:length(tt_all_MCTS)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_all_MCTS(i)); 
    trajectory_in_MCTS(i, :) = (R_body2in*real_trajectory_total_MCTS(i, 1:3)')';
end
trajectory_in_GREEDY = zeros(size(real_trajectory_total_GREEDY(:, 1:3)));
for i = 1:length(tt_all_GREEDY)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_all_GREEDY(i)); 
    trajectory_in_GREEDY(i, :) = (R_body2in*real_trajectory_total_GREEDY(i, 1:3)')';
end

subplot(1,2,1)
plotEllipsoidWithKnownRegion(F, V, final_scores_MCTS, final_known_map_MCTS);
plot3(trajectory_in_MCTS(:, 1), trajectory_in_MCTS(:, 2), trajectory_in_MCTS(:, 3), 'b', 'LineWidth', 1.5)
hold on
for i = 1:length(action_times_list_MCTS)
    pos = find(abs(tt_all_MCTS - action_times_list_MCTS(i)) < 1e-9, 1, 'first');
    man_MCTS = plot3(trajectory_in_MCTS(pos, 1), trajectory_in_MCTS(pos, 2), trajectory_in_MCTS(pos, 3), 'b.', 'MarkerSize', 15);
end
title('Trajectory in ECLIPJ2000 - MCTS')
grid on
grid minor
xlabel('X [km]')
zlabel('Z [km]')
ylabel('Y [km]')
legend(man_MCTS, 'Manoeuvring Point')

subplot(1,2,2)
plotEllipsoidWithKnownRegion(F, V, final_scores_GREEDY, final_known_map_GREEDY);
plot3(trajectory_in_GREEDY(:, 1), trajectory_in_GREEDY(:, 2), trajectory_in_GREEDY(:, 3), 'b', 'LineWidth', 1.5)
hold on
for i = 1:length(action_times_list_GREEDY)
    pos = find(abs(tt_all_GREEDY - action_times_list_GREEDY(i)) < 1e-9, 1, 'first');
    man_GREEDY = plot3(trajectory_in_GREEDY(pos, 1), trajectory_in_GREEDY(pos, 2), trajectory_in_GREEDY(pos, 3), 'b.', 'MarkerSize', 15);
end
title('Trajectory in ECLIPJ2000 - GREEDY')
grid on
grid minor
xlabel('X [km]')
zlabel('Z [km]')
ylabel('Y [km]')
legend(man_GREEDY, 'Manoeuvring Point')

%Plot errori velocità
% figure(6)
% subplot(1,2,1)
% determs = [];
% hold on
% for i = 1:size(P_all_MCTS, 3)
%     DU = 40;
%     TU = sqrt( DU^3/(astroConstants(1)*spacecraft_data.data_asteroids.mass) );
%     VU = DU/TU + omega_body(3)*DU;
% 
%     P_adim = zeros(size(P_all_MCTS(1:6, 1:6)));
%     P_adim(1:3, 1:3) = P_all_MCTS(1:3, 1:3, i)/(DU*DU);
%     P_adim(4:6, 1:3) = P_all_MCTS(4:6, 1:3, i)/(DU*VU);
%     P_adim(1:3, 4:6) = P_all_MCTS(1:3, 4:6, i)/(DU*VU);
%     P_adim(4:6, 4:6) = P_all_MCTS(4:6, 4:6, i)/(VU*VU);
%     determs = [determs, log10(det(P_adim))];
% end
% for i = 1:length(action_times_MCTS)
%     act = xline( (action_times_MCTS(i)-tt_all_MCTS(1) )/3600, 'k--', 'LineWidth', 1);
% end
% deter_plot = semilogy(( tt_all_MCTS(1:end)-tt_all_MCTS(1) )/3600, determs, 'g', 'LineWidth',1.5);
% grid on
% grid minor
% xlabel('Time [h]')
% ylabel('[-]')
% legend([deter_plot, act], 'log10(det(P)) [-]', 'action')
% title('MCTS')
% 
% subplot(1,2,2)
% determs = [];
% hold on
% for i = 1:size(P_all_GREEDY, 3)
%     DU = 40;
%     TU = sqrt( DU^3/(astroConstants(1)*spacecraft_data.data_asteroids.mass) );
%     VU = DU/TU + omega_body(3)*DU;
% 
%     P_adim = zeros(size(P_all_GREEDY(1:6, 1:6)));
%     P_adim(1:3, 1:3) = P_all_GREEDY(1:3, 1:3, i)/(DU*DU);
%     P_adim(4:6, 1:3) = P_all_GREEDY(4:6, 1:3, i)/(DU*VU);
%     P_adim(1:3, 4:6) = P_all_GREEDY(1:3, 4:6, i)/(DU*VU);
%     P_adim(4:6, 4:6) = P_all_GREEDY(4:6, 4:6, i)/(VU*VU);
%     determs = [determs, log10(det(P_adim))];
% end
% for i = 1:length(action_times_GREEDY)
%     act = xline( (action_times_GREEDY(i)-tt_all_GREEDY(1) )/3600, 'k--', 'LineWidth', 1);
% end
% deter_plot = semilogy(( tt_all_GREEDY(1:end)-tt_all_GREEDY(1) )/3600, determs, 'g', 'LineWidth',1.5);
% grid on
% grid minor
% xlabel('Time [h]')
% ylabel('[-]')
% legend([deter_plot, act], 'log10(det(P)) [-]', 'action')
% title('GREEDY')

%Confronto dispersion e knowledge
figure(8)
dispersion_posMCTS = vecnorm(ref_trajectory_total_MCTS(:, 1:3)-real_trajectory_total_MCTS(:, 1:3), 2, 2);
knowledge_posMCTS = vecnorm(filter_trajectory_total_MCTS(:, 1:3)-real_trajectory_total_MCTS(:, 1:3), 2, 2);

dispe = plot((tt_all_MCTS-tt_all_MCTS(1))/3600, dispersion_posMCTS, 'r', 'LineWidth', 1.5)
hold on
knowl = plot((tt_all_MCTS-tt_all_MCTS(1))/3600, knowledge_posMCTS, 'b', 'LineWidth', 1.5)
grid on
grid minor
xlabel('Time [h]')
ylabel('Val [Km]')
legend([dispe, knowl], 'Dispersion', 'Knowledge')
title('Dispersion vs Knowledge MCTS')

magn_maneuvres = {};
magn_corrections1 = {};
magn_corrections2 = {};
for i = 1:size(actual_maneuvres, 1)
    for j = 1:size(actual_maneuvres, 2)
        magn_maneuvres{i,j} = norm(actual_maneuvres{i,j});
        magn_corrections1{i,j} = norm(correction1{i, j});
        magn_corrections2{i,j} = norm(correction2{i, j});
        
    end
end

%Arc duration
arc_duration_MCTS = [];
arc_duration_GREEDY = [];
for i = 1:length(action_times_list_MCTS)-1
    arc_duration_MCTS(i) = (action_times_list_MCTS(i+1)-action_times_list_MCTS(i))/3600;
    arc_duration_GREEDY(i) = (action_times_list_GREEDY(i+1)-action_times_list_GREEDY(i))/3600;
end
n_arcs = length(arc_duration_MCTS);
x = 1:n_arcs;
Y = [arc_duration_MCTS(:), arc_duration_GREEDY(:)];

figure(9)
bar(x, Y)  % colonnine affiancate
xlabel('Numero arco');
ylabel('Durata arco [ore]');
title('Durata degli archi MCTS vs GREEDY');
legend('MCTS', 'GREEDY');
grid on
grid minor
yline((spacecraft_data.data_guidance.Th_max-spacecraft_data.data_guidance.safety_margin)/3600, '--', 'LineWidth', 1.5)
yline(spacecraft_data.data_guidance.DeltaT_after_man/3600, '--', 'LineWidth', 1.5)
