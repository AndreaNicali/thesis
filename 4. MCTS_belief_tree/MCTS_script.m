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

%% Perform an iteration of reachability

%Epoch
t0 = cspice_str2et('2000-08-01 T01:00:00');
tf = cspice_str2et('2000-08-01 T11:00:00');
tt = t0:100:tf;

%Initial Condition
r0 = [8.5855;   44.6428;   -4.4817];
v0 = [0.0164;   -0.0032;    0.0018];
sigma_pos = 0.01;
sigma_vel = 0.03*10^-3;
P0 = [eye(3)*sigma_pos^2, zeros(3); zeros(3),  eye(3)*sigma_vel^2];

%Data for body frame propagation
T_rotation_eros = 5.27025547*3600; 
omega_body = 2*pi/T_rotation_eros*[0; 0; 1];
mass_eros = 6.687e15;

% Semiaxis (in km) of the asteroid model 
a = 20.591;
b = 5.711;
c = 5.332;
n1 = 30;
n2 = 20;

% Features positioning
centres = [120, 965];
vals = [1, 2];

% Asteroid model generation
[F, V, N, C20, C22, A, score] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals);

%Astroid knowledge

face_centroids = (V(F(:,1), :) + V(F(:,2), :) + V(F(:,3), :)) / 3;
known_map = double( ...
    face_centroids(:,1) > 0 & ...
    face_centroids(:,2) > 0 & ...
    face_centroids(:,3) > 0 );

%known_map = ones(size(F, 1), 1);

%Features for navigation
nav_index = [   7,   12,   17,   24,   53,   62,   79,   86,   88,  116,  120, ...
  123,  162,  166,  170,  174,  188,  201,  217,  223,  236,  246, ...
  247,  270,  307,  326,  338,  340,  362,  370,  397,  405,  409, ...
  413,  414,  429,  465,  513,  536,  576,  588,  593,  622,  642, ...
  657,  674,  707,  708,  719,  749,  761,  763,  781,  786,  817, ...
  819,  836,  847,  849,  850,  867,  881,  901,  907,  911,  912, ...
  934,  953,  999, 1004, 1009, 1025, 1040, 1041, 1052, 1066, 1090, ...
 1102]; 
 %1105, 1106, 1127, 1205, 1223, 1226, 1227, 1264, 1266, 1282, ...
 %1338, 1343, 1348, 1358, 1385, 1402, 1410, 1415, 1431, 1447, 1456, ...
 %1466];

%Put data in struct
spacecraft_data = struct();
spacecraft_data.data_guidance.ReachabilityScoreComputation = 1;
spacecraft_data.data_guidance.ReachabilityExplorationScheme = 2;
spacecraft_data.data_guidance.DeltaV_max = 2e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 1*3600;
spacecraft_data.data_guidance.DeltaT_after_man = 0.3*3600;
spacecraft_data.data_guidance.r_impact = 24;
spacecraft_data.data_guidance.r_escape = 100;
spacecraft_data.data_guidance.fov1 = 8*pi/180;
spacecraft_data.data_guidance.fov2 = 8*pi/180;

spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.mass = 6.687e15;
spacecraft_data.data_asteroids.omega = 2*pi/(5.27025547*3600)*[0; 0; 1];
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;
spacecraft_data.data_asteroids.mapping.known_map = known_map;
spacecraft_data.data_asteroids.mapping.incidence = [0, 85]*pi/180;
spacecraft_data.data_asteroids.mapping.emission = [0, 85]*pi/180;
spacecraft_data.data_asteroids.features.score = score;
spacecraft_data.data_asteroids.features.known_map_features = known_map;
spacecraft_data.data_asteroids.navigation = nav_index;

model_dyn = @(t, x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22);


%%
profile clear
profile on
[uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU([r0; v0], P0 ,t0,spacecraft_data);
profile off
profile viewer
[xx_opt,tt_opt] = integrateODE([r0; v0+uu_opt], t0:100:th_opt);

[J_of_t_opt, dJdt_opt, new_scores_opt, new_known_map_opt, mapping_score_t_opt, exploit_score_t_opt, nav_score_opt]  = total_score(xx_opt, tt_opt, P0, spacecraft_data);

figure(2)
plot3(xx_opt(:, 1), xx_opt(:, 2), xx_opt(:, 3), 'b', 'LineWidth', 1.5);
hold on
plot3(xx_opt(1, 1), xx_opt(1, 2), xx_opt(1, 3), 'g.', 'MarkerSize', 12);
plot3(xx_opt(end, 1), xx_opt(end, 2), xx_opt(end, 3), 'r.', 'MarkerSize', 12);
plotEllipsoidWithKnownRegion(F, V, new_scores_opt, new_known_map_opt)
plot_sun_pointing_vector(tt_opt);

[xx0,tt] = integrateODE([r0; v0], t0:100:th_opt);
[J_of_t, dJdt, new_scores, new_known_map, mapping_score_t, exploit_score_t, nav_score]  = total_score(xx0, tt, P0, spacecraft_data);

figure(1)
plot3(xx0(:, 1), xx0(:, 2), xx0(:, 3), 'b', 'LineWidth', 1.5);
hold on
plot3(xx0(1, 1), xx0(1, 2), xx0(1, 3), 'g.', 'MarkerSize', 12);
plot3(xx0(end, 1), xx0(end, 2), xx0(end, 3), 'r.', 'MarkerSize', 12);
plotEllipsoidWithKnownRegion(F, V, score, new_known_map)
plot_sun_pointing_vector(tt);

%% MCTS

profile clear
profile on

spacecraft_data_new = spacecraft_data;
r_start = r0;
v_start = v0;
t_start = t0;
P_start = P0;

all_trees = {};
initial_tree = {};
truth_dyn = @(t, x) dynamics15(t, x, mass_eros, omega_body);
for i = 1:1
    max_iter = 50;
    tree = MCTS_belief_based_reduced([r_start; v_start], P_start, t_start,  max_iter, spacecraft_data_new, initial_tree);
    [best_path, best_actions, best_final_times] = find_best_path(tree);

    [TT_real, XX_real] = compute_trajectory(tree{1}.state, best_final_times, best_actions, truth_dyn, options);
    [J_of_t, dJdt, new_scores_real, new_known_map_real, mapping_score_t, exploit_score_t, nav_score] = total_score(XX_real, TT_real, P0, spacecraft_data);
    
    [TT_model, XX_model] = compute_trajectory(tree{1}.state, best_final_times, best_actions, model_dyn, options);
    [J_of_t, dJdt, new_scores_model, new_known_map_model, mapping_score_t, exploit_score_t, nav_score] = total_score(XX_model, TT_model, P0, spacecraft_data);

    all_trees{i} = tree;
end

profile off
profile viewer

slider(F, V, score, new_known_map, XX_model, TT_model)
%%
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_model);
hold on
plot3(XX_real(:, 1), XX_real(:, 2), XX_real(:, 3), 'b', 'LineWidth', 1.5)
plot3(XX_model(:, 1), XX_model(:, 2), XX_model(:, 3), 'r', 'LineWidth', 1.5)
for i = 1:length(best_final_times)
    man_point = find(TT_model == best_final_times(i));
    plot3(XX_model(man_point, 1), XX_model(man_point, 2), XX_model(man_point, 3), 'r.', 'MarkerSize', 15);
end