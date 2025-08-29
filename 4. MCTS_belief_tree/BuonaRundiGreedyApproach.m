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

%% GREEDY APPROACH ONLY MAPPING
% Initial conditions

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
vals = [0, 0];

% Asteroid model generation
[F, V, N, C20, C22, A, score] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals);

%Astroid knowledge

face_centroids = (V(F(:,1), :) + V(F(:,2), :) + V(F(:,3), :)) / 3;
known_map = double( ...
    face_centroids(:,1) > 0 & ...
    face_centroids(:,2) > 0 & ...
    face_centroids(:,3) > 0 );

known_map = zeros(size(F, 1), 1);

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
spacecraft_data.data_guidance.DeltaV_max = 3e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 2*3600;
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
truth_dyn = @(t, x) dynamics15(t, x, mass_eros, omega_body);


%% 
profile clear
profile on
n = 15;
spacecraft_data_new = spacecraft_data;

r_start = r0;
v_start = v0;
t_start = t0;
P_start = P0;

scores = {};
known_maps = {};
trajectory = [];
covs = P0;

map_scores = [];
exploit_scores = [];
nav_scores = [];
tt_filt = [];
tt_total = [];
t_scores = [];
man_point = [];
arc_length = [];

for i = 1:n
    [uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU([r_start; v_start], P_start , t_start ,spacecraft_data_new);
    
    [tt_opt, xx_opt] = compute_trajectory([r_start; v_start], [t_start, th_opt], uu_opt, truth_dyn, spacecraft_data_new);

    [J_of_t_opt, dJdt_opt, new_scores_opt, new_known_map_opt, mapping_score_t_opt, exploit_score_t_opt, nav_score_opt]  = total_score(xx_opt, tt_opt, P_start, spacecraft_data_new);

    spacecraft_data_new.data_asteroids.features.score = new_scores_opt;
    spacecraft_data_new.data_asteroids.features.known_map_features = new_known_map_opt;
    spacecraft_data_new.data_asteroids.mapping.known_map = new_known_map_opt;

    r_start = xx_opt(end, 1:3)';
    v_start = xx_opt(end, 4:6)';
    t_start = tt_opt(end);

    trajectory = [trajectory; xx_opt];
    man_point = [man_point; xx_opt(end, 1:3)];

    map_scores = [map_scores; mapping_score_t_opt];

    tt_total = [tt_total; tt_opt];

    arc_length = [arc_length; (tt_opt(end)-tt_opt(1))/3600];
    
end
profile off
profile viewer


figure(1)
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_opt);
hold on
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'b');
for i = 1:size(man_point, 1)
    plot3(man_point(i, 1), man_point(i, 2), man_point(i, 3), 'b.', 'MarkerSize', 12)
end

figure(2)
plot((tt_total-tt_total(1))/(3600), cumsum(map_scores), 'b', 'LineWidth', 1.5)
hold on
yline(1)

trajectory_in = zeros(size(trajectory(:, 1:3)));
for i = 1:length(tt_total)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_total(i)); 
    trajectory_in(i, :) = (R_body2in*trajectory(i, 1:3)')';
end
figure(3)
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_opt);
plot3(trajectory_in(:, 1), trajectory_in(:, 2), trajectory_in(:, 3), 'b')

figure(4)
plot(arc_length, 'b', 'Marker', '*')

figure(5)
plot(vecnorm(trajectory(:, 1:3), 2, 2))


%% GREEDY APPROACH ONLY EXPLOIING 
% Initial conditions

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
centres = [120, 965, 640, 810];
vals = [1, 2, 1, 2, 1];

% Asteroid model generation
[F, V, N, C20, C22, A, score] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals);

%Astroid knowledge

face_centroids = (V(F(:,1), :) + V(F(:,2), :) + V(F(:,3), :)) / 3;
known_map = double( ...
    face_centroids(:,1) > 0 & ...
    face_centroids(:,2) > 0 & ...
    face_centroids(:,3) > 0 );

known_map = ones(size(F, 1), 1);

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
spacecraft_data.data_guidance.safety_margin = 2*3600;
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
truth_dyn = @(t, x) dynamics15(t, x, mass_eros, omega_body);


%% 
profile clear
profile on
n = 5;
spacecraft_data_new = spacecraft_data;

r_start = r0;
v_start = v0;
t_start = t0;
P_start = P0;

scores = {};
known_maps = {};
trajectory = [];
covs = P0;

map_scores = [];
exploit_scores = [];
nav_scores = [];
tt_filt = [];
tt_total = [];
t_scores = [];
man_point = [];
arc_length = [];

for i = 1:n
    [uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU([r_start; v_start], P_start , t_start ,spacecraft_data_new);
    
    [tt_opt, xx_opt] = compute_trajectory([r_start; v_start], [t_start, th_opt], uu_opt, truth_dyn, spacecraft_data_new);

    [J_of_t_opt, dJdt_opt, new_scores_opt, new_known_map_opt, mapping_score_t_opt, exploit_score_t_opt, nav_score_opt]  = total_score(xx_opt, tt_opt, P_start, spacecraft_data_new);

    spacecraft_data_new.data_asteroids.features.score = new_scores_opt;
    spacecraft_data_new.data_asteroids.features.known_map_features = new_known_map_opt;
    spacecraft_data_new.data_asteroids.mapping.known_map = new_known_map_opt;

    r_start = xx_opt(end, 1:3)';
    v_start = xx_opt(end, 4:6)';
    t_start = tt_opt(end);

    trajectory = [trajectory; xx_opt];
    man_point = [man_point; xx_opt(end, 1:3)];

    map_scores = [map_scores; mapping_score_t_opt];
    exploit_scores = [exploit_scores; exploit_score_t_opt] 

    tt_total = [tt_total; tt_opt];

    arc_length = [arc_length; (tt_opt(end)-tt_opt(1))/3600];
    
end
profile off
profile viewer


figure(1)
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_opt);
hold on
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'b');
for i = 1:size(man_point, 1)
    plot3(man_point(i, 1), man_point(i, 2), man_point(i, 3), 'b.', 'MarkerSize', 12)
end

figure(2)
plot((tt_total-tt_total(1))/(3600), cumsum(exploit_scores), 'b', 'LineWidth', 1.5)
hold on
yline(1)

trajectory_in = zeros(size(trajectory(:, 1:3)));
for i = 1:length(tt_total)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_total(i)); 
    trajectory_in(i, :) = (R_body2in*trajectory(i, 1:3)')';
end
figure(3)
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_opt);
plot3(trajectory_in(:, 1), trajectory_in(:, 2), trajectory_in(:, 3), 'b')

figure(4)
plot(arc_length, 'b', 'Marker', '*')

figure(5)
plot(vecnorm(trajectory(:, 1:3), 2, 2))

%% GREEDY APPROACH MAPPING + EXPLOIING 
% Initial conditions

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
centres = [120, 965, 637, 810];
vals = [1, 2, 1, 2, 1];

% Asteroid model generation
[F, V, N, C20, C22, A, score] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals);

%Astroid knowledge

face_centroids = (V(F(:,1), :) + V(F(:,2), :) + V(F(:,3), :)) / 3;
known_map = double( ...
    face_centroids(:,1) > 0 & ...
    face_centroids(:,2) > 0 & ...
    face_centroids(:,3) > 0 );

known_map = zeros(size(F, 1), 1);

%Features for navigation
nav_index = [
  7,  12,  17,  24,  35,  53,  58,  62,  70,  79, ...
   86,  88,  97,  116,  120,  123,  140,  150,  162,  166, ...
   170,  174,  180,  188,  201,  217,  223,  230,  236,  246, ...
   247,  260,  270,  290,  307,  315,  326,  338,  340,  362, ...
   370,  383,  397,  405,  409,  413,  414,  420,  429,  465, ...
   490,  513,  536,  560,  576,  588,  593,  610,  622,  642, ...
   657,  674,  690,  707,  708,  719,  730,  749,  761,  763, ...
   781,  786,  800,  817,  819,  836,  847,  849,  850,  867, ...
   881,  890,  901,  907,  911,  912,  934,  953,  999,  1001, ...
   1004,  1009,  1025,  1030,  1040,  1041,  1052,  1066,  1080,  1090, ...
   1102
]; 
 

%Put data in struct
spacecraft_data = struct();
spacecraft_data.data_guidance.ReachabilityScoreComputation = 1;
spacecraft_data.data_guidance.ReachabilityExplorationScheme = 2;
spacecraft_data.data_guidance.DeltaV_max = 2e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 2*3600;
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
truth_dyn = @(t, x) dynamics15(t, x, mass_eros, omega_body);


%% 
profile clear
profile on
n = 1;
spacecraft_data_new = spacecraft_data;

r_start = r0;
v_start = v0;
t_start = t0;
P_start = P0;

scores = {};
known_maps = {};
trajectory = [];
covs = P0;

map_scores = [];
exploit_scores = [];
nav_scores = [];
tt_filt = [];
tt_total = [];
t_scores = [];
man_point = [];
arc_length = [];

for i = 1:n
    [uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU([r_start; v_start], P_start , t_start ,spacecraft_data_new);
    
    [tt_opt, xx_opt] = compute_trajectory([r_start; v_start], [t_start, th_opt], uu_opt, truth_dyn, spacecraft_data_new);

    [J_of_t_opt, dJdt_opt, new_scores_opt, new_known_map_opt, mapping_score_t_opt, exploit_score_t_opt, nav_score_opt]  = total_score(xx_opt, tt_opt, P_start, spacecraft_data_new);

    spacecraft_data_new.data_asteroids.features.score = new_scores_opt;
    spacecraft_data_new.data_asteroids.features.known_map_features = new_known_map_opt;
    spacecraft_data_new.data_asteroids.mapping.known_map = new_known_map_opt;

    r_start = xx_opt(end, 1:3)';
    v_start = xx_opt(end, 4:6)';
    t_start = tt_opt(end);

    trajectory = [trajectory; xx_opt];
    man_point = [man_point; xx_opt(end, 1:3)];

    map_scores = [map_scores; mapping_score_t_opt];
    exploit_scores = [exploit_scores; exploit_score_t_opt];

    tt_total = [tt_total; tt_opt];

    arc_length = [arc_length; (tt_opt(end)-tt_opt(1))/3600];
    
end
profile off
profile viewer
%%

figure(1)
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_opt);
hold on
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'b');
for i = 1:size(man_point, 1)
    plot3(man_point(i, 1), man_point(i, 2), man_point(i, 3), 'b.', 'MarkerSize', 12)
end
title('Trajectory in body fixed rotating frame')

figure(2)
expl = plot((tt_total-tt_total(1))/(3600), cumsum(exploit_scores), 'b', 'LineWidth', 1.5)
hold on
yline(1)
mapp = plot((tt_total-tt_total(1))/(3600), cumsum(map_scores), 'g', 'LineWidth', 1.5)
title('Cumulative score obtained')
xlabel('time [h]')
ylabel('score [-]')
legend([expl, mapp], 'Exploiting Score', 'Mapping Score')

trajectory_in = zeros(size(trajectory(:, 1:3)));
for i = 1:length(tt_total)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_total(i)); 
    trajectory_in(i, :) = (R_body2in*trajectory(i, 1:3)')';
end
figure(3)
plotEllipsoidWithKnownRegion(F, V, score, new_known_map_opt);
plot3(trajectory_in(:, 1), trajectory_in(:, 2), trajectory_in(:, 3), 'b')
plot_sun_pointing_vector(tt_total, 'ECLIPJ2000')
title('Trajectory in ECLIPJ2000')

figure(4)
plot(arc_length, 'b', 'Marker', '*')
title('Temporal lenght of each arc')
xlabel('Arc number')
ylabel('Temporal length [h]')

figure(5)
plot(tt_total, vecnorm(trajectory(:, 1:3), 2, 2))
title('Range over time')
xlabel('Time [h]')
ylabel('Range [km]')


velocity_in = zeros(size(trajectory(:, 4:6)));
for i = 1:length(tt_total)
    R_body2in = cspice_pxform('IAU_EROS', 'ECLIPJ2000', tt_total(i)); 
    velocity_in(i, :) = (R_body2in*trajectory(i, 4:6)')' + ( R_body2in*cross(omega_body, trajectory(i, 1:3)') )';
end
figure(6)
plot(tt_total, vecnorm(velocity_in, 2, 2)*1000)
title('Velocity magnitude over time in ECLIPJ2000')
xlabel('Time [h]')
ylabel('Velocity [m/s]')

