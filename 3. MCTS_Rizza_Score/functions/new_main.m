clc
clear
close all

addpath(genpath('mice\'))
addpath(genpath('sgp4\'))
addpath(genpath('kernels\'))
addpath(genpath('functions\'))
addpath(genpath('3dModel\'))
addpath(genpath('dynamics\'))
addpath(genpath('Reachability'))
addpath(genpath('MCTS_functions'))
addpath(genpath('newScore'))

cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\erosephem_1999004_2002181.bsp')
cspice_furnsh('kernels\erosatt_1998329_2001157_v01.bpc');

%% Perform an iteration of reachability

%Epoch
t0 = cspice_str2et('2000-08-01 T01:00:00');
tf = cspice_str2et('2000-08-01 T05:00:00');

r0 = [8.5855;   44.6428;   -4.4817];
v0 = [0.0164;   -0.0032;    0.0018];

%Data for body frame propagation
T_rotation_eros = 5.27025547*3600; 
omega_body = 2*pi/T_rotation_eros*[0; 0; 1];
mass_eros = 6.687e15;

% FOV of the camera
fov1 = 20*pi/180;
fov2 = 20*pi/180;

% Semiaxis (in km) of the asteroid model 
a = 20.591;
b = 5.711;
c = 5.332;
n1 = 50;
n2 = 30;

centres = [100, 700, 1100, 300];
vals = [20, 50, 30, 40];

% Asteroid model generation
[F, V, N, C20, C22, A, score_map] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals);

%Assumed asteroid fully known
known_map = zeros(size(F, 1), 1);

zones.faces = F;
zones.vertexes = V;
zones.normals = N;
zones.known_map = known_map;
zones.score_map = score_map;

tt = t0:100:tf;
[xx,tt] = integrateODE([r0; v0], tt);

spacecraft_data = struct();
spacecraft_data.data_guidance.ReachabilityScoreComputation = 1;
spacecraft_data.data_guidance.ReachabilityExplorationScheme = 2;
spacecraft_data.data_guidance.DeltaV_max = 2e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 1*3600;
spacecraft_data.data_guidance.DeltaT_after_man = 0.3*3600;
spacecraft_data.data_guidance.r_impact = 24;
spacecraft_data.data_guidance.r_escape = 100;
spacecraft_data.data_guidance.fov1 = 5*pi/180;
spacecraft_data.data_guidance.fov2 = 5*pi/180;

spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.mapping.known_map = known_map;
spacecraft_data.data_asteroids.mapping.incidence = [0, 70]*pi/180;
spacecraft_data.data_asteroids.mapping.emission = [0, 70]*pi/180;
spacecraft_data.data_asteroids.mass = 6.687e15;
spacecraft_data.data_asteroids.omega = 2*pi/(5.27025547*3600)*[0; 0; 1];
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;
spacecraft_data.data_asteroids.features.score_map = score_map;
spacecraft_data.data_asteroids.features.n_observations = zeros(size(F, 1), 1);
spacecraft_data.data_asteroids.features.max_observations = 3*ones(size(F, 1), 1);
spacecraft_data.data_asteroids.features.incidence = [0, 70]*pi/180;
spacecraft_data.data_asteroids.features.emission = [0, 70]*pi/180;

[J_of_t, dJdt, n_observations_1, known_map_1] = score_new(xx, tt, spacecraft_data);

figure(1)
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'b', 'LineWidth', 1.5)
hold on
plotEllipsoidWithKnownRegion(F, V, known_map_1, score_map)
plot_sun_pointing_vector(tt)

spacecraft_data.data_asteroids.mapping.known_map = known_map_1;
spacecraft_data.data_asteroids.features.n_observations = n_observations_1;
[J_of_t, dJdt, n_observations_2, known_map_2] = score_new(xx, tt, spacecraft_data);

spacecraft_data.data_asteroids.mapping.known_map = known_map_2;
spacecraft_data.data_asteroids.features.n_observations = n_observations_2;
[J_of_t, dJdt, n_observations_3, known_map_3] = score_new(xx, tt, spacecraft_data);

spacecraft_data.data_asteroids.mapping.known_map = known_map_3;
spacecraft_data.data_asteroids.features.n_observations = n_observations_3;
[J_of_t, dJdt, n_observations_4, known_map_4] = score_new(xx, tt, spacecraft_data);

spacecraft_data.data_asteroids.mapping.known_map = known_map_3;
spacecraft_data.data_asteroids.features.n_observations = n_observations_4;
[J_of_t, dJdt, n_observations_5, known_map_5] = score_new(xx, tt, spacecraft_data);
%%
%Create spacecraft data
spacecraft_data = struct();
spacecraft_data.data_guidance.ReachabilityScoreComputation = 1;
spacecraft_data.data_guidance.ReachabilityExplorationScheme = 2;
spacecraft_data.data_guidance.DeltaV_max = 2e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 1*3600;
spacecraft_data.data_guidance.DeltaT_after_man = 0.3*3600;
spacecraft_data.data_guidance.r_impact = 24;
spacecraft_data.data_guidance.r_escape = 100;
spacecraft_data.data_guidance.fov1 = 5*pi/180;
spacecraft_data.data_guidance.fov2 = 5*pi/180;

spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.features = features;
spacecraft_data.data_asteroids.known_map = known_map;
spacecraft_data.data_asteroids.mapping.incidence = [0, 90]*pi/180;
spacecraft_data.data_asteroids.mapping.emission = [0, 80]*pi/180;
spacecraft_data.data_asteroids.mapping.range = [25, 50];
spacecraft_data.data_asteroids.mass = 6.687e15;
spacecraft_data.data_asteroids.omega = 2*pi/(5.27025547*3600)*[0; 0; 1];
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;

%Data for the propagation of n optimal arcs from the initial conditions
n = 3;
x_start = [r0; v0];
uu_opt_vec = zeros(3, n);
th_opt_vec = zeros(n+1, 1);

th_opt_vec(1) = t0;

h = waitbar(0, 'Processing...');

%Run reachability n times
for i = 1:n

    waitbar(i / n, h, sprintf('Progress: %d%%', round(100 * i / n)));

    [uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU(x_start , th_opt_vec(i), spacecraft_data);

    %Optimal Initial condition
    uu_opt_vec(:, i) = uu_opt;
    th_opt_vec(i+1) = th_opt;
    
    %Propagate dynamics with the optimal control
    x_opt = x_start + [zeros(3, 1); uu_opt];
    [tt_prop, xx_prop] = ode78(@(t,x) model_rotating_dynamics(t, x, mass_eros, omega_body, C20, C22), linspace(th_opt_vec(i), th_opt_vec(i+1), 200), x_opt, options);
    [~, ~, features_new, known_map_new] = score(xx_prop, tt_prop', spacecraft_data);
    
    %Update feature score
    spacecraft_data.data_asteroids.features = features_new;
    spacecraft_data.data_asteroids.known_map = known_map_new;

    x_start = xx_prop(end, :)';
end

close(h);


model_dyn = @(t, x) model_rotating_dynamics(t, x, mass_eros, omega_body, C20, C22);
[TT_classic, XX_classic] = compute_trajectory([r0; v0], th_opt_vec, uu_opt_vec, model_dyn, options);
