%In this script the reachability analyis is performed for n consecutive
%arcs. The score function is updated after each arcs to verify goal
%completeness state after the process. The score function is very similar
%to the one used by A. Rizza in his PhD thesis.
%The asteroid is assumed to be completely known as well as the feature position.
%The total score obtained with the n consecutive call to the reachability 
%function (greedy approach) is compared with the score obtained with a MCTS
%based approach.

clc
clear
close all

addpath(genpath('kernels\'))
addpath(genpath('functions\'))
addpath(genpath('3dModel\'))
addpath(genpath('dynamics\'))
addpath(genpath('Uncertainties\'))
addpath(genpath('Reachability\'))
addpath(genpath('MCTS_functions\'))

cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\erosephem_1999004_2002181.bsp')
cspice_furnsh('kernels\erosatt_1998329_2001157_v01.bpc');
%% Perform an iteration of reachability

%Initial Epoch
t0 = cspice_str2et('2000-08-01 T00:00:00');

%Initial conditions in fixed body frame
r0 = [6.85775061407898; 
      -40.0375655948265;
  	    8.70203252796758];
v0 = [-0.0152783171329873;
     -0.00131932988095069;
 	  0.00165722605576695];
r0 =  [-29.6631;
      -30.0088;
      -15.6549];
v0 = [-0.0119;
      0.0104;
      0.0014];
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
n1 = 40;
n2 = 20;

% Asteroid model generation
[F, V, N, C20, C22, A] = Ellipsoid(a, b, c, n1, n2);

%Options for propagation in the model dynamic
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

%Feature struct, contains feature index, total score achievable, level of
%completeness, actual score, available score, and requirements of
%incidence, emission and relative angle between the two

features.index = [86 116 1041 536 622 867 1052 1227 1402 24 201 413 749 907 1102 1358]';

features.totalscore = [840, 840, 840, 840, 840, 840, 840, 840, 840, 840, 840, 840, 840, 840, 840, 840]';
features.completeness = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
features.score = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
features.availablescore = features.totalscore - features.score;

features.visualization.emission = [0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10; 0, 10]*pi/180;
features.visualization.incidence = [0, 90; 0, 90; 0, 90; 0, 90; 0, 90; 0, 90; 0, 90; 0, 90; 0, 90;  0, 90; 0, 90; 0, 90; 0, 90; 0, 90; 0, 90; 0, 90]*pi/180;
features.visualization.range = [35, 41; 35, 41; 35, 41; 35, 41; 35, 41; 35, 41;35, 41; 35, 41; 35, 41; 35, 41; 35, 41; 35, 41; 35, 41;35, 41; 35, 41; 35, 41];

%Assumed asteroid fully known
known_map = ones(size(F, 1), 1);

%Identify each feature as known or unknown
for i = 1:length(features.index)
    idx = features.index(i);
    % Calcola il centroide della faccia
    face_verts = V(F(idx, :), :);
    centroid = mean(face_verts, 1);
    if known_map(features.index(i))
        features.known(i) = 1;
    else 
        features.known(i) = 0;
    end
end

%Create spacecraft data struct
spacecraft_data = struct();
spacecraft_data.data_guidance.ReachabilityScoreComputation = 1;
spacecraft_data.data_guidance.ReachabilityExplorationScheme = 2;
spacecraft_data.data_guidance.DeltaV_max = 2e-3;
spacecraft_data.data_guidance.Th_max = 7*3600;
spacecraft_data.data_guidance.safety_margin = 1*3600;
spacecraft_data.data_guidance.DeltaT_after_man = 0.3*3600;
spacecraft_data.data_guidance.r_impact = 24;
spacecraft_data.data_guidance.r_escape = 100;
spacecraft_data.data_guidance.fov1 = fov1;
spacecraft_data.data_guidance.fov2 = fov2;

spacecraft_data.data_asteroids.Faces = F;
spacecraft_data.data_asteroids.Vertexes = V;
spacecraft_data.data_asteroids.Normals = N;
spacecraft_data.data_asteroids.features = features;
spacecraft_data.data_asteroids.known_map = known_map;
spacecraft_data.data_asteroids.mapping.incidence = [0, 90];
spacecraft_data.data_asteroids.mapping.emission = [0, 80];
spacecraft_data.data_asteroids.mapping.range = [25, 50];
spacecraft_data.data_asteroids.mass = 6.687e15;
spacecraft_data.data_asteroids.omega = 2*pi/(5.27025547*3600)*[0; 0; 1];
spacecraft_data.data_asteroids.C20 = C20;
spacecraft_data.data_asteroids.C22 = C22;

%% PROPAGATE n OPTIMAL ARCS: GREEDY APPROACH
%Data for the propagation of n optimal arcs from the initial conditions
n = 10;
uu_opt_vec = zeros(3, n);
th_opt_vec = zeros(n+1, 1);
th_opt_vec(1) = t0;

x_start = [r0; v0];
spacecraft_data_upd = spacecraft_data;
h = waitbar(0, 'Running Greedy approach...');

%Run reachability n times
for i = 1:n
    
    [uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU(x_start , th_opt_vec(i), spacecraft_data_upd);

    %Optimal Initial condition
    uu_opt_vec(:, i) = uu_opt;
    th_opt_vec(i+1) = th_opt;
    
    %Propagate dynamics with the optimal control
    x_opt = x_start + [zeros(3, 1); uu_opt];
    [tt_prop, xx_prop] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), th_opt_vec(i):100:th_opt_vec(i+1), x_opt, options);
    [~, ~, features_new, known_map_new] = score_rizza(xx_prop, tt_prop', spacecraft_data_upd);
    
    %Update feature score
    spacecraft_data_upd.data_asteroids.features = features_new;
    spacecraft_data_upd.data_asteroids.known_map = known_map_new;

    x_start = xx_prop(end, :)';

    waitbar(i / n, h, sprintf('Progress: %d%%', round(100 * i / n)));

end

close(h);

model_dyn = @(t, x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22);
[TT_GREEDY, XX_GREEDY] = compute_trajectory([r0; v0], th_opt_vec, uu_opt_vec, model_dyn, options);

%% MCTS
profile clear
profile on
x0 = [r0;v0];
iter = 50;

tree = MCTS(x0, t0, iter, spacecraft_data);
profile off
profile viewer

[best_path, best_actions, best_times] = find_best_path(tree);

[TT_MCTS, XX_MCTS] = compute_trajectory(x0, best_times, best_actions, model_dyn, options);

%% PLOT
figure(1)
plotEllipsoidWithKnownRegion(F, V, known_map, known_map)
hold on

%Plot feature position
for i = 1:length(features.index)
    idx = features.index(i);
    face_verts = V(F(idx, :), :);
    centroid = mean(face_verts, 1);
    if known_map(features.index(i))
        plot3(centroid(1), centroid(2), centroid(3), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    else 
        plot3(centroid(1), centroid(2), centroid(3), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    end
end

plot_visibility_volumes(features, F, V, N)

plot3(XX_GREEDY(:, 1), XX_GREEDY(:, 2), XX_GREEDY(:, 3), 'b', 'LineWidth', 1.5);

plot3(XX_MCTS(:, 1), XX_MCTS(:, 2), XX_MCTS(:, 3), 'r', 'LineWidth', 1.5);

[J_classic, ~, features_classic, ~] = score_rizza(XX_GREEDY, TT_GREEDY', spacecraft_data);
[J_MCTS, ~, features_MCTS, ~] = score_rizza(XX_MCTS, TT_MCTS', spacecraft_data);

figure(2);

% Assumendo che entrambe abbiano la stessa lunghezza
data = [features_classic.completeness(:), features_MCTS.completeness(:)];

bar(data, 'grouped');
xlabel('Feature index');
ylabel('Completeness level');
title('Bar chart of completeness values');
legend('Classic', 'MCTS');

figure(3)

plot(TT_GREEDY, J_classic, 'b', 'LineWidth', 1.5);
hold on
plot(TT_MCTS, J_MCTS, 'r', 'LineWidth', 1.5);
xline(th_opt_vec, 'b')
xline(best_times, 'r')

score_per_impulse_classic = [];
for i = 1:length(th_opt_vec)
    val = J_classic(TT_GREEDY == th_opt_vec(i));
    if ~isscalar(val)
        val = val(end);
    end

    score_per_impulse_classic(i) = val;
end

score_per_impulse_MCTS = [];
for i = 1:length(best_times)
    
    val = J_MCTS(TT_MCTS == best_times(i));
    if ~isscalar(val)
        val = val(end);
    end
    score_per_impulse_MCTS(i) = val;
end

figure(4)
plot( linspace(0, length(score_per_impulse_MCTS),  length(score_per_impulse_MCTS)+1), [0, score_per_impulse_MCTS], 'r' );
hold on
plot( linspace(0, length(score_per_impulse_classic),  length(score_per_impulse_classic)+1), [0, score_per_impulse_classic], 'b' );

