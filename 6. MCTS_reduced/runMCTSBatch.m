function [all_trees, real_trajectory, filter_trajectory, P_all, tt_all, ...
          total_mapping_score, total_exploiting_score, total_nav_score, ...
          action_times, all_flag, spacecraft_data_out] = ...
    runMCTSBatch(spacecraft_data, r0, v0, t0, P0, iterations, n_trees, options)

%Questa funzione esegue n_trees alberi MCTS da iterations iterazioni in sequenza e aggrega i risultati.

%Set up initial data
spacecraft_data_new = spacecraft_data;
truth_dyn = spacecraft_data.data_guidance.trueDynamics;

real_r_start = r0(:);
real_v_start = v0(:);
plan_r_start = r0(:);
plan_v_start = v0(:);
t_start = t0;
P_start = P0;

all_trees = {};
initial_tree = {};              
all_flag = 1;                  

real_trajectory   = [r0', v0'];
filter_trajectory = [r0', v0'];
P_all  = P0;
tt_all = t0;

total_mapping_score    = [];
total_exploiting_score = [];
total_nav_score        = [];
action_times           = [];

sigma_magn = spacecraft_data.data_guidance.sigma_magn;
sigma_align= spacecraft_data.data_guidance.sigma_align;

for i = 1:n_trees
    % 1) Tree search e best path
    tree = MctsBeliefBased([plan_r_start; plan_v_start], P_start, t_start, iterations, spacecraft_data_new, initial_tree);
    [best_path, planned_best_actions, best_final_times] = find_best_path(tree); 

    % 2) Introduce control error
    Nact = size(planned_best_actions, 2);
    real_best_actions = zeros(size(planned_best_actions));
    P_man_cells = zeros(3,3,Nact);
    for j = 1:Nact
        [uu_real, P_man_j] = pertThrust(planned_best_actions(:, j), sigma_magn, sigma_align);
        real_best_actions(:, j) = uu_real;
        P_man_cells(:,:,j) = P_man_j;
    end

    % 3) Compute real trajectory
    [TT_real_cells, XX_real_cells] = compute_trajectory([real_r_start; real_v_start], best_final_times, real_best_actions, truth_dyn, options);

    % 4) Compute estimated trajectroy
    XX_real = [];
    XX_filter = [];
    P_filtered = [];
    TT_real = [];
    cov_div = tree{1}.cov;           % cov di partenza
    flags = [];
    y0 = [plan_r_start', plan_v_start'];

    for j = 1:length(TT_real_cells)
        % Initial state
        y0 = y0 + [zeros(1,3), planned_best_actions(:, j)'];
        % Add manoeuvre covariance
        P_man_j = P_man_cells(:,:,j);
        cov_div(4:6, 4:6) = cov_div(4:6, 4:6) + P_man_j;

        % EKF
        [P_filtered_parz, XX_filter_parz, eta_f, flag] = navigationFilter(y0, XX_real_cells{j}, cov_div, TT_real_cells{j}, spacecraft_data_new);

        add_x = XX_real_cells{j};
        add_t = TT_real_cells{j};

        if j == 1
            XX_filter = [XX_filter; XX_filter_parz];
            XX_real   = [XX_real;   add_x];
            P_filtered= cat(3, P_filtered, P_filtered_parz);
            TT_real   = [TT_real;   add_t];
            flags     = [flags,     flag];
        else
            XX_filter = [XX_filter; XX_filter_parz(2:end, :)];
            XX_real   = [XX_real;   add_x(2:end, :)];
            P_filtered= cat(3, P_filtered, P_filtered_parz(:, :, 2:end));
            TT_real   = [TT_real;   add_t(2:end)];
            flags     = [flags,     flag(2:end)];
        end

        cov_div = P_filtered(:, :, end);
        y0 = XX_filter(end, :);

        spacecraft_data_new.data_guidance.eta0 = eta_f;

    end

    % 5) Compute obtained score
    parz_mapping_score = [];
    parz_exploiting_score = [];
    parz_nav_score = [];
    
    for j = 1:length(TT_real_cells)
        timeVec = TT_real_cells{j};
        trajVec = XX_real_cells{j};
        ind = find(TT_real == timeVec(1));
        
        P_start_scoring = P_filtered(:, :, ind);
        [J_of_t, dJdt, new_scores_real, new_known_map_real, mapping_score_t, exploit_score_t, nav_score_t] = ...
            total_score(trajVec, timeVec, P_start_scoring, spacecraft_data_new);

        % Aggiorna conoscenza a mappa
        spacecraft_data_new.data_asteroids.features.known_map_features = new_known_map_real;
        spacecraft_data_new.data_asteroids.mapping.known_map           = new_known_map_real;
        spacecraft_data_new.data_asteroids.features.score              = new_scores_real;

        if j == 1
            parz_mapping_score    = [parz_mapping_score;    mapping_score_t];
            parz_exploiting_score = [parz_exploiting_score; exploit_score_t];
            parz_nav_score        = [parz_nav_score;        nav_score_t];
        else
            parz_mapping_score    = [parz_mapping_score;    mapping_score_t(2:end)];
            parz_exploiting_score = [parz_exploiting_score; exploit_score_t(2:end)];
            parz_nav_score        = [parz_nav_score;        nav_score_t(2:end)];
        end
    end

    % 6) Creates vector with the result in time
    if i == 1
        total_mapping_score    = [total_mapping_score;    parz_mapping_score];
        total_exploiting_score = [total_exploiting_score; parz_exploiting_score];
        total_nav_score        = [total_nav_score;        parz_nav_score];
    else
        total_mapping_score    = [total_mapping_score;    parz_mapping_score(2:end)];
        total_exploiting_score = [total_exploiting_score; parz_exploiting_score(2:end)];
        total_nav_score        = [total_nav_score;        parz_nav_score(2:end)];
    end

    P_start  = reshape(P_filtered(:, :, end), 9, 9);
    t_start  = TT_real(end);

    real_trajectory   = [real_trajectory;   XX_real(2:end, :)];
    filter_trajectory = [filter_trajectory; XX_filter(2:end, :)];
    P_all  = cat(3, P_all, P_filtered(:, :, 2:end));
    tt_all = [tt_all; TT_real(2:end)];
    all_flag = [all_flag, flags(2:end)];

    real_r_start = XX_real(end, 1:3)';
    real_v_start = XX_real(end, 4:6)';

    plan_r_start = XX_filter(end, 1:3)';
    plan_v_start = XX_filter(end, 4:6)';

    action_times = [action_times, best_final_times(:)'];

    all_trees{i} = tree;
end

spacecraft_data_out = spacecraft_data_new;
end
