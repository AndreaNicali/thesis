function [real_trajectory, filter_trajectory, P_all, tt_all, ...
          total_mapping_score, total_exploiting_score, total_nav_score, ...
          action_times, all_flag, action_list, spacecraft_data_out] = ...
    greedyApproach(spacecraft_data, r0, v0, t0, P0, n_actions, n_set, options)


model_dyn = spacecraft_data.data_guidance.modelDynamics;
truth_dyn = spacecraft_data.data_guidance.trueDynamics;

%Set up initial data
spacecraft_data_plan = spacecraft_data;
spacecraft_data_real = spacecraft_data;

real_r_start = r0(:);
real_v_start = v0(:);
plan_r_start = r0(:);
plan_v_start = v0(:);
filter_r_start = r0(:);
filter_v_start = v0(:);
t_start = t0;
plan_P_start = P0;
real_P_start = P0;

all_flag = 1;                  

real_trajectory   = [r0', v0'];
filter_trajectory = [r0', v0'];
P_all  = P0;
tt_all = t0;

total_mapping_score    = 0;
total_exploiting_score = 0;
total_nav_score        = 0;
action_times           = t_start;
action_times_log       = action_times;

planned_actions = [];
action_list = {};

sigma_magn = spacecraft_data.data_guidance.sigma_magn;
sigma_align= spacecraft_data.data_guidance.sigma_align;

if isscalar(n_actions)
    n_actions = n_actions*ones(1, n_set);
end

for i = 1:n_set

    % 1) Tree search e best path
    for j = 1:n_actions(i)
        [~,~,~,U,J_val,T,~,I] = exploreU([plan_r_start; plan_v_start], plan_P_start, t_start, spacecraft_data_plan);
        [~, idBest] = max(J_val);
        action = U(:, idBest);
        tf = T(idBest);
        
        [~, P_man] = pertThrust(action, sigma_magn, sigma_align);
        plan_P_start(4:6, 4:6) = plan_P_start(4:6, 4:6) + P_man;
    
        %Propagate trajectory
        [tt, xx] = ode78(@(t,x) model_dyn(t, x), t_start:100:tf, [plan_r_start; plan_v_start+action], options);
        
        %Propagate Uncertainties
        [P_filtered, xx_filt, eta_f] = navigationFilter([plan_r_start', plan_v_start'+action'], xx, plan_P_start, tt, spacecraft_data_plan);
        Pf = P_filtered(:, :, end);
    
        [~ , ~, new_scores, new_known_map] = total_score(xx_filt, tt, plan_P_start, spacecraft_data_plan);
        
        %Put updated features and map data
        spacecraft_data_plan.data_asteroids.features.score = new_scores;
        spacecraft_data_plan.data_asteroids.mapping.known_map = new_known_map;
        spacecraft_data_plan.data_guidance.eta0 = eta_f;
        
        planned_actions = [planned_actions, action];
        action_times = [action_times, tf];
        action_times_log = [action_times_log, tf];

        plan_r_start = xx_filt(end, 1:3).';
        plan_v_start = xx_filt(end, 4:6).';
        plan_P_start = Pf;
        t_start = tf;
        
    end
    
    
    for j = 1:size(planned_actions, 2)
        [uu_real, P_man] = pertThrust(planned_actions(:, j), sigma_magn, sigma_align);
        t_start = action_times(j);
        tf = action_times(j+1);
        real_P_start(4:6, 4:6) = real_P_start(4:6, 4:6) + P_man; 

        %Propagate trajectory
        [tt, xx] = ode78(@(t,x) truth_dyn(t, x), t_start:100:tf, [real_r_start; real_v_start+uu_real], options);
        
        %Propagate Uncertainties
        [P_filtered, xx_filt, eta_f] = navigationFilter([filter_r_start', filter_v_start' + planned_actions(:, j)'], xx, real_P_start, tt, spacecraft_data_real);
        Pf = P_filtered(:, :, end);
    
        [~ , ~, new_scores, new_known_map, mapping_score_t, exploit_score_t, nav_score_t] = total_score(xx, tt, real_P_start, spacecraft_data_real);
        
        %Put updated features and map data
        spacecraft_data_real.data_asteroids.features.score = new_scores;
        spacecraft_data_real.data_asteroids.mapping.known_map = new_known_map;
        spacecraft_data_real.data_asteroids.features.known_map_features = new_known_map;        
        spacecraft_data_real.data_guidance.eta0 = eta_f;
        
        real_r_start = xx(end, 1:3).';
        real_v_start = xx(end, 4:6).';
        real_P_start = Pf;
        t_start = tf;

        filter_r_start = xx_filt(end, 1:3).';
        filter_v_start = xx_filt(end, 4:6).';

        real_trajectory   = [real_trajectory;   xx(:,1:6)];
        filter_trajectory = [filter_trajectory; xx_filt(:,1:6)];
        P_all             = cat(3, P_all, P_filtered);
        tt_all            = [tt_all; tt(:)];
        total_mapping_score    = [total_mapping_score;    mapping_score_t(:)];
        total_exploiting_score = [total_exploiting_score; exploit_score_t(:)];
        total_nav_score        = [total_nav_score;        nav_score_t(:)];
    end

    plan_r_start = filter_r_start;
    plan_v_start = filter_v_start;
    t_start      = action_times(end);   
    action_list{i} = planned_actions;
    spacecraft_data_plan = spacecraft_data_real;

    planned_actions = [];
    action_times    = t_start;
end

action_times = action_times_log;
spacecraft_data_out = spacecraft_data_real;
