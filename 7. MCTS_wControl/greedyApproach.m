function [real_trajectory, filter_trajectory, ref_trajectory, P_all, real_times, ...
          mapping_scores, exploiting_scores, nav_scores, ...
          action_times, all_flag, planned_actions, real_trajectory_total, ref_trajectory_total, filter_trajectory_total, ...
          tt_all, mapping_score_total, exploiting_score_total, nav_score_total, P_total, spacecraft_data_out] = ...
    greedyApproach(spacecraft_data, r0, v0, t0, P0, n_actions, n_set, options)

%Upload data for propagations
model_dyn = spacecraft_data.data_guidance.modelDynamics;
truth_dyn = spacecraft_data.data_guidance.trueDynamics;
t_plan = spacecraft_data.data_guidance.planningTime;
sigma_magn = spacecraft_data.data_guidance.sigma_magn;
sigma_align= spacecraft_data.data_guidance.sigma_align;

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

real_trajectory   = cell(n_set, max(n_actions)+1);
filter_trajectory = cell(n_set, max(n_actions)+1);
ref_trajectory = cell(n_set, max(n_actions)+1);
ref_times = cell(n_set, max(n_actions)+1);
P_all  = cell(n_set, max(n_actions)+1);
real_times = cell(n_set, max(n_actions)+1);

mapping_scores    = cell(n_set, max(n_actions)+1);
exploiting_scores = cell(n_set, max(n_actions)+1);
nav_scores        = cell(n_set, max(n_actions)+1);
action_times           = cell(n_set, max(n_actions)+1);
final_times       = cell(n_set, max(n_actions)+1);
planned_actions = cell(n_set, max(n_actions)+1);

if isscalar(n_actions)
    n_actions = n_actions*ones(1, n_set);
end

for i = 1:n_set
    
    %PLANNING PHASE greedy approach
    for j = 1:n_actions(i)

        %Compute j-th actions for the i-th set
        [~,~,~,U,J_val,T,~,I] = exploreU([plan_r_start; plan_v_start], plan_P_start, t_start, spacecraft_data_plan);
        [~, idBest] = max(J_val);
        action = U(:, idBest);
        tf = T(idBest);
        
        %Add maneuver perturb covariance
        [~, P_man1] = pertThrust(action, sigma_magn, sigma_align);
        plan_P_start(4:6, 4:6) = plan_P_start(4:6, 4:6) + P_man1;
    
        %Propagate and store reference trajectory
        [tt, xx] = ode78(@(t,x) model_dyn(t, x), t_start:100:tf, [plan_r_start; plan_v_start+action], options);
        ref_trajectory{i, j} = xx;
        ref_times{i, j} = tt;

        %Propagate Uncertainties
        [P_filtered, xx_filt, eta_f] = navigationFilter([plan_r_start', plan_v_start'+action'], xx, plan_P_start, tt, spacecraft_data_plan);
        Pf = P_filtered(:, :, end);
        
        %Compute expected score
        [~ , ~, new_scores, new_known_map] = total_score(xx_filt, tt, plan_P_start, spacecraft_data_plan);
        
        %Put updated features and map data
        spacecraft_data_plan.data_asteroids.features.score = new_scores;
        spacecraft_data_plan.data_asteroids.mapping.known_map = new_known_map;
        spacecraft_data_plan.data_guidance.eta0 = eta_f;
        
        %Register planned actions data
        planned_actions{i, j} = action;
        action_times{i, j} = t_start;
        final_times{i, j} = tf;
        
        %Initialize next iteration value
        plan_r_start = xx_filt(end, 1:3).';
        plan_v_start = xx_filt(end, 4:6).';
        plan_P_start = Pf;
        t_start = tf;
        
    end
    
    dvCorr2 = zeros(3,1);

    %COASTING PHASE greedy approach
    for j = 1:n_actions(i)

        %Compute the correction to the planned action
        dvj = planned_actions{i,j} + dvCorr2; %sum second correction of the previous iteration
        [dvCorr1, dvCorr2] = refTargeting([filter_r_start; filter_v_start+dvj], ref_times{i,j}, ref_trajectory{i,j}, spacecraft_data);        
        actCorr1 = dvj + dvCorr1;

        %Perform (perturbed) planned action
        [uu_real1, P_man1] = pertThrust(actCorr1, sigma_magn, sigma_align);
        t_start = action_times{i,j};
        tf = final_times{i, j};
        real_P_start(4:6, 4:6) = real_P_start(4:6, 4:6) + P_man1; 

        %Propagate real trajectory
        [tt, xx] = ode78(@(t,x) truth_dyn(t, x), t_start:100:tf, [real_r_start; real_v_start+uu_real1], options);
        
        %Propagate Uncertainties
        [P_filtered, xx_filt, eta_f] = navigationFilter([filter_r_start', filter_v_start' + actCorr1'], xx, real_P_start, tt, spacecraft_data_real);
        Pf = P_filtered(:, :, end);
        
        %Compute real score
        [~ , ~, new_scores, new_known_map, mapping_score_t, exploit_score_t, nav_score_t] = total_score(xx, tt, real_P_start, spacecraft_data_real);
        
        %Put updated features and map data
        spacecraft_data_real.data_asteroids.features.score = new_scores;
        spacecraft_data_real.data_asteroids.mapping.known_map = new_known_map;
        spacecraft_data_real.data_asteroids.features.known_map_features = new_known_map;        
        spacecraft_data_real.data_guidance.eta0 = eta_f;
        
        %Store data
        real_trajectory{i, j}   = xx;
        filter_trajectory{i, j} = xx_filt(:,1:6);
        P_all{i, j}             = P_filtered;
        real_times{i, j}            = tt(:);

        mapping_scores{i, j}    = mapping_score_t(:);
        exploiting_scores{i, j} = exploit_score_t(:);
        nav_scores{i, j}        = nav_score_t(:);

        real_r_start = xx(end, 1:3).';
        real_v_start = xx(end, 4:6).';
        real_P_start = Pf;

        filter_r_start = xx_filt(end, 1:3).';
        filter_v_start = xx_filt(end, 4:6).';
        
    end
                                                                                                                                        
    %Correction of the last arc of the batch
    [uu_real, P_man] = pertThrust(dvCorr2, sigma_magn, sigma_align);

    real_P_start(4:6, 4:6) = real_P_start(4:6, 4:6) + P_man;

    t0Planning = tf;
    tfPlanning = tf+t_plan;

    %Account for PLANNING PHASE time
    [tt, xx_real_plan] = ode78(@(t,x) truth_dyn(t, x), t0Planning:100:tfPlanning, [real_r_start; real_v_start+uu_real], options);
    [tt, xx_ref_plan] = ode78(@(t,x) model_dyn(t, x), t0Planning:100:tfPlanning, [plan_r_start; plan_v_start+dvCorr2], options);
    [P_filtered, xx_filt_plan, eta_f] = navigationFilter([filter_r_start', filter_v_start' + dvCorr2'], xx_real_plan, real_P_start, tt, spacecraft_data_real);
    
    real_trajectory{i,n_actions(i)+1} = xx_real_plan;
    ref_trajectory{i,n_actions(i)+1} = xx_ref_plan;
    filter_trajectory{i, n_actions(i)+1} = xx_filt_plan(:, 1:6);
    P_all{i, n_actions(i)+1} = P_filtered;
    ref_times{i, n_actions(i)+1} = tt;
    real_times{i, n_actions(i)+1} = tt;
    mapping_scores{i, n_actions(i)+1}    = zeros(size(tt));
    exploiting_scores{i, n_actions(i)+1} = zeros(size(tt));
    nav_scores{i, n_actions(i)+1}        = zeros(size(tt));
    
    t_start = tt(end);
    real_r_start = xx_real_plan(end, 1:3)';
    real_v_start = xx_real_plan(end, 4:6)';
    filter_r_start = xx_filt_plan(end, 1:3)';
    filter_v_start = xx_filt_plan(end, 4:6)';
    plan_r_start = xx_filt_plan(end, 1:3)';
    plan_v_start = xx_filt_plan(end, 4:6)';
    spacecraft_data_plan = spacecraft_data_real;
end

%Also store results in a full vector

[real_trajectory_total,  tt_all] = concatTrajTimeCells(real_trajectory,  real_times,  n_actions+1);
[ref_trajectory_total,   ~     ] = concatTrajTimeCells(ref_trajectory,   real_times,  n_actions+1);
[filter_trajectory_total,~     ] = concatTrajTimeCells(filter_trajectory, real_times,  n_actions+1);

mapping_score_total    = concatScalarCellsAligned(mapping_scores,   real_times, real_trajectory,  n_actions+1);
exploiting_score_total = concatScalarCellsAligned(exploiting_scores,real_times, real_trajectory,  n_actions+1);
nav_score_total        = concatScalarCellsAligned(nav_scores,       real_times, real_trajectory,  n_actions+1);

P_total = concatCovCellsAligned(P_all, real_times, real_trajectory, n_actions+1);

spacecraft_data_out = spacecraft_data_real;

end

function [M, T] = concatTrajTimeCells(C, Ct, n_actions)
% C:  cell(n_set, max(n_actions)) con traiettorie Nx6
% Ct: cell(n_set, max(n_actions)) con vettori tempo
% n_actions: vettore n_setx1 con # archi per set
% Restituisce: M (Nx6) e T (Nx1) con stesse lunghezze e giunzioni coerenti.

    tolX   = 1e-12;   % toll. continuità stato
    tolT   = 1e-9;    % toll. uguaglianza tempi
    epsT   = 1e-6;    % bump per mantenere T strettamente crescente quando c’è un impulso

    chunksX = {};
    chunksT = {};

    prevX = [];
    prevT = [];

    for i = 1:size(C,1)
        for j = 1:n_actions(i)
            X  = C{i,j};
            tt = Ct{i,j};
            if isempty(X) || isempty(tt), continue; end
            tt = tt(:);
            % (opzionale) sanity: stessa lunghezza
            if size(X,1) ~= numel(tt)
                error('Lunghezze mismatch in cella (%d,%d): size(X,1)=%d, numel(tt)=%d', i,j,size(X,1),numel(tt));
            end

            if ~isempty(chunksX)
                % Giunto con il segmento precedente
                if ~isempty(prevT) && abs(prevT(end) - tt(1)) < tolT
                    % Se lo stato è continuo, togli il duplicato di testa
                    if size(prevX,2)==size(X,2) && all(abs(prevX(end,:) - X(1,:)) < tolX)
                        X  = X(2:end,:);
                        tt = tt(2:end);
                    else
                        % Impulso: tieni entrambi e rendi T strettamente crescente
                        tt(1) = tt(1) + epsT;
                    end
                end
            end
            chunksX{end+1} = X;   %#ok<AGROW>
            chunksT{end+1} = tt;  %#ok<AGROW>
            prevX = X; prevT = tt;
        end
    end

    if isempty(chunksX)
        M = []; T = [];
    else
        M = vertcat(chunksX{:});
        T = vertcat(chunksT{:});
    end
end

function S = concatScalarCellsAligned(Cs, Ct, Cstate, n_actions)
% Cs:   cell con vettori scalari (stessa lunghezza di Ct)
% Ct:   cell con vettori tempo
% Cstate: cell con traiettorie (per decidere se rimuovere il primo campione)
% n_actions: vettore n_setx1

    tolX = 1e-12; tolT = 1e-9;

    chunksS = {};
    prevX = []; prevT = [];

    for i = 1:size(Cs,1)
        for j = 1:n_actions(i)
            s  = Cs{i,j};
            tt = Ct{i,j};
            X  = Cstate{i,j};
            if isempty(s) || isempty(tt) || isempty(X), continue; end
            s  = s(:);
            tt = tt(:);

            if numel(s) ~= numel(tt) || size(X,1) ~= numel(tt)
                error('Mismatch lunghezze in cella (%d,%d).', i,j);
            end

            if ~isempty(chunksS)
                if ~isempty(prevT) && abs(prevT(end) - tt(1)) < tolT
                    if size(prevX,2)==size(X,2) && all(abs(prevX(end,:) - X(1,:)) < tolX)
                        % stato continuo: rimuovi duplicato di testa anche nello score
                        s  = s(2:end);
                        % (il tempo è rimosso nella funzione tempi; qui allineiamo solo s)
                    else
                        % impulso: teniamo entrambi -> niente da rimuovere
                    end
                end
            end

            chunksS{end+1} = s; %#ok<AGROW>
            prevX = X; prevT = tt;
        end
    end

    if isempty(chunksS)
        S = [];
    else
        S = vertcat(chunksS{:});
    end
end

function Ptot = concatCovCellsAligned(Pcell, Ct, Cstate, n_actions)
% Pcell : cell(n_set, max(n_actions)) con covarianze [n x n x Ni]
% Ct    : cell con vettori tempo (Ni x 1)
% Cstate: cell con traiettorie (Ni x m), usate per decidere continuità
% n_actions : vettore (n_set x 1) con # archi per set
%
% Output:
%   Ptot : [n x n x Ntot] concatenazione lungo la 3^ dimensione

    tolX = 1e-12;   % toll. continuità stato
    tolT = 1e-9;    % toll. uguaglianza tempi

    chunks = {};
    prevX = [];
    prevT = [];

    for i = 1:size(Pcell,1)
        for j = 1:n_actions(i)
            P  = Pcell{i,j};
            tt = Ct{i,j};
            X  = Cstate{i,j};

            if isempty(P) || isempty(tt) || isempty(X)
                continue;
            end
            tt = tt(:);

            % sanity checks
            [n1,n2,Ni] = size(P);
            if n1 ~= n2
                error('P_all{%d,%d} non è quadrata: size=%dx%dx%d', i,j,n1,n2,Ni);
            end
            if Ni ~= numel(tt)
                error('Mismatch lunghezze in P_all{%d,%d}: size(P,3)=%d, numel(tt)=%d', i,j,Ni,numel(tt));
            end
            if size(X,1) ~= numel(tt)
                error('Mismatch lunghezze in Cstate{%d,%d}: size(X,1)=%d, numel(tt)=%d', i,j,size(X,1),numel(tt));
            end

            % gestione giunto col segmento precedente
            if ~isempty(chunks)
                if ~isempty(prevT) && abs(prevT(end) - tt(1)) < tolT
                    % stato continuo → rimuovi il primo slice del nuovo segmento
                    if size(prevX,2)==size(X,2) && all(abs(prevX(end,:) - X(1,:)) < tolX)
                        if Ni >= 2
                            P = P(:,:,2:end);
                        else
                            % segmento di lunghezza 1: si elimina del tutto
                            P = [];
                        end
                    else
                        % impulso → tieni tutto (nessuna rimozione)
                    end
                end
            end

            if ~isempty(P)
                chunks{end+1} = P; %#ok<AGROW>
                prevX = X; prevT = tt;
            end
        end
    end

    if isempty(chunks)
        % restituisci array vuoto ben formato
        Ptot = zeros(0,0,0);
    else
        % verifica dimensioni coerenti tra i blocchi
        n = size(chunks{1},1);
        for k = 1:numel(chunks)
            if size(chunks{k},1)~=n || size(chunks{k},2)~=n
                error('Dimensioni non coerenti tra blocchi P: trovato %dx%d al blocco %d (atteso %dx%d).', ...
                      size(chunks{k},1), size(chunks{k},2), k, n, n);
            end
        end
        Ptot = cat(3, chunks{:});
    end
end
