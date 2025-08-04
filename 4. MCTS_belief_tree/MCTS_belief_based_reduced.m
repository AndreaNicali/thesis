function tree = MCTS_belief_based_reduced(x0, P0, t, max_iterations, spacecraft_data, initial_tree)

%Performs Monte Carlo Tree Search from a initial state and time. The struct
%spacecraft_data contains all the reachability parameters togheter with the
%asteroid discretization and the feature characteristics.

%Set options and extract some needed data
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
mass_eros = spacecraft_data.data_asteroids.mass;
omega_body = spacecraft_data.data_asteroids.omega;
C20 = spacecraft_data.data_asteroids.C20;
C22 = spacecraft_data.data_asteroids.C22;

if length(initial_tree) == 0
    % Root
    tree{1} = create_observation_node(1, x0, P0, 0, 0, t, spacecraft_data);
else
    tree = initial_tree;
end

%MCTS parameters
node_ids = length(tree)+1;
gamma = 0.2;
max_depth = 100;
iter = 0;
ka = 3;
alpha_a = 0.1;
c = sqrt(2)/(2);

ko = 3;
alpha_o = 0.1;

if ~any(spacecraft_data.data_asteroids.features.known_map_features == 1)
    max_iterations = 1;
end

h = waitbar(0, 'MCTS...');

%Cicle on the full iterations (each starting from root)
for i  = 1:max_iterations

    waitbar(i/max_iterations, h); 

    depth = 0;
    new_node = 0;
    
    %Cycle on each tree branch
    while depth < max_depth && new_node == 0
        
        %CURRENT PARENT NODE SELECTION
        if depth == 0 %Select root
            id_parent = 1;

            if tree{id_parent}.expanded == 0 %Node is expanded i.e. actions from that node are computed
                [~,~,~,U,J_val,T,S,I] = exploreU(tree{id_parent}.state,tree{id_parent}.cov, tree{id_parent}.time, tree{id_parent}.data);
                tree{id_parent}.actions = [U(:, I==1); J_val(I==1)'; T(I==1)'];
                tree{id_parent}.expanded = 1;
                tree{id_parent}.tried_actions = zeros(1, size(tree{id_parent}.actions, 2));

            end

            parent_node = tree{id_parent};
            parent_data = parent_node.data;

            
        else %Parent node is the one chosen at the end of the previous cycle

            if tree{id_parent}.expanded == 0 %Node is expanded i.e. actions from that node are computed
                [~,~,~,U,J_val,T,S,I] = exploreU(tree{id_parent}.state, tree{id_parent}.cov, tree{id_parent}.time, tree{id_parent}.data);
                tree{id_parent}.actions = [U(:, I==1); J_val(I==1)'; T(I==1)'];
                tree{id_parent}.expanded = 1;
                tree{id_parent}.tried_actions = zeros(1, size(tree{id_parent}.actions, 2));

                bro_idx = cellfun(@(s) s.parent == tree{id_parent}.parent, tree);
                if sum(bro_idx) > 1
                    for j = find(bro_idx)
                        tree{j}.actions = [U(:, I==1); J_val(I==1)'; T(I==1)'];
                        tree{j}.expanded = 1;
                        tree{j}.tried_actions = zeros(1, size(tree{j}.actions, 2));
                    end
                end

            end
            parent_node = tree{id_parent};
            parent_data = parent_node.data;

        end
        
        %PROGRESSIVE WIDENING  (to regolate the node expansion)
        N_children_max = max( [ceil(ka * parent_node.visits ^ alpha_a), 1.1]);
        
        %If an i-th node can be added, select the i-th best action to add the node
        if length(parent_node.children)+1 < N_children_max
            parent_node.n_tried_actions = parent_node.n_tried_actions + 1; 
            untried_not_found = 1;
            ind = 1;

            while untried_not_found
                list = parent_node.actions(4, :);
                [~, id_actions] = maxk(list, ind);
                if ~isscalar(id_actions)
                    id_actions = id_actions(end);
                end

                if parent_node.tried_actions(id_actions)
                    ind = ind+1;
                else
                    untried_not_found = 0;
                end
            end

            t0 = parent_node.time;
            
            action = parent_node.actions(1:3, id_actions);
            parent_node.tried_actions_list = [parent_node.tried_actions_list, action];
            tf = parent_node.actions(5, id_actions);

            parent_node.tried_actions(id_actions) = 1;
            parent_node.final_times(parent_node.n_tried_actions) = tf;
            parent_node.children = [parent_node.children, node_ids];

            action_node = create_action_node(node_ids, action, parent_node.id, parent_node.depth+1, t0, tf);
            node_ids = node_ids+1;

            %CREATION OF THE NEW OBSERVATION NODE FROM THE NEW ACTION NODE
            P = parent_node.cov;
            if ~issymmetric(P)
                P = (P+P')/2;
            end
            state = mvnrnd(parent_node.state, P)' + [zeros(3, 1); action];
            
            %Propagate trajectory
            [tt, xx] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), t0:100:tf, state, options);
            
            %Propagate Uncertainties
            [P_filtered, filt_time] = navigation_for_MCTS(xx, P, tt, parent_data);
            Pf = P_filtered(:, :, end);

            [J_of_t, dJdt, new_scores, new_known_map] = total_score(xx, tt, Pf, parent_data);
            
            %Put updated features and map data in the child node data
            child_data = parent_data;
            child_data.data_asteroids.features.score = new_scores;
            child_data.data_asteroids.mapping.known_map = new_known_map;
            child_data.data_asteroids.features.known_map_features = parent_data.data_asteroids.features.known_map_features; %known features are not updated inside MCTS because the system can't know a feature is there
            
            % Update action node and create a new observation node
            action_node.children = [action_node.children, node_ids];
            action_node.initial_states = [action_node.initial_states, state];
            action_node.n_initial_states = action_node.n_initial_states+1;

            child_node = create_observation_node(node_ids, xx(end, :)', Pf, action_node.id, action_node.depth+1, tf, child_data);
            child_node.score = J_of_t(end);
            child_node.single_score = J_of_t(end);

            %UPDATE NODE VALUES
            tree{parent_node.id} = parent_node;
            tree{action_node.id} = action_node;
            tree{child_node.id} = child_node;
            
            node_ids = node_ids+1;
            new_node = 1;
            
        else %if i can't expand, choose the best child node from the current parent
           
            depth = depth+1;
            
            %Choose the next node
            [id_action, ~] = node_selection(tree, id_parent, c);

            %Choose the next observation node from the previous chosen
            %action
            visits = tree{id_action}.visits;
            N_obs_max = max( [ceil(ko * visits ^ alpha_o), 1.1]);
            action = tree{id_action}.action;
            P = parent_node.cov;

            action_node = tree{id_action};

            if visits+1 < N_obs_max %If an i-th node can be added, select the i-th best action to add the node
                if ~issymmetric(P)
                   s = 1;
                end
                state = mvnrnd(parent_node.state, P)' + [zeros(3, 1); action];
                t0 = action_node.time;
                tf = action_node.final_time;

                %Propagate trajectory
                [tt, xx] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), t0:100:tf, state, options);
                
                %Propagate Uncertainties
                [P_filtered, filt_time] = navigation_for_MCTS(xx, P, tt, parent_data);
                Pf = P_filtered(:, :, end);
    
                [J_of_t, dJdt, new_scores, new_known_map] = total_score(xx, tt, Pf, parent_data);
                
                %Put updated features and map data in the child node data
                child_data = parent_data;
                child_data.data_asteroids.features.score = new_scores;
                child_data.data_asteroids.mapping.known_map = new_known_map;
                child_data.data_asteroids.features.known_map_features = parent_data.data_asteroids.features.known_map_features; %known features are not updated inside MCTS because the system can't know a feature is there
                
                % Update action node and create a new observation node
                action_node.children = [action_node.children, node_ids];
                action_node.initial_states = [action_node.initial_states, state];
                action_node.n_initial_states = action_node.n_initial_states+1;
    
                child_node = create_observation_node(node_ids, xx(end, :)', Pf, action_node.id, action_node.depth+1, tf, child_data);
                child_node.score = J_of_t(end);
                child_node.single_score = J_of_t(end);
                
    
                %UPDATE NODE VALUES
                tree{parent_node.id} = parent_node;
                tree{action_node.id} = action_node;
                tree{child_node.id} = child_node;
                
                node_ids = node_ids+1;
                new_node = 1;
            else
                depth = depth+1;

                [id_parent, ~] = node_selection(tree, id_action, c);
            end
        end
    end
    
    %BACKPROPAGATION
    tree = backpropagation(tree, child_node, gamma);
    iter = iter + 1;
end

close(h); 

end
