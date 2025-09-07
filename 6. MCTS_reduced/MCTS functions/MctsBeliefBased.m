function tree = MctsBeliefBased(x0, P0, t, max_iterations, spacecraft_data, initial_tree)

%Performs Monte Carlo Tree Search from a initial state and time. The struct
%spacecraft_data contains all the reachability parameters togheter with the
%asteroid discretization and the feature characteristics. In this version
%EKF is used for belief update and reachability analysis is the same for
%all observation nodes to reduce computational effort

%Set options and extract some needed data
% options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
% mass_eros = spacecraft_data.data_asteroids.mass;
% omega_body = spacecraft_data.data_asteroids.omega;
% C20 = spacecraft_data.data_asteroids.C20;
% C22 = spacecraft_data.data_asteroids.C22;

if isempty(initial_tree)
    % Root
    tree{1} = createBeliefNode(1, x0, P0, 0, 0, t, spacecraft_data);
else
    tree = initial_tree;
end

%Set MCTS parameters
node_ids = length(tree)+1;
gamma = 1; %discount factor
iter = 0;
ka = 3; % Progressive widening coeff for action nodes
alpha_a = 0.1; % Progressive widening exp for action nodes
c = sqrt(2)/(3); % Exploration Constant

ko = 3; % Progressive widening coeff for belief nodes
alpha_o = 0.1; % Progressive widening exp for belief nodes

%if known map is empty, perform only one iteration of planning
if ~any(spacecraft_data.data_asteroids.features.known_map_features == 1)
    max_iterations = 1;
end

h = waitbar(0, 'MCTS...');

%Cycle on the full iterations (one iteration starts from the root and ends when a new belief node is created)
for i  = 1:max_iterations

    waitbar(i/max_iterations, h); 

    depth = 0;
    new_node = 0;
    
    %Cycle on each tree branch
    while new_node == 0
        
        %CURRENT PARENT NODE SELECTION
        if depth == 0 %Select root
            id_parent = 1;

            if tree{id_parent}.expanded == 0 %Node is expanded i.e. actions from that node are computed
               [tree] = ComputeActions(tree, id_parent);
            end

            parent_node = tree{id_parent};
            
        else %Parent node is the one chosen at the end of the previous cycle

            if tree{id_parent}.expanded == 0 %Node is expanded i.e. actions from that node are computed
                [tree] = ComputeActions(tree, id_parent);
            end
            parent_node = tree{id_parent};

        end
        
        %PROGRESSIVE WIDENING  (to regolate the node expansion)
        N_children_max = max( [ceil(ka * parent_node.visits ^ alpha_a), 1.1]);
        
        %If an action node can be added, select the best untried action to add to the node
        if length(parent_node.children)+1 < N_children_max
            
            [action, tf, parent_node, num] = SelectAction(parent_node, node_ids);
            t0 = parent_node.time; 
            action_node = createActionNode(node_ids, action, parent_node.id, parent_node.depth+1, t0, tf, num);
            node_ids = node_ids+1;

            %CREATION OF THE NEW OBSERVATION NODE FROM THE NEW ACTION NODE
            [child_node, action_node] = PropagateFromActionNode(parent_node, action_node, node_ids, 1);

            %UPDATE NODE VALUES
            tree{parent_node.id} = parent_node;
            tree{action_node.id} = action_node;
            tree{child_node.id} = child_node;
            
            node_ids = node_ids+1;
            new_node = 1;
            
        else %if i can't expand, choose the best child node from the current parent
           
            depth = depth+1;
            
            %Choose the next node
            [id_action, ~] = nodeSelection(tree, id_parent, c);

            %Choose the next observation node from the previous chosen
            %action
            visits = tree{id_action}.visits;
            N_obs_max = max( [ceil(ko * visits ^ alpha_o), 1.1]); %Progressive Widening
            action_node = tree{id_action};

            if length(action_node.children)+1 < N_obs_max %If an i-th node can be added, select the i-th best action to add the node

                [child_node, action_node] = PropagateFromActionNode(parent_node, action_node, node_ids, 0);
                
                %UPDATE NODE VALUES
                tree{parent_node.id} = parent_node;
                tree{action_node.id} = action_node;
                tree{child_node.id} = child_node;
                
                node_ids = node_ids+1;
                new_node = 1;
            else
                depth = depth+1;

                [id_parent, ~] = nodeSelection(tree, id_action, c);
            end
        end
    end
    
    %BACKPROPAGATION
    tree = backpropagation(tree, child_node, gamma);
    iter = iter + 1;
end

close(h); 

end

%% LOCAL FUNCTIONS

function [tree] = ComputeActions(tree, id_parent)
%This function compute the reachabilty set from a belief node. It is done
%on one belief node and it assign the same action list to all the belief
%nodes "brothers" of the original one, (i.e. all the belief nodes that have
%the same parent of the first).
    
    seq = tree{id_parent}.sequence;
    cand = find(cellfun(@(s) isfield(s,'type') && strcmp(s.type,'observation') && isequal(s.sequence, seq) && s.expanded == 1, tree));

    cand(cand == id_parent) = [];    
    
    if isempty(cand)
        %Start reachability
        [~,~,~,U,J_val,T,~,I] = exploreU(tree{id_parent}.state, tree{id_parent}.cov, tree{id_parent}.time, tree{id_parent}.data);
        actions = [U(:, I==1); J_val(I==1)'; T(I==1)'];
        
    else
        actions = tree{cand(1)}.actions;

    end
    
    tree{id_parent}.actions = actions;
    tree{id_parent}.expanded = 1;
    tree{id_parent}.tried_actions = zeros(1, size(tree{id_parent}.actions, 2));
    
end

function [action, tf, parent_node, ind] = SelectAction(parent_node, node_ids)
%This function selects the best action to add to the tree from the
%reachable set of a belief node. It works in a "best first" way: it selects
%the action with the best score that has not been added yet.

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
    
    parent_node.n_tried_actions = parent_node.n_tried_actions + 1; 
    action = parent_node.actions(1:3, id_actions);
    parent_node.tried_actions_list = [parent_node.tried_actions_list, action];
    tf = parent_node.actions(5, id_actions);
    
    parent_node.tried_actions(id_actions) = 1;
    parent_node.final_times(parent_node.n_tried_actions) = tf;
    parent_node.children = [parent_node.children, node_ids];

end

function [child_node, action_node] = PropagateFromActionNode(parent_node, action_node, node_ids, main)
%This function propagate an action node and generates a new belief node.

    mass_eros = parent_node.data.data_asteroids.mass;
    omega_body = parent_node.data.data_asteroids.omega;
    C20 = parent_node.data.data_asteroids.C20;
    C22 = parent_node.data.data_asteroids.C22;
    sigma_magn = parent_node.data.data_guidance.sigma_magn;
    sigma_align = parent_node.data.data_guidance.sigma_align;
    options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
    parent_data = parent_node.data;
    
    action = action_node.action;
    t0 = action_node.time;
    tf = action_node.final_time;

    P = parent_node.cov;
    if ~issymmetric(P)
        P = (P+P')/2;
    end

    if main
        state = parent_node.state + [zeros(3, 1); action];
    else
        state = mvnrnd(parent_node.state, P(1:6, 1:6))' + [zeros(3, 1); action];
    end

    [~, P_man] = pertThrust(action, sigma_magn, sigma_align);
    P(4:6, 4:6) = P(4:6, 4:6) + P_man;

    %Propagate trajectory
    [tt, xx] = ode78(@(t,x) dynamicsModel(t, x, mass_eros, omega_body, C20, C22), t0:100:tf, state, options);
    
    %Propagate Uncertainties
    [P_filtered, xx_filt, eta_f] = navigationFilter(state', xx, P, tt, parent_data);
    Pf = P_filtered(:, :, end);

    [J_of_t, ~, new_scores, new_known_map] = total_score(xx_filt, tt, Pf, parent_data);
    
    %Put updated features and map data in the child node data
    child_data = parent_data;
    child_data.data_asteroids.features.score = new_scores;
    child_data.data_asteroids.mapping.known_map = new_known_map;
    child_data.data_asteroids.features.known_map_features = parent_data.data_asteroids.features.known_map_features; %known features are not updated inside MCTS because the system can't know a feature is there
    child_data.data_guidance.eta0 = eta_f;

    % Update action node and create a new observation node
    action_node.children = [action_node.children, node_ids];
    action_node.initial_states = [action_node.initial_states, state];
    action_node.n_initial_states = action_node.n_initial_states+1;

    child_node = createBeliefNode(node_ids, xx_filt(end, :)', Pf, action_node.id, action_node.depth+1, tf, child_data);
    child_node.score = J_of_t(end);
    child_node.single_score = J_of_t(end);

    child_node.sequence = [parent_node.sequence, action_node.number];


end

function node = createActionNode(id, action, parent_id, depth, t0, tf, num)

    node = struct( ...
        'type', 'action', ...           %Type of node (action or observation)
        'id', id, ...               % ID node
        'action', action, ...        
        'cov', 0, ...
        'parent', parent_id, ...    % ID parent
        'children', [], ...         % ID children
        'visits', 0, ...            % Visits
        'score', 0, ...             % Total score
        'initial_states', [], ...        
        'n_initial_states', 0, ...      
        'depth', depth, ...         % Profondità del nodo
        'time', t0, ...              % Time of the node
        'final_time', tf,...
        'number', num...
    );

end

function node = createBeliefNode(id, state, P, parent_id, depth, t, spacecraft_data)

    node = struct( ...
        'type', 'observation', ...           %Type of node (action or observation)
        'id', id, ...               % ID node
        'state', state, ...         % State
        'cov', P, ...
        'parent', parent_id, ...    % ID parent
        'children', [], ...         % ID children
        'visits', 0, ...            % Visits
        'score', 0, ...             % Total score
        'actions', [], ...          % Available actions
        'n_tried_actions', 0, ...   % Tried actions 
        'depth', depth, ...         % Profondità del nodo
        'time', t, ...              % Time of the node
        'data', spacecraft_data, ...% Data 
        'expanded', 0, ...          % expanded or not expanded
        'single_score', 0, ...
        'tried_actions', [], ...
        'tried_actions_list', [], ...
        'final_times', [], ...
        'sequence', [] ...
    );

end


function tree = backpropagation(tree, child_node, gamma)

%Perform backpropagation after an iteration is completed

adding_score = child_node.score;
parent_id = child_node.parent;
flag = 1;

while flag
    
    if parent_id == 0 %If the root has been reached
        return
    end
    %Otherwise keep updating parent node
    tree{parent_id}.score = tree{parent_id}.score + gamma*adding_score;
    tree{parent_id}.visits = tree{parent_id}.visits + 1;

    child_node = tree{parent_id};
    parent_id = child_node.parent;
end

end



function [id_node, UCB_score] = nodeSelection(tree, id_parent, c)

%Uses Upper Confidence Bound formula to choose which node to expand. If a
%node hasn't been expanded yet, its value is infinite hence it is always 
%expanded

%Find all the children of the parent node
all_id = cellfun(@(n) n.id, tree);
parent = find(all_id == id_parent);

%For each child, compute UCB
layer = tree{parent}.children;
N_parents = tree{id_parent}.visits;

for i = 1:length(layer)

    if tree{layer(i)}.visits == 0
        UCB_score(i) = inf;

    else
        average_score = tree{layer(i)}.score / tree{layer(i)}.visits;
        exploration = c * sqrt(log(N_parents) / tree{layer(i)}.visits);
        UCB_score(i) = average_score + exploration;

    end

end

%Find max UCB node
[~, id] = max(UCB_score);
id_node = layer(id);

end