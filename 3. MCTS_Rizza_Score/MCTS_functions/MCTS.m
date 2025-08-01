function tree = MCTS(x0, t,  max_iterations, spacecraft_data)

%Performs Monte Carlo Tree Search from a initial state and time. The struct
%spacecraft_data contains all the reachability parameters togheter with the
%asteroid discretization and the feature characteristics.

%OUTPUT: tree, cell of structs, each struct represent a node of the tree

%Set options and extract needed data
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
mass_eros = spacecraft_data.data_asteroids.mass;
omega_body = spacecraft_data.data_asteroids.omega;
C20 = spacecraft_data.data_asteroids.C20;
C22 = spacecraft_data.data_asteroids.C22;

% Initialize tree
tree = {};

% Root
tree{1} = create_node(1, x0,0, 0, 0, t, spacecraft_data);

%MCTS parameters
gamma = 0.2; %discount factor
max_depth = 100; %maximum depth (not important since the chosen ending condition is the number of iterations)
ka = 3; %Data for progressive widening
alpha_a = 0.1;
c = sqrt(2)/100; %Exploration constant

%Initialize parameters
node_idx = 2;
iter = 0;


h = waitbar(0, 'MCTS...');

%Cicle on the full iteration. A full iteration starts from the root and
%ends when an unexpanded node is reached

for i  = 1:max_iterations

    depth = 0;
    new_node = 0;
    
    %Cycle on each tree branch
    while depth < max_depth && new_node == 0
        
        %CURRENT PARENT NODE SELECTION
        if depth == 0 %Select root
            id_parent = 1;

            if tree{id_parent}.expanded == 0 %If not already expanded, node is expanded i.e. actions from that node are computed
                [~,~,~,U,J_val,T,S,I] = exploreU(tree{id_parent}.state, tree{id_parent}.time, tree{id_parent}.data);
                tree{id_parent}.actions = [U(:, I==1); J_val(I==1)'; T(I==1)'];
                tree{id_parent}.expanded = 1;

            end

            parent_node = tree{id_parent};
            parent_data = parent_node.data;

            
        else %If it's not the root, parent node is the one chosen at the end of the previous cycle

            if tree{id_parent}.expanded == 0 %If not already expanded, node is expanded i.e. actions from that node are computed
                [~,~,~,U,J_val,T,S,I] = exploreU(tree{id_parent}.state, tree{id_parent}.time, tree{id_parent}.data);
                tree{id_parent}.actions = [U(:, I==1); J_val(I==1)'; T(I==1)'];
                tree{id_parent}.expanded = 1;

            end
            parent_node = tree{id_parent};
            parent_data = parent_node.data;

        end
        
        %CHILD SELECTION
        %Progressive widening  (to regolate the node expansion)
        N_children_max = max( [ceil(ka * parent_node.visits ^ alpha_a), 1.1]);

        if length(parent_node.children)+1 < N_children_max %If an i-th node can be added, select the i-th best action to add the node
            parent_node.n_tried_actions = parent_node.n_tried_actions + 1; 
            [~, id_actions] = maxk(parent_node.actions(4, :), parent_node.n_tried_actions);
            t0 = parent_node.time;
            
            if isscalar(id_actions)
                action = parent_node.actions(1:3, id_actions);
                tf = parent_node.actions(5, id_actions);

                parent_node.tried_actions(:, parent_node.n_tried_actions) = action;
                parent_node.final_times(parent_node.n_tried_actions) = tf;

            else
                action = parent_node.actions(1:3, id_actions(end));
                tf = parent_node.actions(5, id_actions(end));

                parent_node.tried_actions(:, parent_node.n_tried_actions) = action;
                parent_node.final_times(parent_node.n_tried_actions) = tf;

            end

            %CREATION OF THE NEW NODE FROM THE NEW ACTION
            state = parent_node.state + [zeros(3, 1); action];

            [tt, xx] = ode78(@(t,x) dynamicsEllipsoid(t, x, mass_eros, omega_body, C20, C22), t0:100:tf, state, options);
            [J, ~, features_post, known_map_post] = score_rizza(xx, tt', parent_data); %Obtain updated features and map data
            
            features_post.known = parent_data.data_asteroids.features.known;

            %Put updated features and map data in the child node data
            child_data = parent_data;
            child_data.data_asteroids.features = features_post;
            child_data.data_asteroids.known_map = known_map_post;
            
            %Create children node with updated data
            tree{node_idx} = create_node(node_idx, xx(end, :)', 0, parent_node.id, parent_node.depth+1, tf, child_data);
            
            tree{node_idx}.score = J(end);
            tree{node_idx}.single_score = J(end);

            %UPDATE PARENT NODE VALUES

            parent_node.children = [parent_node.children, node_idx];
            tree{id_parent} = parent_node;
            child_node = tree{node_idx};

            node_idx = node_idx + 1;
            new_node = 1;
            
        else %if i can't expand, choose the best child node from the current parent
           
            depth = depth+1;
            
            %Choose the next node
            [id_parent, ~] = node_selection(tree, id_parent, c);

        end
    end
    
    %BACKPROPAGATION of data up to the tree
    tree = backpropagation(tree, child_node, gamma);
    iter = iter + 1;

    waitbar(i/max_iterations, h); 

end

close(h); 

end