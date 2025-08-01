function node = create_observation_node(id, state, P, parent_id, depth, t, spacecraft_data)

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
        'final_times', [] ...
    );




end


