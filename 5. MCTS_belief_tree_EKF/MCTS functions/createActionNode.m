function node = createActionNode(id, action, parent_id, depth, t0, tf)

    node = struct( ...
        'type', 'action', ...           %Type of node (action or observation)
        'id', id, ...               % ID node
        'action', action, ...         % State
        'cov', 0, ...
        'parent', parent_id, ...    % ID parent
        'children', [], ...         % ID children
        'visits', 0, ...            % Visits
        'score', 0, ...             % Total score
        'initial_states', [], ...        
        'n_initial_states', 0, ...      
        'depth', depth, ...         % Profondità del nodo
        'time', t0, ...              % Time of the node
        'final_time', tf...
    );




end

