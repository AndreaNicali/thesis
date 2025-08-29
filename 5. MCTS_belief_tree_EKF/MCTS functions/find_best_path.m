function [best_path, best_actions, best_final_times] = find_best_path(tree)

    % Mappa ID → indice nella cella
    id_to_index = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(tree)
        id_to_index(tree{i}.id) = i;
    end

    % Trova radice
    root_idx = find(cellfun(@(n) n.parent == 0, tree), 1);
    root_node = tree{root_idx};

    % Inizializza lo score della radice solo se è 'observation'
    if strcmp(root_node.type, 'observation')
        initial_score = root_node.single_score;
    else
        initial_score = 0;
    end

    stack = struct( ...
        'idx', root_idx, ...
        'path', root_node.id, ...
        'actions', [], ...
        'final_times', [], ...
        'score', initial_score ...
    );
    stack = stack(:);

    % Output iniziali
    best_score = -inf;
    best_path = [];
    best_actions = [];
    best_final_times = [];

    while ~isempty(stack)
        current = stack(end);
        stack(end) = [];

        node = tree{current.idx};

        % Nodo foglia
        if isempty(node.children)
            if current.score > best_score
                best_score = current.score;
                best_path = current.path;
                best_actions = current.actions;
                best_final_times = current.final_times;
            end
        else
            for i = 1:length(node.children)
                child_id = node.children(i);
                if isKey(id_to_index, child_id)
                    child_idx = id_to_index(child_id);
                    child_node = tree{child_idx};

                    % Score solo se observation
                    new_score = current.score;
                    if strcmp(child_node.type, 'observation')
                        new_score = new_score + child_node.single_score;
                    end

                    new_path = [current.path, child_id];
                    new_actions = current.actions;
                    new_final_times = current.final_times;

                    % Se il nodo attuale è un 'action', prendi dati
                    if strcmp(node.type, 'action')
                        new_actions = [new_actions, node.action];
                        new_final_times = [new_final_times, node.final_time];
                    end

                    % Aggiungi allo stack
                    stack(end+1) = struct( ...
                        'idx', child_idx, ...
                        'path', new_path, ...
                        'actions', new_actions, ...
                        'final_times', new_final_times, ...
                        'score', new_score ...
                    );
                end
            end
        end
    end

    % Aggiungi tempo iniziale (radice)
    if isfield(tree{root_idx}, 'time')
        best_final_times = [tree{root_idx}.time, best_final_times];
    end
end
