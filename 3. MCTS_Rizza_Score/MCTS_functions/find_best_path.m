function [best_path, best_actions, best_final_times] = find_best_path(tree)
    % Mappa da ID a indice della cella
    id_to_index = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(tree)
        id = double(tree{i}.id);
        id_to_index(id) = i;
    end

    % Trova radice (parent == 0)
    root_idx = find(cellfun(@(n) n.parent == 0, tree), 1);

    % Stack per DFS: struct con percorso, azioni, tempi e score
    stack = struct( ...
        'idx', root_idx, ...
        'path', tree{root_idx}.id, ...
        'actions', [], ...           % matrice 3×N
        'final_times', [], ...
        'score', tree{root_idx}.single_score ...
    );
    stack = stack(:);

    % Risultato
    best_path = [];
    best_actions = [];
    best_final_times = [];
    best_score = -inf;

    while ~isempty(stack)
        current = stack(end);
        stack(end) = [];

        node = tree{current.idx};

        if isempty(node.children)
            if current.score > best_score
                best_score = current.score;
                best_path = current.path;
                best_actions = current.actions;
                best_final_times = current.final_times;
            end
        else
            for i = 1:length(node.children)
                child_id = double(node.children(i));
                if isKey(id_to_index, child_id)
                    child_idx = id_to_index(child_id);
                    child_node = tree{child_idx};

                    new_path = [current.path, child_id];
                    new_score = current.score + child_node.single_score;

                    % Azione (colonna 3x1)
                    action_vec = node.tried_actions(:, i);
                    new_actions = [current.actions, action_vec];

                    % Tempo finale
                    new_times = [current.final_times, node.final_times(i)];

                    % Push nel stack
                    stack(end+1) = struct( ...
                        'idx', child_idx, ...
                        'path', new_path, ...
                        'actions', new_actions, ...
                        'final_times', new_times, ...
                        'score', new_score ...
                    );
                end
            end
        end
    end
    best_final_times = [tree{1}.time, best_final_times];
end
