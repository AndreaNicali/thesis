function [best_path, best_actions, best_action_times, best_final_times] = find_best_path(tree)

    % Mappa ID -> indice in "tree"
    id_to_index = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(tree)
        id_to_index(tree{i}.id) = i;
    end

    % Trova la radice
    root_idx  = find(cellfun(@(n) n.parent == 0, tree), 1);
    root_node = tree{root_idx};

    % Score iniziale (solo se radice 'observation')
    if strcmp(root_node.type, 'observation')
        initial_score = root_node.single_score;
    else
        initial_score = 0;
    end

    % Stack per DFS: actions / action_times / final_times sono CELLE
    stack = struct( ...
        'idx', root_idx, ...
        'path', root_node.id, ...
        'actions', {{}}, ...
        'action_times', {{}}, ...
        'final_times', {{}}, ...
        'score', initial_score ...
    );
    stack = stack(:);

    % Migliori risultati
    best_score        = -inf;
    best_path         = [];
    best_actions      = {};
    best_action_times = {};
    best_final_times  = {};

    while ~isempty(stack)
        current   = stack(end);
        stack(end) = [];

        node = tree{current.idx};

        % Nodo foglia
        if isempty(node.children)
            if current.score > best_score
                best_score        = current.score;
                best_path         = current.path;
                best_actions      = current.actions;       % cell 1xN
                best_action_times = current.action_times;  % cell 1xN
                best_final_times  = current.final_times;   % cell 1xN
            end
        else
            % Espandi i figli
            for i = 1:length(node.children)
                child_id = node.children(i);
                if isKey(id_to_index, child_id)
                    child_idx  = id_to_index(child_id);
                    child_node = tree{child_idx};

                    % Aggiorna score (solo se figlio 'observation')
                    new_score = current.score;
                    if strcmp(child_node.type, 'observation')
                        new_score = new_score + child_node.single_score;
                    end

                    % Propaga path / azioni / tempi
                    new_path         = [current.path, child_id];
                    new_actions      = current.actions;
                    new_action_times = current.action_times;
                    new_final_times  = current.final_times;

                    % Se il NODO ATTUALE è 'action', appendo tempi di inizio/fine
                    if strcmp(node.type, 'action')
                        % Azione
                        new_actions{end+1} = node.action; %#ok<AGROW>

                        % Tempo di INIZIO azione:
                        % - se il nodo 'action' ha un campo 'time' lo uso
                        % - altrimenti, se la radice ha 'time' e non abbiamo ancora tempi, userò root.time dopo
                        if isfield(node, 'time')
                            new_action_times{end+1} = node.time; %#ok<AGROW>
                        else
                            new_action_times{end+1} = []; %#ok<AGROW>
                        end

                        % Tempo di FINE azione (richiesto)
                        if isfield(node, 'final_time')
                            new_final_times{end+1} = node.final_time; %#ok<AGROW>
                        else
                            new_final_times{end+1} = []; %#ok<AGROW>
                        end
                    end

                    % Push sullo stack
                    stack(end+1) = struct( ... %#ok<AGROW>
                        'idx', child_idx, ...
                        'path', new_path, ...
                        'actions', {new_actions}, ...
                        'action_times', {new_action_times}, ...
                        'final_times', {new_final_times}, ...
                        'score', new_score ...
                    );
                end
            end
        end
    end

    % Allineamento tempi secondo la tua definizione:
    % - best_action_times: include l'istante iniziale dell'albero, NON l'ultimo istante finale
    % - best_final_times : include gli istanti finali delle azioni, NON l'istante iniziale
    if isfield(root_node, 'time')
        % Se il primo tempo d'azione non è la radice (o è vuoto), prepend root.time
        prepend_root = true;
        if ~isempty(best_action_times) && ~isempty(best_action_times{1}) ...
                && isnumeric(best_action_times{1}) && isnumeric(root_node.time)
            prepend_root = abs(best_action_times{1} - root_node.time) > 1e-12;
        end
        if prepend_root
            best_action_times = [{root_node.time}, best_action_times];
        end
    end

    % Rimuovi eventuali placeholder vuoti rimasti (se qualche nodo non aveva campi tempo)
    best_action_times = best_action_times(~cellfun(@isempty, best_action_times));
    best_final_times  = best_final_times(~cellfun(@isempty, best_final_times));
end
