function pruned_tree = pruneTreeById(tree, target_id)
% Estrae il sottoalbero radicato in target_id e rimappa gli ID da 1 in poi

    % Crea mappa da ID → indice nella cella originale
    id_to_index = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(tree)
        id_to_index(tree{i}.id) = i;
    end

    % Ricorsivamente raccoglie tutti gli ID dei discendenti
    function indices = collectSubtreeIndices(current_id)
        idx = id_to_index(current_id);
        node = tree{idx};
        indices = idx;  % includi il nodo corrente
        for c = node.children
            indices = [indices, collectSubtreeIndices(c)];
        end
    end

    % Step 1: Estrai gli indici dei nodi da mantenere
    original_indices = unique(collectSubtreeIndices(target_id));
    original_tree = tree(original_indices);

    % Step 2: Crea una nuova mappa da vecchi ID a nuovi ID
    old_ids = cellfun(@(n) n.id, original_tree);
    new_ids = 1:length(old_ids);
    id_map = containers.Map(old_ids, new_ids);  % old_id → new_id

    % Step 3: Costruisci il nuovo albero con ID rimappati
    pruned_tree = cell(1, length(original_tree));
    for i = 1:length(original_tree)
        node = original_tree{i};

        % Rimappa ID e parent
        node.id = id_map(node.id);
        if node.parent == 0 || ~isKey(id_map, node.parent)
            node.parent = 0;  % nuova radice
        else
            node.parent = id_map(node.parent);
        end

        % Rimappa children
        if ~isempty(node.children)
            node.children = arrayfun(@(cid) id_map(cid), node.children);
        end

        pruned_tree{i} = node;
    end
end
