function plot_MCTS_tree(tree)
    % Estrai gli ID dei nodi e le relazioni genitore-figlio
    edges = [];
    node_types = containers.Map('KeyType', 'double', 'ValueType', 'char');

    for i = 1:length(tree)
        node = tree{i};
        node_types(node.id) = node.type;

        for j = 1:length(node.children)
            child_id = node.children(j);
            edges = [edges; node.id, child_id];
        end
    end

    % Assegna nomi espliciti ai nodi
    node_ids = unique(edges(:));
    node_names = string(node_ids);
    G = digraph(edges(:,1), edges(:,2), [], table(node_names, 'VariableNames', {'Name'}));

    % Disegna il grafo
    figure;
    h = plot(G, 'Layout', 'layered', 'Direction', 'down', ...
        'NodeLabel', {}, 'Marker', 'o', 'NodeColor', [0.7 0.7 0.7]);

    % Categorie dei nodi
    obs_nodes = [];
    act_nodes = [];

    for i = 1:numnodes(G)
        id = str2double(G.Nodes.Name{i});
        if isKey(node_types, id)
            if strcmp(node_types(id), 'observation')
                obs_nodes(end+1) = i;
            elseif strcmp(node_types(id), 'action')
                act_nodes(end+1) = i;
            end
        end
    end

    % Evidenzia tipi
    highlight(h, obs_nodes, 'NodeColor', 'b', 'Marker', 's', 'MarkerSize', 8);
    highlight(h, act_nodes, 'NodeColor', 'r', 'Marker', 'o', 'MarkerSize', 8);

    % Etichette
    labelnode(h, 1:numnodes(G), G.Nodes.Name);

    title('MCTS Tree');
    axis off;
end
