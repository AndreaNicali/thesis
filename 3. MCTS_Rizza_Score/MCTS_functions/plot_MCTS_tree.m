function plot_MCTS_tree(tree)
    % Estrai gli ID dei nodi e le relazioni genitore-figlio
    edges = [];

    for i = 1:length(tree)
        node = tree{i};
        for j = 1:length(node.children)
            child_id = node.children(j);
            edges = [edges; node.id, child_id];
        end
    end

    % Crea grafo orientato
    G = digraph(edges(:,1), edges(:,2));

    % Disegna il grafo
    figure;
    h = plot(G, 'Layout', 'layered', 'Direction', 'down', 'MarkerSize', 7, ...
        'NodeColor', 'c', 'EdgeAlpha', 0.6, 'LineWidth', 1.5);

    % Aggiungi etichette ai nodi (ID o altri dati se vuoi)
    labelnode(h, 1:numnodes(G), string(1:numnodes(G)));

    title('MCTS Tree');
    axis off;
end
