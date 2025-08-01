function [id_node, UCB_score] = node_selection(tree, id_parent, c)

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