function tree = backpropagation(tree, child_node, gamma)

%Perform backpropagation after an iteration is completed

adding_score = child_node.score;
parent_id = child_node.parent;
flag = 1;

while flag
    
    if parent_id == 0 %If the root has been reached
        return

    end
    %Otherwise keep updating parent node
    tree{parent_id}.score = tree{parent_id}.score + gamma*adding_score;
    tree{parent_id}.visits = tree{parent_id}.visits + 1;

    child_node = tree{parent_id};
    parent_id = child_node.parent;
end



