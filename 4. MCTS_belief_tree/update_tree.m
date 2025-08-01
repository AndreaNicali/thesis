function new_tree = update_tree(tree, id, known_map, gamma)

    pruned_tree = pruneTreeById(tree, id);
    new_tree = {};
    root_obs = pruned_tree{1};

    % aggiorna la known map nel nodo radice
    root_obs.data.data_asteroids.features.known_map_features = known_map;
    root_obs.data.data_asteroids.mapping.known_map = known_map;
    new_tree{1} = root_obs;

    % mappa ID → indice nuovo albero
    id_to_new_idx = containers.Map(root_obs.id, 1);

    count = 2;
    for i = 1:length(pruned_tree)
        node = new_tree{i};
        if strcmp(node.type, 'observation')
            % nodo observation: valutiamo azioni
            for j = 1:length(node.children)
                id_action = node.children(j);
                action_node = pruned_tree{id_action};
                action = action_node.action;
                final_time = action_node.final_time;
                initial_time = action_node.time;

                new_tree{id_action} = action_node;
                new_tree{id_action}.score = 0;

                for k = 1:length(action_node.children)
                    id_child = action_node.children(k);
                    state = action_node.initial_states(:, k);
                    tt = initial_time:100:final_time;

                    [xx, tt] = integrateODE(state, tt);
                    P0 = node.cov;
                    % UKF da node.state → stima finale
                    [P_filtered, filt_time] = navigation_for_MCTS(xx, P0, tt', node.data);
                    P_final = P_filtered(:, :, end);

                    % nuova score function
                    [J_of_t, dJdt, new_scores, new_known_map] = total_score(xx, tt, P0, node.data);

                    %Put updated features and map data in the child node data
                    child_data = node.data;
                    child_data.data_asteroids.features.score = new_scores;
                    child_data.data_asteroids.mapping.known_map = new_known_map;

                    new_tree{id_child} = create_observation_node(id_child, xx(end, :)', P_final, id_action, action_node.depth+1, tt(end), child_data);

                    new_tree{id_child}.score = J_of_t(end);
                    new_tree{id_child}.single_score = J_of_t(end);

                    new_tree = backpropagation(new_tree, new_tree{id_child}, gamma);
                end
            end
        end
    end
end
