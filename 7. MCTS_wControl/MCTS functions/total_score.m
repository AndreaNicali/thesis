function [J_of_t, dJdt, new_scores, new_known_map, mapping_score_t, exploit_score_t, nav_score] = total_score(y, t, P0, spacecraft_data)

    alpha1 = spacecraft_data.data_guidance.alpha(1);
    alpha2 = spacecraft_data.data_guidance.alpha(2);
    alpha3 = spacecraft_data.data_guidance.alpha(3);

    t = reshape(t, [length(t), 1]);
    
    %Scientific score contains both mapping and expoitation score
    [mapping_score_t, exploit_score_t, new_scores, new_known_map, penalty] = scientificScore(y, t, spacecraft_data);
    [nav_score] = navigationScore(y, P0, t, spacecraft_data);
    
    J_of_t = cumsum(alpha1*mapping_score_t + alpha2*exploit_score_t + alpha3*nav_score + penalty);

    dJdt = diff(J_of_t) ./ diff(t);

    dJdt = [dJdt; dJdt(end)];


end