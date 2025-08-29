function [J_of_t, dJdt, new_scores, new_known_map, mapping_score_t, exploit_score_t, nav_score] = total_score(y, t, P0, spacecraft_data)

    alpha1 = 0.5;
    alpha2 = 0.5;
    alpha3 = 1-alpha1-alpha2;

    t = reshape(t, [length(t), 1]);

    [mapping_score_t, exploit_score_t, new_scores, new_known_map, penalty] = score(y, t, spacecraft_data);
    [time] = feature_visibility(y, t, spacecraft_data);
    
    nav_bool = zeros(size(t));
    for i = 1:length(t)
        list = find(time == t(i));
        if length(list)>8
            nav_bool(i) = 1;
        end
    end
    detPMin = -48;

    if log10(det(P0)) < detPMin
        nav_score = 0* nav_bool;
    else
        nav_score = (1 - ( log10(det(P0) ) / detPMin) ) * nav_bool;
    end
    
    nav_score = 0;

    J_of_t = cumsum(alpha1*mapping_score_t + alpha2*exploit_score_t + alpha3*nav_score + penalty);

    dJdt = diff(J_of_t) ./ diff(t);

    dJdt = [dJdt; dJdt(end)];


end