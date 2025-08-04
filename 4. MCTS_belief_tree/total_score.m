function [J_of_t, dJdt, new_scores, new_known_map, mapping_score_t, exploit_score_t, nav_score] = total_score(y, t, P0, spacecraft_data)

    alpha1 = 1;
    alpha2 = 1;
    alpha3 = 1;

    t = reshape(t, [length(t), 1]);

    [mapping_score_t, exploit_score_t, new_scores, new_known_map, penalty] = score(y, t, spacecraft_data);
    [time] = feature_visibility(y, t, spacecraft_data);
    
    nav_bool = zeros(size(t));
    for i = 1:length(t)
        list = find(time == t(i));
        if length(list)>2
            nav_bool(i) = 1;
        end
    end
    detPMin = -51;

    nav_score = (1 - ( log10(det(P0) ) / detPMin) ) * nav_bool;

    if any( size(mapping_score_t) ~= size(exploit_score_t)) || any (size(mapping_score_t) ~= size(nav_score)) || any(size(mapping_score_t) ~= size(penalty))
        a = 0;
    end

    J_of_t = cumsum(alpha1*mapping_score_t + alpha2*exploit_score_t + alpha3*nav_score + penalty);

    dJdt = diff(J_of_t) ./ diff(t);

    dJdt = [dJdt; dJdt(end)];


end