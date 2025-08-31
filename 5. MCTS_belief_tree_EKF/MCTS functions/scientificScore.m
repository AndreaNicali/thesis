function [mapping_score_t, exploit_score_t, new_scores, new_known_map, penalty_of_t] = scientificScore(y, t, spacecraft_data)

% Estrai dati di guida e sensore
data = spacecraft_data.data_guidance;
r_impact = data.r_impact;
r_escape = data.r_escape;
fov1 = data.scientificFov1;
fov2 = data.scientificFov2;

% Estrai geometria dell’asteroide
F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;

% Mappatura e feature
mapping = spacecraft_data.data_asteroids.mapping;
scores = spacecraft_data.data_asteroids.features.score;
known_map_features = spacecraft_data.data_asteroids.features.known_map_features;
known_map = spacecraft_data.data_asteroids.mapping.known_map;

% Inizializza mappe aggiornate
new_known_map = known_map;
new_scores = scores;
single_score_mapping = 1 / length(known_map);

% Preallocazione
r = y(:, 1:3);
n_time = size(y, 1);
n_faces = size(F, 1);

mapping_score_t = zeros(n_time, 1);
exploit_score_t = zeros(n_time, 1);
penalty_of_t = zeros(n_time, 1);
r_sun_vec = zeros(n_time, 3);
R_body_cam_all = zeros(3, 3, n_time);
R_body2topo_all = zeros(3, 3, n_faces);
centroid_vec = zeros(n_faces, 3);

% Pre-calcolo: posizione Sole e rotazione corpo-camera
for i = 1:n_time
    r_sun = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    r_sun = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i)) * (-r_sun(1:3));
    r_sun_vec(i, :) = r_sun;
    R_body_cam_all(:, :, i) = body2camera(y(i, 1:3), y(i, 4:6));
end

% Pre-calcolo: centroidi e sistemi topo-centrici
for j = 1:n_faces
    centroid = mean(V(F(j, :), :), 1);
    centroid_vec(j, :) = centroid;
    Z = N(:, j);
    X = [-centroid(2); centroid(1); 0];
    X = X / norm(X);
    Y = cross(Z, X);
    R_body2topo_all(:, :, j) = [X'; Y'; Z'];
end

% Ciclo temporale
for i = 1:n_time
    r_now = r(i, :);
    R_cam = R_body_cam_all(:, :, i);
    r_sun = r_sun_vec(i, :);

    for j = 1:n_faces
        centroid = centroid_vec(j, :);
        r_relative = r_now - centroid;
        norm_rrel = norm(r_relative);
        normal_j = N(:, j);

        % === EXPLOITATION ===
        if known_map_features(j) == 1 && new_scores(j).score > 0
            for k = 1:length(new_scores(j).completeness)
                if new_scores(j).completeness(k) ~= 1
                    is_emission = strcmp(new_scores(j).type, 'emission');
                    is_relative = strcmp(new_scores(j).type, 'relative');

                    if is_emission
                        angle_emission = acos(dot(r_relative, normal_j) / (norm_rrel * norm(normal_j)));
                        in_range = angle_emission > new_scores(j).angles(k, 1) && angle_emission < new_scores(j).angles(k, 2);

                        if in_range
                            inc = dot(r_sun, normal_j) / (norm(r_sun) * norm(normal_j));
                            if inc > 0
                                % === CHECK FOV ===
                                inside1 = check_FOV(V(F(j, 1), :), r_now, fov1, fov2, R_cam);
                                if inside1
                                    inside2 = check_FOV(V(F(j, 2), :), r_now, fov1, fov2, R_cam);
                                    if inside2
                                        inside3 = check_FOV(V(F(j, 3), :), r_now, fov1, fov2, R_cam);
                                        if inside3
                                            if isnan(new_scores(j).actual_range(k)) || abs(norm_rrel-new_scores(j).ideal_range) < abs(new_scores(j).actual_range(k)-new_scores(j).ideal_range)
                                                new_scores(j).actual_range(k) = norm_rrel;
                                                space = 2;
                                                if norm_rrel < new_scores(j).ideal_range+space && norm_rrel > new_scores(j).ideal_range-space
                                                    new_scores(j).completeness(k) = 1;
                                                else 
                                                    if norm_rrel < new_scores(j).ideal_range
                                                        new_scores(j).completeness(k) = exp( -(norm_rrel - (new_scores(j).ideal_range - space))^2 / 2*0.5 );
                                                    else
                                                        new_scores(j).completeness(k) = exp( -(norm_rrel - (new_scores(j).ideal_range + space))^2 / 2*0.5 );
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    elseif is_relative
                        angle_rel = acos(dot(r_relative, r_sun) / (norm_rrel * norm(r_sun)));
                        in_range = angle_rel > new_scores(j).angles(k, 1) && angle_rel < new_scores(j).angles(k, 2);

                        if in_range
                            emission_dot = dot(r_relative, normal_j) / (norm_rrel * norm(normal_j));
                            inc = dot(r_sun, normal_j) / (norm(r_sun) * norm(normal_j));

                            if emission_dot > 0 && inc > 0
                                %=== CHECK FOV  ===
                                inside1 = check_FOV(V(F(j, 1), :), r_now, fov1, fov2, R_cam);
                                if inside1
                                    inside2 = check_FOV(V(F(j, 2), :), r_now, fov1, fov2, R_cam);
                                    if inside2
                                        inside3 = check_FOV(V(F(j, 3), :), r_now, fov1, fov2, R_cam);
                                        if inside3
                                            if isnan(new_scores(j).actual_range(k)) || abs(norm_rrel-new_scores(j).ideal_range) < abs(new_scores(j).actual_range(k)-new_scores(j).ideal_range)
                                                new_scores(j).actual_range(k) = norm_rrel;
                                                space = 2;
                                                if norm_rrel < new_scores(j).ideal_range+space && norm_rrel > new_scores(j).ideal_range-space
                                                    new_scores(j).completeness(k) = 1;
                                                else 
                                                    if norm_rrel < new_scores(j).ideal_range
                                                        new_scores(j).completeness(k) = exp( -(norm_rrel - (new_scores(j).ideal_range - space))^2 / 2*0.5 );
                                                    else
                                                        new_scores(j).completeness(k) = exp( -(norm_rrel - (new_scores(j).ideal_range + space))^2 / 2*0.5 );
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

        % === ZONE DISCOVERY ===
        elseif known_map(j) == 0
            inc_angle = acos(dot(r_sun, normal_j) / (norm(r_sun) * norm(normal_j)));
            if inc_angle > mapping.incidence(1) && inc_angle < mapping.incidence(2)
                angle_emission = acos(dot(r_relative, normal_j) / (norm_rrel * norm(normal_j)));
                if angle_emission > mapping.emission(1) && angle_emission < mapping.emission(2)
                    if norm_rrel < 50
                    if norm_rrel > 0
                        inside1 = check_FOV(V(F(j, 1), :), r_now, fov1, fov2, R_cam);
                        if inside1
                            inside2 = check_FOV(V(F(j, 2), :), r_now, fov1, fov2, R_cam);
                            if inside2
                                inside3 = check_FOV(V(F(j, 3), :), r_now, fov1, fov2, R_cam);
                                if inside3
                                    new_known_map(j) = 1;
                                end
                            end
                        end
                    end
                    end
                end
            end
        end
    end

    % === MAPPING SCORE ===
    if i == 1
        mapping_score_t(i) = sum(new_known_map - known_map) * single_score_mapping;
        old_known_map = new_known_map;
    else
        mapping_score_t(i) = sum(new_known_map - old_known_map) * single_score_mapping;
        old_known_map = new_known_map;
    end

    % === PENALTY ===
    norm_r = norm(r_now);
    if norm_r < r_impact || norm_r > r_escape
        penalty_of_t(i) = -5000;
    end

    % === EXPLOIT SCORE ===
    if i == 1
        old_scores = scores;
    end

    for j = 1:n_faces
        for k = 1:length(new_scores(j).completeness)
            if new_scores(j).completeness(k)
                delta = new_scores(j).completeness(k) - old_scores(j).completeness(k);
                s = delta * old_scores(j).score / length(old_scores(j).completeness);
                exploit_score_t(i) = exploit_score_t(i) + s;
            end
        end
    end
    old_scores = new_scores;
end

end
