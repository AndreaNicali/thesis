function [mapping_score_t, exploit_score_t, new_scores, new_known_map, penalty_of_t] = score(y, t, spacecraft_data)

%Extract data from spacecraft struct
r_impact = spacecraft_data.data_guidance.r_impact;
r_escape = spacecraft_data.data_guidance.r_escape;
fov1 = spacecraft_data.data_guidance.fov1;
fov2 = spacecraft_data.data_guidance.fov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;

known_map = spacecraft_data.data_asteroids.mapping.known_map;
new_known_map = known_map;
mapping = spacecraft_data.data_asteroids.mapping;

scores = spacecraft_data.data_asteroids.features.score;
new_scores = scores;
known_map_features = spacecraft_data.data_asteroids.features.known_map_features;

% SET UP RELATIVE SCORE
single_score_mapping = 1/length(known_map);

%Initialize variables
r = y(:, 1:3);

obs_emission = zeros(size(F, 1), size(y, 1));
sun_incidence = zeros(size(F, 1), size(y, 1));
obs_relative = zeros(size(F, 1), size(y, 1));

exploit_score_t = zeros(size(t));
mapping_score_t = zeros(size(t));
penalty_of_t = zeros(size(t));

r_sun_vec = zeros(size(y,1), 3);
R_body_cam_all = zeros(3, 3, size(y,1));
R_body2topo_all = zeros(3, 3, size(F,1));

%Extract rotation matrices before
for i = 1:size(y, 1)
    r_sun = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    r_sun = r_sun(1:3);
    eclip2eros = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i));
    r_sun = eclip2eros*(-r_sun);
    r_sun_vec(i, :) = r_sun;
    R_body_cam_all(:, :, i) = body2camera(y(i,1:3), y(i,4:6));
end

centroid_vec = zeros(size(F, 1), 3);
for j = 1:size(F, 1)
        centroid_vec(j, :) = ( V(F(j, 1), :) + V(F(j, 2), :) +  V(F(j, 3), :) )/3;
        Z_topo = N(:, j);
        X_topo = [-centroid(2); centroid(1); 0];
        X_topo = X_topo/norm(X_topo);
        Y_topo = cross(Z_topo, X_topo);
        R_body2topo = [X_topo'; Y_topo'; Z_topo'];
        R_body2topo_all(:, :, j) = R_body2topo;
end

for i = 1:size(y, 1)

    for j = 1:size(F, 1)
                
        %FIRST PART: EXPLOITATION
        if known_map_features(j) == 1
            if new_scores(j).score > 0
                centroid = centroid_vec(j, :);
                r_relative = r(i, 1:3) - centroid;

                for k = 1:length(new_scores(j).completeness)
                    if ~new_scores(j).completeness(k)
                        if strcmp(new_scores(j).type, 'emission')

                            obs_emission(j,i) = acos((dot(r_relative, N(:, j)))/( norm(r_relative) * norm(N(:, j)) )) ;
                            emission_check = obs_emission(j,i) < new_scores(j).angles(k, 2) & obs_emission(j, i) > new_scores(j).angles(k, 1);                            
            
                            %If the emission angle is inside the limit
                            if emission_check
                                r_sun = r_sun_vec(i,:);
                                sun_incidence(j, i) = (dot(r_sun, N(:, j)))/( norm(r_sun) * norm(N(:, j)));
                                incidence_check = sun_incidence(j, i) > 0;
                            
                                %If Incidence angle is inside the limit
                                if incidence_check
                                    R_body_cam = R_body_cam_all(:,:,i);
                                    inside1 = check_FOV(V(F(j, 1), :), r(i, :), fov1, fov2, R_body_cam);
                                    
                                    %If FOV is ok
                                    if inside1
                                        inside2 = check_FOV(V(F(j, 2), :), r(i, :), fov1, fov2, R_body_cam);
            
                                        if inside2
                                            inside3 = check_FOV(V(F(j, 3), :), r(i, :), fov1, fov2, R_body_cam);
            
                                            if inside3
                                                if j == 1067
                                                    s = 2;
                                                end
                                                new_scores(j).completeness(k) = 1;
                                                new_scores(j).actual_range(k) = norm(r_relative);
                                                                                        
                                            end
                                        end
                                    end
                                end
                            end
                        else 
                            if strcmp(new_scores(j).type, 'relative')
                                r_sun = r_sun_vec(i,:);
                                obs_relative(j,i) = acos((dot(r_relative, r_sun))/( norm(r_relative) * norm(r_sun) )) ;
                                relative_check = obs_relative(j,i) > new_scores(j).angles(k, 1) & obs_relative(j,i) < new_scores(j).angles(k, 2);                            
                
                                %If the emission angle is inside the limit
                                if relative_check
                                    obs_emission(j,i) = (dot(r_relative, N(:, j)))/( norm(r_relative) * norm(N(:, j)) ) ;
                                    emission_check = obs_emission(j,i) > 0;                            
                
                                    %If the emission angle is inside the limit
                                    if emission_check
                                        sun_incidence(j, i) = (dot(r_sun, N(:, j)))/( norm(r_sun) * norm(N(:, j)) );
                                        incidence_check = sun_incidence(j, i) > 0;
                                    
                                        %If Incidence angle is inside the limit
                                        if incidence_check
                                            R_body_cam = R_body_cam_all(:,:,i);
                                            inside1 = check_FOV(V(F(j, 1), :), r(i, :), fov1, fov2, R_body_cam);
                                            
                                            %If FOV is ok
                                            if inside1
                                                inside2 = check_FOV(V(F(j, 2), :), r(i, :), fov1, fov2, R_body_cam);
                    
                                                if inside2
                                                    inside3 = check_FOV(V(F(j, 3), :), r(i, :), fov1, fov2, R_body_cam);
                    
                                                    if inside3
                                                        new_scores(j).completeness(k) = 1;
                                                        new_scores(j).actual_range(k) = norm(r_relative);
                                                
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
        %SECOND PART: ZONE DISCOVERY
        else
            if known_map(j) == 0
                r_sun = r_sun_vec(i,:);
                sun_incidence(j, i) = acos((dot(r_sun, N(:, j)))/( norm(r_sun) * norm(N(:, j)) ));
                incidence_check = sun_incidence(j, i) < mapping.incidence(2) & sun_incidence(j, i) > mapping.incidence(1);
                
                
                %If is inside incidence emission and range and in the FOV
                if incidence_check
                    r_relative = r(i, 1:3) - centroid;
                    obs_emission(j,i) = acos((dot(r_relative, N(:, j)))/( norm(r_relative) * norm(N(:, j)) )) ;
                    emission_check = obs_emission(j,i) < mapping.emission(2) & obs_emission(j, i) > mapping.emission(1);
    
                    if emission_check
                        R_body_cam = R_body_cam_all(:,:,i);
                        inside1 = check_FOV(V(F(j, 1), :), r(i, :), fov1, fov2, R_body_cam);
                        
                        if inside1
                            inside2 = check_FOV(V(F(j, 2), :), r(i, :), fov1, fov2, R_body_cam);
    
                            if inside2
                                inside3 = check_FOV(V(F(j, 3), :), r(i, :), fov1, fov2, R_body_cam);
    
                                if inside3
                                    %Set that map zone as known
                                    new_known_map(j) = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %insert scores in a score function
    if i == 1
        mapping_score_t(i) = sum(new_known_map - known_map) * single_score_mapping;
        old_new_known_map = new_known_map;
    else
        mapping_score_t(i) = sum(new_known_map - old_new_known_map) * single_score_mapping;
        old_new_known_map = new_known_map;
    end
    
    %Penalty
    if norm(r(i, :)) < r_impact || norm(r(i, :)) > r_escape
        penalty_of_t(i) = - 5000;
    end
    
    if i == 1
        for j = 1:size(F, 1)
            for k = 1:length(new_scores(j).completeness)
                if new_scores(j).completeness(k) == 1
                    s = ( new_scores(j).completeness(k) - scores(j).completeness(k) )* scores(j).score / length(scores(j).completeness);
                    s = s * exp ( - (new_scores(j).actual_range(k) - new_scores(j).ideal_range )^2 / (2*100) );
                    exploit_score_t(i) = exploit_score_t(i) + s;
                end
            end
        end
        old_scores = new_scores;
    else
        for j = 1:size(F, 1)
            for k = 1:length(new_scores(j).completeness)
                if new_scores(j).completeness(k) == 1
                    s = ( new_scores(j).completeness(k) - old_scores(j).completeness(k) )* old_scores(j).score / length(old_scores(j).completeness);
                    s = s * exp ( - (new_scores(j).actual_range(k) - new_scores(j).ideal_range )^2 / (2*100) );
                    exploit_score_t(i) = exploit_score_t(i) + s;
                end
            end
        end
        old_scores = new_scores;
    end

end

end