function [J_of_t, dJdt, new_features, new_known_map] = score_features(y, t, spacecraft_data)

%Extract data from spacecraft struct
r_impact = spacecraft_data.data_guidance.r_impact;
r_escape = spacecraft_data.data_guidance.r_escape;
fov1 = spacecraft_data.data_guidance.fov1;
fov2 = spacecraft_data.data_guidance.fov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;

known_map = spacecraft_data.data_asteroids.mapping.known_map;
mapping = spacecraft_data.data_asteroids.mapping;
new_known_map = known_map;

score_map = spacecraft_data.data_asteroids.features.score_map;
n_observations = spacecraft_data.data_asteroids.features.n_observations;
max_observations = spacecraft_data.data_asteroids.features.max_observations;
features = spacecraft_data.data_asteroids.features;
new_features = features;

observations_run = zeros(length(t), length(max_observations));
observations_score = zeros(length(t), length(max_observations));
mapping_score_t = zeros(size(t));
penalty_of_t = zeros(size(t));
% SET UP RELATIVE SCORE
initial_score_mapping = 100;
n_zones_known = sum(known_map);
n_zones_unknown = length(known_map) - n_zones_known;

abs_score_mapping = initial_score_mapping - initial_score_mapping/(n_zones_unknown + n_zones_known) * n_zones_known;

abs_score_features = sum(score_map(known_map == 1 & n_observations < max_observations));

abs_score_total = abs_score_mapping + abs_score_features;

single_score_mapping = abs_score_mapping/n_zones_unknown;

%Initialize variables
r = y(:, 1:3);

obs_emission = zeros(size(F, 1), size(y, 1));

sun_incidence = zeros(size(F, 1), size(y, 1));

J_of_t = zeros(size(t));
J_features = zeros(size(t));
J_mapping = zeros(size(t));

r_sun_vec = zeros(size(y,1), 3);
R_body_cam_all = zeros(3, 3, size(y,1));

%Extract rotation matrices before
for i = 1:size(y, 1)
    r_sun = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    r_sun = r_sun(1:3);
    eclip2eros = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i));
    r_sun = eclip2eros*(-r_sun);
    r_sun_vec(i, :) = r_sun;
    R_body_cam_all(:, :, i) = body2camera(y(i,1:3), y(i,4:6));
end


for i = 1:size(y, 1)

    for j = 1:size(F, 1)
        
        %Compute topocentric frame
        centroid = ( V(F(j, 1), :) + V(F(j, 2), :) +  V(F(j, 3), :) )/3;
        Z_topo = N(:, j);
        X_topo = [-centroid(2); centroid(1); 0];
        X_topo = X_topo/norm(X_topo);
        Y_topo = cross(Z_topo, X_topo);
        R_body2topo = [X_topo'; Y_topo'; Z_topo'];
        
        %FIRST PART: ZONE OBSERVATION
        if known_map(j) == 1
            if score_map(j) > 0
                r_relative = r(i, 1:3) - centroid;

                obs_emission(j,i) = acos((dot(r_relative, N(:, j)))/( norm(r_relative) * norm(N(:, j)) )) ;
                emission_check = obs_emission(j,i) < features.emission(2) & obs_emission(j, i) > features.emission(1);
            
                %If the emission angle is inside the limit
                if emission_check
                    r_sun = r_sun_vec(i,:);
                    sun_incidence(j, i) = acos((dot(r_sun, N(:, j)))/( norm(r_sun) * norm(N(:, j)) ));
                    incidence_check = sun_incidence(j, i) < features.incidence(2) & sun_incidence(j, i) > features.incidence(1);
                
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

                                    if n_observations(j) < max_observations(j)
                                        observations_run(i, j) = 1;
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        
        %SECOND PART: ZONE DISCOVERY
        elseif known_map(j) == 0

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
    
    %insert scores in a score function
    if i == 1
        mapping_score_t(i) = sum(new_known_map - known_map) * single_score_mapping;
        old_new_known_map = new_known_map;
    else
        mapping_score_t(i) = sum(new_known_map - old_new_known_map) * single_score_mapping;
        old_new_known_map = new_known_map;
    end

    % if i == 1
    %     J_of_t(i) =  mapping_score(i) + observations_score(i);
    % else
    %     J_of_t(i) = J_of_t(i-1) + mapping_score(i) + observations_score(i);
    % end

    %Penalty
    if norm(r(i, :)) < r_impact || norm(r(i, :)) > r_escape
        penalty_of_t(i) = - 5000;
    end

end

for j = 1:size(observations_run, 2)
    for i = 1:size(observations_run, 1)
        if observations_run(i, j) == 1
            observations_run(:, j ) = 0;
            observations_run(i, j) = 1;
            observations_score(i, j) = observations_run(i, j)*score_map(j);
            break
        end
    end
end

observations_score_t = sum(observations_score, 2);

n_observations = n_observations + sum(observations_run, 1)';

J_per_t = observations_score_t' + mapping_score_t + penalty_of_t;
J_of_t = cumsum(J_per_t);

%Compute numerical derivative
dJdt = diff(J_of_t) ./ diff(t);

dJdt = [dJdt, dJdt(end)];

new_features.n_observations = n_observations;
end