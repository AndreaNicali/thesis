function [J_of_t, dJdt, features_post, known_map] = score_rizza(y, t, spacecraft_data)

%SCORE FUNCTION FOR RIZZA'S INSPIRED PROBLEM

%Compute score as time passed inside a visibility region of a feature.
% INPUT:
% y: [NX6] trajectory of the s/c
% t: [NX1] propagation time in ephemeris time
% spacecraft data: struct containing a set of required values

% OUTPUT:
% J_of_t: [NX1] score gained at time up to time t
% dJdt: [NX1] numerically computed derivative of score with time
% feature_post: struct containing the updated state of feature visited time
% known_map: map containing 1 if a face is known, 0 if unknown (for this
%            problem, all the asteroid is unknown.

%Extract data from spacecraft struct
r_impact = spacecraft_data.data_guidance.r_impact;
r_escape = spacecraft_data.data_guidance.r_escape;
fov1 = spacecraft_data.data_guidance.fov1;
fov2 = spacecraft_data.data_guidance.fov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;
features_pre = spacecraft_data.data_asteroids.features;

known_map = spacecraft_data.data_asteroids.known_map;
mapping = spacecraft_data.data_asteroids.mapping;

%Initialize variables
r = y(:, 1:3);

obs_coaltitude = zeros(size(F, 1), size(y, 1));
t_obs = zeros(size(F, 1), size(y, 1));
obs_distance = zeros(size(F, 1), size(y, 1));
sun_incidence = zeros(size(F, 1), size(y, 1));
features_post = features_pre;
feature_discovered = zeros(size(features_post.known));

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
        
        %If we are on a feature
        if any(features_post.index == j)
           centroid = ( V(F(j, 1), :) + V(F(j, 2), :) +  V(F(j, 3), :) )/3;
           r_relative = r(i, 1:3) - centroid;

           ind = find(features_post.index == j);
           obs_coaltitude(j,i) = acos((dot(r_relative, N(:, j)))/( norm(r_relative) * norm(N(:, j)) )) ;
           emission_check = obs_coaltitude(j,i) < features_post.visualization.emission(ind, 2) & obs_coaltitude(j, i) > features_post.visualization.emission(ind, 1);
            
           %If the emission angle is inside the limit
            if emission_check
                r_sun = r_sun_vec(i,:);
                sun_incidence(j, i) = acos((dot(r_sun, N(:, j)))/( norm(r_sun) * norm(N(:, j)) ));
                incidence_check = sun_incidence(j, i) < features_post.visualization.incidence(ind, 2) & sun_incidence(j, i) > features_post.visualization.incidence(ind, 1);
                
  
                %If relative angle between incidence and emission is
                %inside the limit
                if incidence_check
                    obs_distance(j, i) = norm(r_relative);
                    range_check = obs_distance(j, i) < features_post.visualization.range(ind, 2) & obs_distance(j, i) > features_post.visualization.range(ind, 1);
                    %If range is inside the limit
                    if range_check
                        %Compute observation time and
                        %update available score and
                        %completeness level
                        if i == size(y, 1)
                            t_obs(j, i) = ( t(i) - t(i-1));
                        else
                            t_obs(j, i) = ( t(i+1) - t(i));
                        end
                
                        if features_post.known(ind)

                            if features_post.availablescore(ind) > 0
                                features_post.score(ind) = features_post.score(ind) + t_obs(j, i);

                                if features_post.score(ind) > features_post.totalscore(ind)
                                    features_post.score(ind) = features_post.totalscore(ind);
                                    features_post.availablescore(ind) = 0;
                                    features_post.completeness(ind) = 1;
                                else
                                    features_post.availablescore(ind) = features_post.totalscore(ind) - features_post.score(ind);
                                    features_post.completeness(ind) = features_post.score(ind)/features_post.totalscore(ind);
                                end
                            end
                        else
                            feature_discovered(ind) = 1;
                        end
                    end
                end
            end
        
        %If that map zone is uknown
        elseif known_map(j) == 0
            r_sun = r_sun_vec(i,:);
            sun_incidence(j, i) = acos((dot(r_sun, N(:, j)))/( norm(r_sun) * norm(N(:, j)) ));
            incidence_check = sun_incidence(j, i) < mapping.incidence(2) & sun_incidence(j, i) > mapping.incidence(1);
            
            %If is inside incidence emission and range and in the FOV
            if incidence_check
                r_relative = r(i, 1:3) - centroid;
                obs_coaltitude(j,i) = acos((dot(r_relative, N(:, j)))/( norm(r_relative) * norm(N(:, j)) )) ;
                emission_check = obs_coaltitude(j,i) < mapping.emission(2) & obs_coaltitude(j, i) > mapping.emission(1);

                if emission_check
                obs_distance(j, i) = norm(r_relative);
                range_check = obs_distance(j, i) < mapping.range(2) & obs_distance(j, i) > mapping.range(1);

                    if range_check
                        R_body_cam = R_body_cam_all(:,:,i);
                        inside1 = check_FOV(V(F(j, 1), :), r(i, :), fov1, fov2, R_body_cam);
                        
                        if inside1
                            inside2 = check_FOV(V(F(j, 2), :), r(i, :), fov1, fov2, R_body_cam);

                            if inside2
                                inside3 = check_FOV(V(F(j, 3), :), r(i, :), fov1, fov2, R_body_cam);

                                if inside3
                                    %Set that map zone as known
                                    known_map(j) = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %insert scores in a score function
    J_features(i) = sum( (features_post.completeness - features_pre.completeness) );

    J_mapping(i) = sum(known_map - spacecraft_data.data_asteroids.known_map);
    
    J_of_t(i) = J_features(i) + J_mapping(i);
    
    %Penalty
    if norm(r(i, :)) < r_impact || norm(r(i, :)) > r_escape
        J_of_t(i) = J_of_t(i) - 5000;
    end

end

for i = 1:length(features_post.known)
    if feature_discovered(i)
        features_post.known(i) = 1;
    end
end

%Compute numerical derivative
dJdt = diff(J_of_t) ./ diff(t);

dJdt = [dJdt, dJdt(end)];

end