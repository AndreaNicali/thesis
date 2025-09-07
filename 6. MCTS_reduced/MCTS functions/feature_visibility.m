function [time] = feature_visibility(y, t, spacecraft_data)

fov1 = spacecraft_data.data_guidance.navigationFov1;
fov2 = spacecraft_data.data_guidance.navigationFov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;
known_map = spacecraft_data.data_asteroids.features.known_map_features;
navigation = spacecraft_data.data_asteroids.navigation_features;

r = y(:, 1:3);
N_steps = size(y, 1);
N_faces = size(F, 1);

r_sun_vec = zeros(N_steps, 3);
R_body_cam_all = zeros(3, 3, N_steps);

% Precompute rotation matrices and Sun vectors
for i = 1:N_steps
    r_sun = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    r_sun = -r_sun(1:3);
    r_sun = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i)) * r_sun;
    r_sun_vec(i, :) = r_sun';
    R_body_cam_all(:, :, i) = body2camera(y(i,1:3), y(i,4:6));
end

% Precompute centroids
centroid_vec = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;

% Prepare logical index for features in navigation
is_in_navigation = ismember(1:N_faces, navigation);

% Allocate output time array
time = zeros(N_steps * N_faces, 1); % preallocate upper bound
iter = 0;

% Main loop
for i = 1:N_steps
    r_i = r(i, :);
    r_sun = r_sun_vec(i, :);
    R_body_cam = R_body_cam_all(:, :, i);

    for j = 1:N_faces
        if ~known_map(j) || ~is_in_navigation(j), continue, end

        r_rel = r_i - centroid_vec(j, :);
        n_j = N(:, j);
        
        % Emission angle
        dot_rn = dot(r_rel, n_j);
        norm_rn = norm(r_rel) * norm(n_j);
        emission_angle = acos(dot_rn / norm_rn);
        if emission_angle >= deg2rad(80), continue, end

        % Incidence angle
        dot_sn = dot(r_sun, n_j);
        norm_sn = norm(r_sun) * norm(n_j);
        incidence_angle = acos(dot_sn / norm_sn);
        if incidence_angle >= deg2rad(80), continue, end

        % FOV check
        verts = V(F(j, :), :);
        visible = all(arrayfun(@(v) check_FOV(verts(v,:), r_i, fov1, fov2, R_body_cam), 1:3));

        if visible
            iter = iter + 1;
            time(iter) = t(i);
        end
    end
end

% Trim to actual size
time = time(1:iter);

end
