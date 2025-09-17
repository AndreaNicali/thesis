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

for i = 1:N_steps
    rs = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    rs = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i)) * (-rs(1:3));
    r_sun_vec(i, :) = rs.';
    R_body_cam_all(:, :, i) = body2camera(y(i,1:3), y(i,4:6));
end
r_sun_vec = r_sun_vec ./ vecnorm(r_sun_vec,2,2);

centroid_vec = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;

is_in_navigation = ismember((1:N_faces).', navigation(:));

time = zeros(N_steps * N_faces, 1);
iter = 0;

thr = deg2rad(80);

for i = 1:N_steps
    r_i = r(i, :);
    rsu = r_sun_vec(i, :);
    R_body_cam = R_body_cam_all(:, :, i);

    r_rel = r_i - centroid_vec;
    nr = vecnorm(r_rel, 2, 2);
    EA = acos( max(-1,min(1, dot(r_rel, N', 2) ./ max(nr, eps))) );
    IA = acos( max(-1,min(1, (rsu*N)')) );

    faceInFov = check_FOV(V, F, r_i, fov1, fov2, R_body_cam);

    mask = known_map & is_in_navigation & (EA < thr) & (IA < thr) & faceInFov;

    nvis = sum(mask);
    if nvis > 0
        time(iter+1:iter+nvis) = t(i);
        iter = iter + nvis;
    end
end

time = time(1:iter);
end
