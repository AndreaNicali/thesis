function [measurements] = measurementsFun(y, t, spacecraft_data)

fov1 = spacecraft_data.data_guidance.navigationFov1;
fov2 = spacecraft_data.data_guidance.navigationFov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;
known_map = spacecraft_data.data_asteroids.features.known_map_features;
navigation = spacecraft_data.data_asteroids.navigation_features;

r = y(:, 1:3);
N_steps = size(y,1);
N_faces = size(F,1);

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

maskNavKnown = known_map & ismember((1:N_faces).', navigation(:));
thr = deg2rad(80);

Kmax = N_steps * N_faces;
measurements.val = zeros(2, Kmax);
measurements.feature = zeros(1, Kmax);
measurements.time = zeros(1, Kmax);
measurements.los = zeros(3, Kmax);
iter = 0;

for i = 1:N_steps
    r_i = r(i, :);
    rsu = r_sun_vec(i, :);
    R_body_cam = R_body_cam_all(:, :, i);

    faceInFov = check_FOV(V, F, r_i, fov1, fov2, R_body_cam);

    r_rel = r_i - centroid_vec;
    nr = vecnorm(r_rel, 2, 2);

    EA = acos(max(-1,min(1, dot(r_rel, N', 2) ./ max(nr, eps))));
    IA = acos(max(-1,min(1, (rsu*N)')));

    mask = maskNavKnown & (nr < 50) & (EA < thr) & (IA < thr) & faceInFov;
    idx = find(mask);
    nvis = numel(idx);
    if nvis == 0
        continue
    end

    rr = r_rel(idx, :);
    nn = nr(idx);

    los = (-rr).';
    measurements.los(:, iter+1:iter+nvis) = los;
    measurements.feature(iter+1:iter+nvis) = idx.';
    measurements.time(iter+1:iter+nvis) = t(i);

    z = max(-1,min(1, (-rr(:,3))./max(nn, eps)));
    x = (-rr(:,1))./max(nn, eps);
    yv = (-rr(:,2))./max(nn, eps);
    measurements.val(1, iter+1:iter+nvis) = asin(z).';
    measurements.val(2, iter+1:iter+nvis) = atan2(yv, x).';

    iter = iter + nvis;
end

measurements.val = measurements.val(:, 1:iter);
measurements.feature = measurements.feature(1, 1:iter);
measurements.time = measurements.time(1, 1:iter);
measurements.los = measurements.los(:, 1:iter);
end
