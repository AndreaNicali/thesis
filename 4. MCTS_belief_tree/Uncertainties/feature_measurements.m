function [measurements] = feature_measurements(y, t, spacecraft_data)

fov1 = spacecraft_data.data_guidance.fov1;
fov2 = spacecraft_data.data_guidance.fov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;
known_map = spacecraft_data.data_asteroids.features.known_map_features;

navigation = spacecraft_data.data_asteroids.navigation;
r = y(:, 1:3);

r_sun_vec = zeros(size(y,1), 3);
R_body_cam_all = zeros(3, 3, size(y,1));

for i = 1:size(y, 1)
    r_sun = cspice_spkgeo(2000433, t(i), 'ECLIPJ2000', 10);
    eclip2eros = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t(i));
    r_sun_vec(i, :) = eclip2eros*(-r_sun(1:3));
    R_body_cam_all(:, :, i) = body2camera(y(i,1:3), y(i,4:6));
end

iter = 1;

measurements.coaltitude.val = [];
measurements.coaltitude.feature = [];
measurements.coaltitude.time = [];

measurements.azimuth.val = [];
measurements.azimuth.feature = [];
measurements.azimuth.time = [];

measurements.range.val = [];
measurements.range.feature = [];
measurements.range.time = [];

for i = 1:size(y, 1)
    R_body_cam = R_body_cam_all(:,:,i);
    r_curr = r(i, :);
    r_sun = r_sun_vec(i,:);
    for j = 1:size(F, 1)
        if known_map(j) && any(navigation == j)
            centroid = ( V(F(j, 1), :) + V(F(j, 2), :) +  V(F(j, 3), :) )/3;
            r_relative = r_curr - centroid;
            norm_r_rel = norm(r_relative);
            Nj = N(:, j);
            obs_coaltitude = acos(dot(r_relative, Nj)/( norm_r_rel * norm(Nj) ));
            if obs_coaltitude < 80*pi/180
                sun_incidence = acos(dot(r_sun, Nj)/( norm(r_sun) * norm(Nj) ));
                if sun_incidence < 80*pi/180
                    v1 = V(F(j,1), :); v2 = V(F(j,2), :); v3 = V(F(j,3), :);
                    if check_FOV(v1, r_curr, fov1, fov2, R_body_cam) && ...
                       check_FOV(v2, r_curr, fov1, fov2, R_body_cam) && ...
                       check_FOV(v3, r_curr, fov1, fov2, R_body_cam)
                        Z_topo = Nj;
                        ref = [0;0;1];
                        if abs(dot(Z_topo, ref)) > 0.99, ref = [1;0;0]; end
                        X_topo = cross(ref, Z_topo); X_topo = X_topo/norm(X_topo);
                        Y_topo = cross(Z_topo, X_topo);
                        R_body2topo = [X_topo'; Y_topo'; Z_topo'];
                        r_topo = R_body2topo*r_relative';
                        measurements.coaltitude.val(iter) = acos(r_topo(3)/norm(r_topo));
                        measurements.coaltitude.feature(iter) = j;
                        measurements.coaltitude.time(iter) = t(i);
                        measurements.azimuth.val(iter) = atan2(r_topo(1), r_topo(2));
                        measurements.azimuth.feature(iter) = j;
                        measurements.azimuth.time(iter) = t(i);
                        measurements.range.val(iter) = norm_r_rel;
                        measurements.range.feature(iter) = j;
                        measurements.range.time(iter) = t(i);
                        iter = iter + 1;
                    end
                end
            end
        end
    end
end

end
