function [measurements] = measurementsFun(y, t, spacecraft_data)

fov1 = spacecraft_data.data_guidance.navigationFov1;
fov2 = spacecraft_data.data_guidance.navigationFov2;

F = spacecraft_data.data_asteroids.Faces;
V = spacecraft_data.data_asteroids.Vertexes;
N = spacecraft_data.data_asteroids.Normals;
known_map = spacecraft_data.data_asteroids.features.known_map_features;

navigation = spacecraft_data.data_asteroids.navigation_features;
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

measurements.val = [];
measurements.feature = [];
measurements.time = [];
measurements.los = [];

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
            if norm_r_rel < 50
            if obs_coaltitude < 80*pi/180
                sun_incidence = acos(dot(r_sun, Nj)/( norm(r_sun) * norm(Nj) ));
                if sun_incidence < 80*pi/180
                    v1 = V(F(j,1), :); v2 = V(F(j,2), :); v3 = V(F(j,3), :);
                    if check_FOV(v1, r_curr, fov1, fov2, R_body_cam) && ...
                       check_FOV(v2, r_curr, fov1, fov2, R_body_cam) && ...
                       check_FOV(v3, r_curr, fov1, fov2, R_body_cam)
                       
                        measurements.los(:, iter) = -r_relative;
                        measurements.feature(iter) = j;
                        measurements.time(iter) = t(i);
                        measurements.val(1, iter) = asin(measurements.los(3, iter)/norm_r_rel);
                        measurements.val(2, iter) = atan2(measurements.los(2, iter)/norm_r_rel, measurements.los(1, iter)/norm_r_rel);
                        iter = iter+1;
                    end
                end
            end
            end
        end
    end
end

end
