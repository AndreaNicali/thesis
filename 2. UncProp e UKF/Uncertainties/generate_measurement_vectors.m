function [origin_vectors, direction_vectors] = generate_measurement_vectors(measurements, r, t, spacecraft_data)

    F = spacecraft_data.data_asteroids.Faces;
    V = spacecraft_data.data_asteroids.Vertexes;

    num_meas = length(measurements.coaltitude.feature);
    origin_vectors = zeros(num_meas, 3);
    direction_vectors = zeros(num_meas, 3);

    for i = 1:num_meas
        feature_id = measurements.coaltitude.feature(i);

        [~, idx_time] = min(abs(t - measurements.coaltitude.time(i)));
        spacecraft_pos = r(idx_time, :);

        triangle = F(feature_id, :);
        centroid = mean(V(triangle, :), 1);

        vec = centroid - spacecraft_pos;

        origin_vectors(i, :) = spacecraft_pos;
        direction_vectors(i, :) = vec / norm(vec);
    end
end
