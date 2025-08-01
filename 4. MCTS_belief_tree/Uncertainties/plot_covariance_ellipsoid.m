function plot_covariance_ellipsoid(mu, Sigma, nsigma, color)
    
    % Decompose the covariance matrix
    [V, D] = eig(Sigma);  % V: eigenvectors, D: eigenvalues

    % Generate points on the unit sphere
    [x, y, z] = sphere(50);  % mesh of unit sphere
    xyz = [x(:) y(:) z(:)]'; % 3xN
    
    % Scale by sqrt of eigenvalues and nsigma
    radii = nsigma * sqrt(diag(D));  % scale axes
    ellipsoid_points = V * diag(radii) * xyz;
    
    % Translate by the mean
    ellipsoid_points = ellipsoid_points + mu(:);
    
    % Reshape for surface plotting
    X = reshape(ellipsoid_points(1,:), size(x));
    Y = reshape(ellipsoid_points(2,:), size(y));
    Z = reshape(ellipsoid_points(3,:), size(z));
    
    % Plot
    axis equal; grid on;
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', color);
    plot3(mu(1), mu(2), mu(3), 'r.', 'MarkerSize', 20); % mean point
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title([num2str(nsigma), '\sigma Covariance Ellipsoid']);
 end
