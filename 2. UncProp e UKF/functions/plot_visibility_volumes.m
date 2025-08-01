function plot_visibility_volumes(features, F, V, N)
% plot_visibility_volumes - Visualizza i volumi di visualizzabilità delle feature come cupole sferiche
%
% Input:
%   features - struct contenente:
%       .index: indici delle facce
%       .visualization.emission: [min max] in radianti (angolo)
%       .visualization.range: [min max] in km
%       .known: vettore logico che indica se la feature è già conosciuta
%   F - facce (nx3)
%   V - vertici (mx3)
%   N - normali (3xn)

hold on; axis equal;

n_angle = 40;  % risoluzione angolare

for i = 1:length(features.index)
    idx = features.index(i);
    centroid = mean(V(F(idx,:), :), 1);
    normal = N(:, idx); normal = normal / norm(normal);

    theta_max = features.visualization.emission(i,2);
    Rmax = features.visualization.range(i,2);
    Rmin = features.visualization.range(i,1);

    if features.known(i)
        col = 'g';
    else
        col = 'r';
    end

    %% Rotazione verso la normale
    axis0 = [0; 0; 1];
    v = cross(axis0, normal);
    if norm(v) < 1e-6
        Rrot = eye(3);
    else
        c = dot(axis0, normal); s = norm(v);
        vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        Rrot = eye(3) + vx + vx*vx*((1-c)/(s^2));
    end

    %% Cupola a Rmax
    phi = linspace(0, 2*pi, n_angle);       % azimut
    theta = linspace(0, theta_max, n_angle);% zenith
    [PHI, THETA] = meshgrid(phi, theta);
    r = Rmax;

    X = r * sin(THETA) .* cos(PHI);
    Y = r * sin(THETA) .* sin(PHI);
    Z = r * cos(THETA);
    pts = Rrot * [X(:)'; Y(:)'; Z(:)'];
    pts = pts + centroid';
    X = reshape(pts(1,:), size(X));
    Y = reshape(pts(2,:), size(Y));
    Z = reshape(pts(3,:), size(Z));
    surf(X, Y, Z, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    %% Disco di chiusura a Rmin
    r = Rmin;
    Rc = r * sin(theta_max);
    Zc = r * cos(theta_max);
    phiC = linspace(0, 2*pi, n_angle);
    Xc = Rc * cos(phiC);
    Yc = Rc * sin(phiC);
    Zc = Zc * ones(size(phiC));
    pts2 = Rrot * [Xc; Yc; Zc];
    pts2 = pts2 + centroid';
    fill3(pts2(1,:), pts2(2,:), pts2(3,:), col, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

    %% Superficie laterale tra Rmin e Rmax
    phi_side = linspace(0, 2*pi, n_angle);
    [PHI, RADIUS] = meshgrid(phi_side, [Rmin, Rmax]);
    Xs = RADIUS .* sin(theta_max) .* cos(PHI);
    Ys = RADIUS .* sin(theta_max) .* sin(PHI);
    Zs = RADIUS .* cos(theta_max);

    pts_side = Rrot * [Xs(:)'; Ys(:)'; Zs(:)'];
    pts_side = pts_side + centroid';
    Xs = reshape(pts_side(1,:), size(Xs));
    Ys = reshape(pts_side(2,:), size(Ys));
    Zs = reshape(pts_side(3,:), size(Zs));
    surf(Xs, Ys, Zs, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    %% ➕ Testo al centro del disco a Rmin
    text_pos = Rrot * [0; 0; Zc(1)];  % centro del disco
    text_pos = text_pos + centroid';
    text(text_pos(1), text_pos(2), text_pos(3), num2str(idx), ...
         'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

camlight; lighting gouraud;
xlabel('X'); ylabel('Y'); zlabel('Z');

end
