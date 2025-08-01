function plot_spherical_grid_cells(grid, center)
% plot_spherical_grid_cells - Plotta i "cubotti sferici" della griglia
%
% Input:
%   grid   - struttura con edges_r, edges_theta, edges_phi
%   center - [1x3] centro globale della griglia

hold on;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Cubotti sferici della griglia');

n_r = length(grid.edges_r) - 1;
n_theta = length(grid.edges_theta) - 1;
n_phi = length(grid.edges_phi) - 1;

res = 5;  % risoluzione triangolazione

for ir = 1:n_r
    r0 = grid.edges_r(ir);
    r1 = grid.edges_r(ir+1);
    for itheta = 1:n_theta
        th0 = grid.edges_theta(itheta);
        th1 = grid.edges_theta(itheta+1);
        for iphi = 1:n_phi
            ph0 = grid.edges_phi(iphi);
            ph1 = grid.edges_phi(iphi+1);

            % superfici radiali (interno ed esterno)
            [PHI, THETA] = meshgrid(linspace(ph0, ph1, res), linspace(th0, th1, res));
            for r = [r0, r1]
                X = r * sin(THETA) .* cos(PHI);
                Y = r * sin(THETA) .* sin(PHI);
                Z = r * cos(THETA);
                surf(X + center(1), Y + center(2), Z + center(3), ...
                     'FaceAlpha', 0.15, 'EdgeColor', 'k', 'FaceColor', 'cyan');
            end

            % superfici laterali in φ (piani azimutali)
            for phi = [ph0, ph1]
                [R, THETA] = meshgrid(linspace(r0, r1, res), linspace(th0, th1, res));
                X = R .* sin(THETA) .* cos(phi);
                Y = R .* sin(THETA) .* sin(phi);
                Z = R .* cos(THETA);
                surf(X + center(1), Y + center(2), Z + center(3), ...
                     'FaceAlpha', 0.15, 'EdgeColor', 'k', 'FaceColor', 'yellow');
            end

            % superfici laterali in θ (piani zenitali), evitiamo i poli
            for theta_val = [th0, th1]
                if abs(theta_val) < 1e-3 || abs(theta_val - pi) < 1e-3
                    continue;  % salta i poli per evitare artefatti
                end
                [R, PHI] = meshgrid(linspace(r0, r1, res), linspace(ph0, ph1, res));
                X = R .* sin(theta_val) .* cos(PHI);
                Y = R .* sin(theta_val) .* sin(PHI);
                Z = R .* cos(theta_val) * ones(size(PHI));
                surf(X + center(1), Y + center(2), Z + center(3), ...
                     'FaceAlpha', 0.15, 'EdgeColor', 'k', 'FaceColor', 'magenta');
            end
        end
    end
end

camlight; lighting gouraud;
end
