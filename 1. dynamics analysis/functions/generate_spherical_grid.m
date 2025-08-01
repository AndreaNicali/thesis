function grid = generate_spherical_grid(F, V, r_min, r_max, n_r, n_theta, n_phi)
% generate_spherical_grid - Crea una griglia di celle sferiche (r, theta, phi)
% attorno all'ellissoide definito da F, V
%
% Output:
%   grid.centers - [N x 3] coordinate xyz dei centri delle celle
%   grid.edges_r - [n_r+1 x 1] bordi dei bin radiali
%   grid.edges_theta - [n_theta+1 x 1] bordi dei bin zenitali
%   grid.edges_phi - [n_phi+1 x 1] bordi dei bin azimutali

% Trova il centro dell’ellissoide
centroid = mean(V, 1);  % centroide dei vertici

% Genera i bordi dei bin
edges_r = linspace(r_min, r_max, n_r+1);
edges_theta = linspace(0, pi, n_theta+1);       % zenith (0 = nord, pi = sud)
edges_phi = linspace(0, 2*pi, n_phi+1);         % azimuth

% Prealloca centri
centers = [];

% Scorri i volumi
for ir = 1:n_r
    for itheta = 1:n_theta
        for iphi = 1:n_phi
            % Calcola centro della cella in coordinate sferiche
            r_c = (edges_r(ir) + edges_r(ir+1)) / 2;
            theta_c = (edges_theta(itheta) + edges_theta(itheta+1)) / 2;
            phi_c = (edges_phi(iphi) + edges_phi(iphi+1)) / 2;

            % Converti in coordinate cartesiane
            x = r_c * sin(theta_c) * cos(phi_c);
            y = r_c * sin(theta_c) * sin(phi_c);
            z = r_c * cos(theta_c);

            % Trasla rispetto al centroide
            centers(end+1, :) = centroid + [x, y, z];
        end
    end
end

% Output
grid.centers = centers;
grid.edges_r = edges_r;
grid.edges_theta = edges_theta;
grid.edges_phi = edges_phi;
end
