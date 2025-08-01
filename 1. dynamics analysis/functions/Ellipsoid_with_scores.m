function [F, V, N, C20, C22, A, scores ] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals)

%Create a triangular mesh for the ellipsoid
%F = [3 X N] indicating vertex indexes of each face
%V = [3 X M] indicating vertex coordinates in body fixed frame
%N = [3 X N] Normal vector of each face
%A = [3 X N] Area of each face
%C20, C22 spherical harmonics coefficients

% Griglia angolare
[u, v] = meshgrid(linspace(0, 2*pi, n1), linspace(0, pi, n2));

% Parametrizzazione ellissoide
X = a * cos(u) .* sin(v);
Y = b * sin(u) .* sin(v);
Z = c * cos(v);

%Spherical Harmonic coefficients
C20 = 1/5 * (2*c^2 - a^2 - b^2) / (a^2 + b^2 + c^2);
C22 = 1/20 * (a^2 - b^2) / (a^2 + b^2 + c^2);
 
% Converto la superficie in triangoli
[F, V] = surf2patch(X, Y, Z, 'triangles');

N = zeros(3, size(F, 1)); 
A = zeros(size(F, 1), 1); 

for i = 1:size(F, 1)
    vert1 = F(i, 1);
    vert2 = F(i, 2);
    vert3 = F(i, 3);
    
    v1 = V(vert2, :) - V(vert1, :);
    v2 = V(vert3, :) - V(vert1, :);
    
    normal = cross(v1, v2);
    area = 0.5 * norm(normal);
    A(i) = area;
    
    n = normal / norm(normal); 
    N(:, i) = n;

    % Orientamento normale uscente
    centroide = (V(vert1, :) + V(vert2, :) + V(vert3, :)) / 3;
    if dot(n, centroide) < 0
        N(:, i) = -n;
    end
end

C20 = -5.24618393097E-02;
C22 = 8.23993879858E-02;

% Genera una mappa di punteggio deterministica sulle facce della mesh V, F

% Calcola i centroidi delle facce
centroids = zeros(size(F,1), 3);
for i = 1:size(F,1)
    centroids(i,:) = mean(V(F(i,:),:), 1);
end

scores = zeros(size(F,1), length(centres));

% Per ogni spot, aggiungi una gaussiana localizzata
for s = 1:length(centres)
    c = centroids(centres(s), :);
    dist = vecnorm(centroids - c, 2, 2);
    
    sigma = 1;
    spot_score = exp(-(dist.^2) / (2*sigma^2)) * vals(s);
    
    scores(:, s) = spot_score;
    scores(scores < 10^-2) = 0;

    scores(:, s) = scores(:, s)/max(scores(:, s));

end

scores = sum(scores, 2);

end