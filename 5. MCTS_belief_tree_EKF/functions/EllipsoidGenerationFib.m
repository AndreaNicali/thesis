function [F, V, N, C20, C22, A, score_struct] = EllipsoidGenerationFib(a, b, c, M, centres, vals, known_map)

%------------------------
% PUNTI con reticolo di Fibonacci
%------------------------
gold = (1+sqrt(5))/2;

k = (0:M-1).';
z = 1 - 2*(k+0.5)/M;
theta = 2*pi*k/gold;
r = sqrt(max(0,1 - z.^2));

xs = r .* cos(theta);
ys = r .* sin(theta);
zs = z;

% scala a ellissoide
V = [a*xs, b*ys, c*zs];

% triangolazione della superficie con l’involucro convesso
K = convhull(V);
F = K;

%------------------------
% COEFFICIENTI (come prima)
%------------------------
C20 = -5.24618393097E-02;
C22 =  8.23993879858E-02;

%------------------------
% NORMALI e AREE per faccia
%------------------------
N = zeros(3, size(F, 1));
A = zeros(size(F, 1), 1);

for i = 1:size(F, 1)
    verts = V(F(i,:), :);
    v1 = verts(2,:) - verts(1,:);
    v2 = verts(3,:) - verts(1,:);
    normal = cross(v1, v2);
    area = 0.5 * norm(normal);
    n = normal / norm(normal);
    centroide = mean(verts, 1);
    if dot(n, centroide) < 0
        n = -n;
    end
    N(:, i) = n;
    A(i) = area;
end

%------------------------
% CENTROIDI per faccia
%------------------------
centroids = zeros(size(F,1), 3);
for i = 1:size(F,1)
    centroids(i,:) = mean(V(F(i,:),:), 1);
end

%------------------------
% SCORE features (come prima, con piccola pulizia)
%------------------------
scores = zeros(size(F,1), length(centres));
id = zeros(size(F,1), 1);

sigma = 0.5;
for s = 1:length(centres)
    c = centroids(centres(s), :);
    dist = vecnorm(centroids - c, 2, 2);
    spot = exp(-(dist.^2) / (2*sigma^2)) * vals(s);
    % soglia “hard” (se vuoi mantenere coda morbida, rimuovi queste due righe)
    spot(spot < 1e-2) = 0;
    spot(spot > 0) = vals(s);

    scores(:, s) = spot;
    id(spot > 0) = centres(s);
end

score_map = sum(scores, 2);
total_score = sum(score_map);

%------------------------
% COSTRUZIONE score_struct (come prima)
%------------------------
score_struct = repmat(struct('score', 0, 'type', '', 'angles', [], ...
    'completeness', [0; 0; 0; 0], 'id', [], 'ideal_range', [], 'actual_range', []), size(F,1), 1);

for f = 1:size(F,1)
    if score_map(f) > 0
        score_struct(f).score = score_map(f)/total_score;
        feature_idx = id(f);
        switch feature_idx
            case {centres(1)}
                score_struct(f).type = 'emission';
                score_struct(f).angles = [0 19; 20 39; 40 59; 60 79]*pi/180;
                score_struct(f).completeness = [0; 0; 0; 0];
                score_struct(f).id = feature_idx;
                score_struct(f).ideal_range = 37;
                score_struct(f).actual_range = nan(size(score_struct(f).completeness));
            case {centres(2)}
                score_struct(f).type = 'relative';
                score_struct(f).angles = [0 19; 20 39; 40 59; 60 79]*pi/180;
                score_struct(f).completeness = [0; 0; 0; 0];
                score_struct(f).id = feature_idx;
                score_struct(f).ideal_range = 37;
                score_struct(f).actual_range = nan(size(score_struct(f).completeness));
        end
    end
end


end
