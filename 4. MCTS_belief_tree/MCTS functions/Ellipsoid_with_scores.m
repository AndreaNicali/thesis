function [F, V, N, C20, C22, A, score_struct] = Ellipsoid_with_scores(a, b, c, n1, n2, centres, vals)

% Griglia angolare
[u, v] = meshgrid(linspace(0, 2*pi, n1), linspace(0, pi, n2));

% Parametrizzazione ellissoide
X = a * cos(u) .* sin(v);
Y = b * sin(u) .* sin(v);
Z = c * cos(v);

% Coefficienti armoniche sferiche
C20 = -5.24618393097E-02;
C22 =  8.23993879858E-02;

% Mesh triangolare
[F, V] = surf2patch(X, Y, Z, 'triangles');

% Normali e aree
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

% Calcolo dei centroidi
centroids = zeros(size(F,1), 3);
for i = 1:size(F,1)
    centroids(i,:) = mean(V(F(i,:),:), 1);
end

% Inizializzazione score e identificativi
scores = zeros(size(F,1), length(centres));
id = zeros(size(F,1), 1);

% Score per ogni feature
sigma = 1;
for s = 1:length(centres)
    c = centroids(centres(s), :);
    dist = vecnorm(centroids - c, 2, 2);
    spot_score = exp(-(dist.^2) / (2*sigma^2)) * vals(s);
    spot_score = spot_score / max(spot_score);
    spot_score(spot_score < 1e-2) = 0;
    
    scores(:, s) = spot_score;
    id(spot_score > 0) = centres(s);
end

% Mappa finale dei punteggi
score_map = sum(scores, 2);

% Costruzione score_struct come vettore di struct
score_struct = repmat(struct('score', 0, 'type', '', 'angles', [], 'completeness', 0), size(F,1), 1);

for f = 1:size(F,1)
    if score_map(f) > 0
        score_struct(f).score = score_map(f);
        feature_idx = id(f);
        
        switch feature_idx
            case {centres(1), centres(3)}
                score_struct(f).type = 'emission';
                score_struct(f).angles = [0 19; 20 39; 40 59; 60 79]*pi/180;
                score_struct(f).completeness = [0; 0; 0; 0];
                score_struct(f).id = feature_idx;
                score_struct(f).ideal_range = 35;
                score_struct(f).actual_range = nan(size(score_struct(f).completeness));

            case {centres(2), centres(4)}
                score_struct(f).type = 'relative';
                score_struct(f).angles = [0 10; 70 80]*pi/180;
                score_struct(f).completeness = [0; 0];   
                score_struct(f).id = feature_idx;
                score_struct(f).ideal_range = 35;
                score_struct(f).actual_range = nan(size(score_struct(f).completeness));

        end
    end
end

end
