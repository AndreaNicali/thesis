function scores = generate_patch_scores(V, F, centres)
% Genera una mappa di punteggio deterministica sulle facce della mesh V, F

% Calcola i centroidi delle facce
centroids = zeros(size(F,1), 3);
for i = 1:size(F,1)
    centroids(i,:) = mean(V(F(i,:),:), 1);
end

scores = zeros(size(F,1), 1);

% Per ogni spot, aggiungi una gaussiana localizzata
for s = 1:length(centres)
    c = centroids(centres(s), :);
    dist = vecnorm(centroids - c, 2, 2);
    
    sigma = 1;
    spot_score = exp(-(dist.^2) / (2*sigma^2));
    
    scores = scores + spot_score;
end

% Normalizza
scores = scores / max(scores);
scores(scores < 10^-2) = 0;
end
