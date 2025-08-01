function plotEllipsoidWithKnownRegion(F, V, known_map, score)
% Plot the ellipsoid:
% - Unknown faces in grey
% - Known faces colored by score using jet colormap (0 = blue, >0 = warmer)
% 
% Inputs:
% - F: faces (Nx3)
% - V: vertices (Mx3)
% - known_map: logical vector (Nx1), 1 if known, 0 if unknown
% - score: score per face (Nx1)

% Indici noti e ignoti
known_index = known_map == 1;
unknown_index = known_map == 0;

% Plot facce sconosciute (in grigio)
patch('Faces', F(unknown_index,:), 'Vertices', V, ...
      'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
hold on;

% Se ci sono facce conosciute...
if any(known_index)
    % Estrai i punteggi solo delle facce conosciute
    known_scores = score(known_index);

    % Normalizza i punteggi > 0 per usare la colormap jet
    max_score = max(known_scores);
    min_score = min(known_scores(known_scores > 0));
    % Evita divisione per zero
    if isempty(min_score) || max_score == 0
        min_score = 0; max_score = 1;
    end

    % Mappa dei colori: 256 colori di jet
    cmap = jet(256);

    % Colori finali per ogni faccia conosciuta
    face_colors = zeros(sum(known_index), 3);
    for i = 1:sum(known_index)
        s = known_scores(i);
        if s <= 0
            face_colors(i,:) = [0 0 0.5];  % blu fisso
        else
            % Normalizza tra 1 e 256
            idx = round(1 + (s - min_score) / (max_score - min_score) * 255);
            idx = max(min(idx, 256), 1);  % clampa
            face_colors(i,:) = cmap(idx,:);
        end
    end

    % Plot delle facce conosciute con colorazione personalizzata
    patch('Faces', F(known_index,:), 'Vertices', V, ...
          'FaceVertexCData', face_colors, ...
          'FaceColor', 'flat', ...
          'EdgeColor', 'k', 'LineWidth', 0.3);
end

axis equal; view(3); grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
colorbar; colormap(jet);
title('Ellipsoid with Known Regions Colored by Score');

end
