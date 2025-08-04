function plotEllipsoidWithKnownRegion(F, V, score_struct, known_map)
% Visualizza un ellissoide con:
% - Facce grigie per zone non conosciute
% - Facce colorate dal giallo (basso score) al rosso (alto score)

    hold on;
    axis equal;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Ellissoide con zone conosciute e punteggiate');
    grid on
    grid minor

    % Estrai punteggi
    scores = arrayfun(@(s) s.score, score_struct);

    % Mappa colore personalizzata
    % Giallo = punteggio 0 -> Rosso = punteggio massimo
    max_score = max(scores);
    if max_score == 0
        max_score = 1; % evita divisione per zero
    end

    face_colors = zeros(size(F,1), 3);

    for i = 1:size(F,1)
        if known_map(i) == 0
            face_colors(i, :) = [0.6 0.6 0.6];  % grigio per sconosciuti
        else
            % Colore da giallo (1,1,0) a rosso (1,0,0)
            s = scores(i) / max_score; % normalizza score
            face_colors(i, :) = [1, 1 - s, 0]; % giallo → rosso
        end
    end

    % Disegna la mesh
    patch('Faces', F, 'Vertices', V, ...
        'FaceVertexCData', face_colors, ...
        'FaceColor', 'flat', ...
        'EdgeColor', 'k'); % Bordi neri

    % Aggiungi legenda colore per le zone conosciute
    colormap(gca, [linspace(1,1,100)' linspace(1,0,100)' zeros(100,1)]); % giallo → rosso
    colorbar;
end

