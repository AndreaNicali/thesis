function plot_ellipsoid_with_scores(F, V, scores)
% Plotta un ellissoide colorato in base ai punteggi per faccia
%
% Input:
%   F      - [n x 3] indici dei vertici delle facce
%   V      - [m x 3] coordinate dei vertici
%   scores - [n x 1] punteggi associati a ogni faccia

    % Crea il plot con trisurf
    figure;
    trisurf(F, V(:,1), V(:,2), V(:,3), scores, ...
            'EdgeColor', 'none', 'FaceAlpha', 1);
        
    colormap jet;
    axis equal;
    view(3);
    lighting gouraud;
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
