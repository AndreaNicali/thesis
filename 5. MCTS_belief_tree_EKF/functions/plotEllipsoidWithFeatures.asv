function plotEllipsoidWithFeatures(F, V, known_map, nav_feat)
% plotEllipsoidWithFeatures - Plotta l'ellissoide con regioni conosciute e le feature
%
% Input:
%   F - facce (Nx3)
%   V - vertici (Mx3)
%   known_map - vettore logico (Nx1), 1 se faccia nota, 0 altrimenti
%   features - struct con:
%       .index: indici delle facce feature
%       .known: vettore logico come sopra

% Indici facce note e ignote
known_index = known_map == 1;
unknown_index = known_map == 0;

% --- Plot superfici ---
hold on;
patch('Faces', F(unknown_index,:), 'Vertices', V, ...
      'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');  % grigio
patch('Faces', F(known_index,:), 'Vertices', V, ...
      'FaceColor', [0 0.45 1], 'EdgeColor', 'k', 'LineWidth', 0.3);  % blu

% --- Plot centri delle feature ---
for i = 1:length(nav_feat)
    idx = nav_feat(i);
    face_verts = V(F(idx,:), :);
    centroid = mean(face_verts, 1);

    if known_map(i)
        plot3(centroid(1), centroid(2), centroid(3), ...
              'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    end
end
view(3); grid on;
grid minor;
axis equal;

end
