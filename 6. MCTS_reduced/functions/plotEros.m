function plotEros
% Caricamento file .obj
[obj] = readObj('3dModel\Eros Gaskell 50k poly.obj');

% Estrarre i vertici e le facce
vertices = obj.v;
faces = obj.f.v;

% Plot
f = trisurf(obj.f.v, obj.v(:,1), obj.v(:,2), obj.v(:,3), ...
        'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
camlight; lighting gouraud;
grid on
grid minor
end
