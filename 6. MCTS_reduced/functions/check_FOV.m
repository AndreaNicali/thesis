function [faceInFov] = check_FOV(V, F, r, fov1, fov2, R_body_cam)

Vcam = (V - r) * R_body_cam.';   % Nx3

z  = Vcam(:,3);
tx = tan(fov1/2);
ty = tan(fov2/2);

inside = (z > 0) & ...
         (abs(Vcam(:,1)) <= tx .* z) & ...
         (abs(Vcam(:,2)) <= ty .* z);

faceInFov = inside(F(:,1)) & inside(F(:,2)) & inside(F(:,3));
end
