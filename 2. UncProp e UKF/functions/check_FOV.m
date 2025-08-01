function [inside] = check_FOV(V, r, fov1, fov2, R_body_cam)

%Checks if a particular vertex is inside the rectangular fov of the camera

%Compute distance of FOV (length of the projection of s/c - vertex vector on the nadir pointing vector)
proj1 = dot(V, r)/dot( r,  r) * r;
dist1 = norm(r - proj1);

%Compute boundaries in camera frame (rectangular shape)
x_max1 = tan(fov1/2)*dist1;
x_min1 = -x_max1;
y_max1 = tan(fov2/2)*dist1;
y_min1 = -y_max1;

%Rotate vertex in camera frame
vert1 = R_body_cam*(V-r)';

%Check if inside the limit
inside_x1 = vert1(1) > x_min1 && vert1(1) < x_max1;
inside_y1 = vert1(2) > y_min1 && vert1(2) < y_max1;

inside = inside_x1 && inside_y1;
end