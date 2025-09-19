function R = body2camera(r, v)
%Rotation matrix from body frame to camera frame. Camera frame is defined
%as: Z axis in Nadir Pointing, Y axis is the projection of velocity on the
%plane perpendicular to Z axis, X axis complete the frame 

%Nadir pointing
z_cam = -r / norm(r);

%Velocity projection on the plane perpendicular to r
v_proj = v - dot(v, z_cam) * z_cam;  
y_cam = v_proj / norm(v_proj);

%complete the frame
x_cam = cross(y_cam, z_cam);

R = [x_cam; y_cam; z_cam];
end