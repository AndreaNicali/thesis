function dy = dynamicsEllipsoidEKF( t, y, mass, omega, C20, C22, tao)

%Position and velocity
r = y(1:3);
v = y(4:6);
etaR = y(7:9);
PHI = reshape(y(10:end), 9,9);
% Distance from the primary

%Compute model spherical armonic perturbation and solar radiation pressure
[acc_shp, G] = rhs_modelSTM(y, mass, C20, C22);

G(1,1) = G(1,1)+omega(2)^2+omega(3)^2;
G(1,2) = G(1,2)-omega(1)*omega(2);
G(1,3) = G(1,3)-omega(1)*omega(3);
G(2,1) = G(1,2);
G(2,2) = G(2,2)+omega(1)^2+omega(3)^2;
G(2,3) = G(2,3)-omega(2)*omega(3);
G(3,1) = G(1,3);
G(3,2) = G(2,3);
G(3,3) = G(3,3) + omega(1)^2+omega(2)^2;

dGdv = zeros(3,3);
dGdv(1,2) = 2*omega(3);
dGdv(1,3) = -2*omega(2);
dGdv(2,1) = -2*omega(3);
dGdv(2,3) = 2*omega(1);
dGdv(3,1) = 2*omega(2);
dGdv(3,2) = -2*omega(1);

A = [zeros(3, 3), eye(3), zeros(3);
     G, dGdv , eye(3);
     zeros(3), zeros(3), -1/tao*eye(3)];

Phidot = A*PHI;

% Set the derivatives of the state
dy = [ v;
 - 2*cross(omega, v) - cross(omega, cross(omega, r)) + acc_shp;
 -etaR/tao;
 Phidot(:)];

%omegaR = mvnrnd(zeros(3, 1), QR(7:9,7:9));

dy = dy + [zeros(3,1);
            etaR;
            zeros(3,1);
            zeros(81, 1)];

end
