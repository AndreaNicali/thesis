function [P] = fibonacci_distribution(x_c,y_c,z_c,a,b,c,N)
% Based on:
% http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
P = zeros(floor(N),3); 


% Define Fibonacci Lattice in 2D
F = (1+sqrt(5))/2;
F_2D_fun = @(ii,n) [ii./(F)-floor(ii./(F));ii./n];

% Compute the Fibonacci Lattice up to order N
ii_vect = 0:N;
F_2D = F_2D_fun(ii_vect(1:end-1)+1,N);

% Map the Fibonacci Lattice onto the Ellipsoid

% Definition of Cylindrical Equal Area Projection
alpha_fun = @(x) 2*pi*x;
delta_fun = @(y) acos(1-2*y);

% Computation of right ascension and (inverse) declination (from North
% Pole)
alpha_vect = alpha_fun(F_2D(1,:));
delta_vect = delta_fun(F_2D(2,:));

% Computation of Fibonacci points distribution on the approximating ellipse
Px = x_c+a*sin(delta_vect).*cos(alpha_vect);
Py = y_c+b*sin(delta_vect).*sin(alpha_vect);
Pz = z_c+c*cos(delta_vect);

P(:,:) = [Px(:) Py(:) Pz(:)];

end