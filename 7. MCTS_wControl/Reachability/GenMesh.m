function [P0,P] = GenMesh(vv0,drmax,nmax,nr)
% This function generates a uniform sampling inside a unitary sphere
% INPUT
%           vv0: non-dimensional center of sphere [3x1]
%           drmax: radius of the sphere [1x1]
%           nmax: maximum number of points on the outer layer
%           nr: number of layers
% OUTPUT
%           P0: all points
%           P: all points inside the unitary sphere
% Author: Antonio Rizza

% P0 and P should be preallocated (issue opened on this)

% Radial discretization
drr = linspace(0,drmax,nr);

% Compute surface density
% rhoS = nmax/(4*pi*drmax^2);
alpha_factor = nmax/((nr-1)^2); 

% P0 preallocation 
numPointsPerLayers = floor(alpha_factor * ((0:nr-1).^2)); % points per layer
totalPoints = sum(numPointsPerLayers); % total number of samples
P0 = zeros(totalPoints,3); 

% preallocated P0
for k = 1:nr
    dr = drr(k); 
    n = alpha_factor*((k-1)^2);
    [p] = fibonacci_distribution(0,0,0,dr,dr,dr,n);

    if numPointsPerLayers(k) ~= 0 
        P0(sum(numPointsPerLayers(1:k-1))+1:sum(numPointsPerLayers(1:k-1))+numPointsPerLayers(k),:) = p;
    end

    
end

% null control sample
P0 = [zeros(1,3); P0];

% Shift in the center position
P0 = vv0(:)'+P0; 

% Cut points outside the sphere
PointsNorm = vecnorm(P0,2,2); 
[id_ext_bound] = find(PointsNorm<=1.001); 
P = P0(id_ext_bound,:); 

end