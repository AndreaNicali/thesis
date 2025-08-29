function [uu_opt,th_opt,J_opt,U,J,T,S,I] = exploreU(xx0,P0,t0,spacecraft_data)
% This function computes the reachable set
% INPUT:
%       - xx0: Cartesian state at the epoch of the first maneuver [6x1] [km km/s]
%       - tt0: Epoch of the first maneuver [1x1]
%       - spacecraft_data: spacecraft_data class
% OUTPUT
%       - uu0_opt: optimal DV [3x1] [km/s]
%       - th_opt: epoch of the next maneuver [1x1] [ET]
% Author: Antonio Rizza


% Extract data guidance
% spacecraft_data = spacecraft.spacecraft_data;
data_guidance = spacecraft_data.data_guidance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid search parameters (SHOULD BE DEFINED IN data_guidance_class)

% number of most promising samples from initial mesh
ns = 3;
% number of refinement iterations
n_it = 1;

% initial sampling parameters
nr0 = 4;
nmax0 = 20; %QUESTO E' CAMBIATO
drmax0 = 0.9;
vv0_0 = [0 0 0];
alpha_factor0 = nmax0/((nr0-1)^2);
% P0 preallocation
numPointsPerLayers0 = floor(alpha_factor0 * ((0:nr0-1).^2)); % points per layer
totalPoints0 = sum(numPointsPerLayers0)+1; % total number of samples

% refined sampling parameters
vv0_refined = [0 0 0];
drmaxRefined = 0.99;
nmaxRefined = 5;
nrRefined = 3;
alpha_factorRefined = nmaxRefined/((nrRefined-1)^2);
% P0 preallocation
numPointsPerLayersRefined = floor(alpha_factorRefined * ((0:nrRefined-1).^2)); % points per layer
totalPointsRefined = sum(numPointsPerLayersRefined); % total number of samples

% total samples considered during refinement
totalPointsRefinedIter = 1+totalPointsRefined;
totalPointsRefined = (totalPointsRefinedIter)*ns*n_it;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output
uu_opt = zeros(3,1);
th_opt = 0;

% preallocation
U = zeros(3,totalPoints0+totalPointsRefined);
J = zeros(totalPoints0+totalPointsRefined,1);
T = zeros(totalPoints0+totalPointsRefined,1);
S = zeros(totalPoints0+totalPointsRefined,1);
I = zeros(ns*n_it,1);

if data_guidance.ReachabilityScoreComputation == 1
    % Evaluate J afterwork

     if data_guidance.ReachabilityExplorationScheme == 2
        % Define cube size and number of points
        L = 1;   % Cube side length
        n = 25;  % Number of grid points along each axis

        % Generate a uniform grid within the cube
        x = linspace(-L, L, n);
        y = linspace(-L, L, n);
        z = linspace(-L, L, n);
        [X, Y, Z] = ndgrid(x, y, z);


        % Flatten the grid into a set of (x, y, z) coordinates
        P = [X(:), Y(:), Z(:)];
        P(vecnorm(P',2,1)>1,:) = [];


        np = size(P,1);
        % Assume P and np to be given !

        % Initialize variables
        I = zeros(np,1);
        J = zeros(np,1);
        T = zeros(np,1);
        S = zeros(np,1);

        % Pick P0
        %     P0 =
        %     random_vector = randn(50,1);

        %     P0 = P_full(random_vector>0,:);
        nmax0 = 20; 
        ns = 3; 
        [~,P0_c] = GenMesh(vv0_0,drmax0,nmax0,nr0); % ... first on the unitary sphere ...
        %     P0 = P0_c;
        for k = 1:size(P0_c,1)
            [~,id_k] = min(vecnorm((P0_c(k,:)-P)',2,1));
            id_0(k) = id_k;
        end

        % assume id_0 to be given !
        I(id_0) = 1;


%         figure()
%         %     plot3(P(:,1), P(:,2), P(:,3), '.b');
%         plot3(P(id_0,1), P(id_0,2), P(id_0,3), 'or');
%         hold on
% 
%         axis equal
%         grid minor
%         box on
%         drawnow

        % Evaluate the initial mesh
        [J0,T0,S0] = EvaluateSet(t0,xx0,P0,data_guidance.DeltaV_max*P(id_0,:)',spacecraft_data);

        % initial sampling
        %     U(id_0,:) = U0;
        J(id_0,:) = J0;  % objective function
        T(id_0,:) = T0;
        S(id_0,:) = S0;  % safety margin

        [~,id_best] = maxk(J,ns);
        counter = 0;
        counter_2 = 0;
        stop_flag = 0;
        while stop_flag == 0
            % Loop over the best nodes
            for l = 1:ns
                Pl = P(id_best(l),:);
                [~,id_l] = mink(vecnorm((P-Pl)'),3); % closest nodes
                if any(I(id_l)) == 1
                    % remove already explored nodes
                    id_l(I(id_l)==1) = [];
                end
                if isempty(id_l)
                    % All the neighbourhood of node id_best(l) has been
                    % explored
                    counter = counter+1;
                else
                    % Evaluate the initial mesh
                    [J(id_l),T(id_l),S(id_l)] = EvaluateSet(t0,xx0,P0,data_guidance.DeltaV_max*P(id_l,:)', spacecraft_data);

                    % Track nodes as explored
                    I(id_l) = 1;
% 
% %                     Plot for debug
%                     plot3(P(id_l,1), P(id_l,2), P(id_l,3), 'or');
%                     drawnow
                end
            end
            if counter == ns % exit condition
                stop_flag = 1;
                U = P'*data_guidance.DeltaV_max;
                [J_opt,id_max] = max(J);
                uu_opt = U(:,id_max);
                th_opt = T(id_max);
            else
                counter = 0;
                [~,id_best] = maxk(J,ns);
            end

            counter_2 = counter_2 + 0;
            if counter_2 == 100
                stop_flag = 1;
                U = P'*data_guidance.DeltaV_max;
                [J_opt,id_max] = max(J);
                uu_opt = U(:,id_max);
                th_opt = T(id_max);
            end

        end

    end




elseif data_guidance.ReachabilityScoreComputation == 2




end

end