function [J,T,S] = EvaluateSet(t0,xx0,P0, U0,spacecraft_data)
% This function performs the evaluation of a given RS
% INPUT
%       -t0: time [s]
%       -xx0: Cartesian state [6x1] [km km/s]
%       -U0: set of nodes to evaluate
% OUTPUT
%       -J: score function at the evalution nodes [s]
%       -T: next maneuver time at the evalutaion nodes [s]
%       -S: safety margin a the evaluation nodes [s]
% References: Rizza PhD thesis, Equations 4.8, 4.12, 4.13
% Author: Antonio Rizza


data_guidance = spacecraft_data.data_guidance;


% Initialize score
J = zeros(size(U0,2),1);
T = zeros(size(U0,2),1);
S = zeros(size(U0,2),1);

% Propagate samples and evaluate scores
tf = t0+data_guidance.Th_max;
tstep = 100;
uu0_1 = U0(:,1);
xx0_1 = xx0+[zeros(3,1);uu0_1];
[xx,tt] = integrate_ode_reachability(xx0_1,t0,tf,tstep);

% Score function
spacecraft_data.data_guidance.r_impact = spacecraft_data.data_guidance.r_impact+2;
spacecraft_data.data_guidance.r_escape = spacecraft_data.data_guidance.r_escape-10;

[J_hist, dJ_dt_hist] = total_score(xx, tt, P0, spacecraft_data); % Can be any user defined function

% Process the flow to compute the next maneuvering point
% Note that flow_processing uses the derivatives of the score
% over time (a different logic must be defined if this variable
% is not availe, or dJ/dt can be computed numerically)

[t_star,~,J_star,s] = flow_processing_reduced(tt, J_hist,dJ_dt_hist,data_guidance);

T(1) = t_star; % Maneuvering time
J(1) = J_star; % Score at the next maneuvring epoch
S(1) = s; % Safety margin

for k = 2:size(U0,2)
    uu0_k = U0(:,k);
    xx0_k = xx0+[zeros(3,1);uu0_k];
    [xx_k,tt_k] = integrate_ode_reachability(xx0_k,t0,tf,tstep); % Propagate sample

    % Score function
    [J_hist, dJ_dt_hist] = total_score(xx_k, tt_k,P0, spacecraft_data); % Can be any user defined function

    % Process the flow to compute the next maneuvering point
    % Note that flow_processing uses the derivatives of the score
    % over time (a different logic must be defined if this variable
    % is not availe, or dJ/dt can be computed numerically)

    [t_star,~,J_star,s] = flow_processing_reduced(tt_k, J_hist, dJ_dt_hist, data_guidance);


    % Assign score
    T(k) = t_star; % Maneuvering time
    J(k) = J_star; % Score at the next maneuvring epoch
    S(k) = s; % Safety margin

end

end