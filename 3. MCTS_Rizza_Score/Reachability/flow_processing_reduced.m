function [t_star,id_t_star,J_star,s] = flow_processing_reduced(tt, J_hist, dJ_dt_hist,data_guidance)



% Filter dJ_dt cutting all small values
dJ_dt_hist(abs(dJ_dt_hist)<1e-13) = 0;

safety_margin_toll = data_guidance.safety_margin;
time_afer_man = data_guidance.DeltaT_after_man; 

% Find time at which you have the last pick
if any(dJ_dt_hist<0)
    % Find the first negative value
    id_first_negative = find(dJ_dt_hist<0,1,'first');
    t_first_negative = tt(id_first_negative);  % first time of constraints violation 
    t_bar = t_first_negative-safety_margin_toll;  
    [~,id_bar] = min(abs(tt-t_bar)); 
    id_bar = id_bar(1); 
    id_t_star = find(dJ_dt_hist(1:id_bar)>0,1,'last');  % last time of positive potential
    if ~isempty(id_t_star)
        % 1) The trajectory is not fully safe and there is some improvement
        % Maneuver after the pick if the safety margin is long enough,
        % otherwise move to another point
    else
        % 2) The trajectory is not fully safe and there is no improvement
        % move to another point
        id_t_star = id_bar; 
    end

    if tt(id_t_star) < tt(1)+time_afer_man
        id_t_star = find(tt>=tt(1)+time_afer_man,1,'first');
    end

    
    s = (tt(id_first_negative)-tt(id_t_star));

else
    % case no safety violation
    t_bar = tt(end)-safety_margin_toll;
    [~,id_bar] = min(abs(tt-t_bar));
    id_bar = id_bar(1);
    id_t_star = find(dJ_dt_hist(1:id_bar)>0,1,'last');
    if ~isempty(id_t_star)
        % 3) the trajectory is always safe on the time horizon and there is
        % some improvement up to a certain extent
        % Maneuvering is performed after the pick has settled if the safety
        % margin is long enough otherwise move to another point
%         if safety_margin < safety_margin_toll
%             % discard solution, not enough safety margin is given
%             id_t_star = []; % !!! Questo è sbagliato, potrebbe finire alla
% %             fine
%         end
    else
        % 4) The trajectory is always safe on the time horizon but does not
        % bring any improvement on the objective function
        % perform the maneuver at the end of the arc
        id_t_star = id_bar; 
    end

    if tt(id_t_star) < tt(1)+time_afer_man
        id_t_star = find(tt>=tt(1)+time_afer_man,1,'first');
    end

    s = (tt(end)-tt(id_t_star));

end
% disp(['Safety margin: ',num2str(s/3600),' h'])


t_star = tt(id_t_star);
J_star = J_hist(id_t_star); 


end