function [TT, XX] = compute_trajectory(x0, action_times, final_times, actions, dyn, options)
n = numel(actions); x = x0(:); TT = cell(1,n); XX = cell(1,n);
for i = 1:n
    a = actions{i};
    if isnumeric(a), dv = a(:);
    elseif isfield(a,'dv'), dv = a.dv(:);
    else, dv = a.delta_v(:);
    end
    x = x + [zeros(3,1); dv];
    t0 = action_times{i}; tf = final_times{i};
    ts = t0:100:tf; if ts(end) ~= tf, ts = [ts tf]; end
    [t_seg, x_seg] = ode78(@(t,y) dyn(t,y), ts, x, options);
    TT{i} = t_seg; XX{i} = x_seg; x = x_seg(end,:)';
end
end
