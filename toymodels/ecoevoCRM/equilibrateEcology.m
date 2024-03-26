function abdFinal = equilibrateEcology(state, T2h, time)
params = state; params.T2h = T2h; params.time = time;
abdFinal = solveodes(state,T2h,time);
end

function abdFinal = solveodes(state,T2h,time)
% state must have fields: environment, sigma, chi, abd
environment = state.environment;
% ode15s is more accurate but slower.
% If the number of species is above this threshold, run a coarser solve first.
ode15s_threshold = 80;
abdThreshold = 1e-4; % relative; Types above 1 in 10000 for speed
%% COARSE
if sum(state.abd>1e-6) > ode15s_threshold
    a2h = @(a)T2h(a' * state.sigma, environment);
    dadt = @(~,a)a.*(state.sigma * a2h(a)' - state.chi);
    
    options = odeset('Events', @(t,a)stopEventFcn(a, ode15s_threshold, abdThreshold));
    sol.xe = [];
    t0 = 0;
    abd0 = state.abd;
    while isempty(sol.xe) && t0<time
        sol = ode45(dadt, [t0 min(t0+1e4, time)], abd0, options);
        abd0 = sol.y(:,end);
        t0 = sol.x(end);
    end
    abd0 = sol.y(:,end);
    t0 = sol.x(end);
    % Remove phenotypes that had died off, to speed up the fine stage
    alive = abd0/sum(abd0)>abdThreshold;
    abd0 = abd0(alive);
    state.sigma = state.sigma(alive, :);
    state.chi = state.chi(alive);
else
    t0 = 0;
    abd0 = state.abd;
    alive = true(size(abd0));
end

if t0<time
    %% FINE
    a2h = @(a)T2h(a' * state.sigma, environment);
    dadt = @(~,a)a.*(state.sigma * a2h(a)' - state.chi);
    opts_ode15s = odeset('NonNegative',1:length(state.chi));
    sol = ode15s(dadt, [t0, t0+(time-t0)/2, time], abd0, opts_ode15s);
    abd0 = sol.y(:,end);
end
abdFinal = zeros(size(alive));
abdFinal(alive) = abd0;
abdFinal(abdFinal<0) = 0;
end

function [value,isterminal,direction] = stopEventFcn(a,N, abdThresh)
value = N-sum(a/sum(a)>abdThresh); 
isterminal = 1;
direction = 0; 
end
