function [trajectory, parent] = runTrajectorySimulation(params)
% Runs an eco-evolutionary trajectory from a specified starting point
% (timepoint, environment condition, phenotype(s), and abundance(s)), using
% a hybrid stochastic-deterministic Gillespie algorithm to numerical solve
% the dynamics.
% See 'computeData' for parameter descriptions.

g0 = params.g0;
a0 = params.a0;
env = params.environment;
rng(params.evoSeed,'combRecursive');
if isempty(g0)
    t0 = a0.time;
else
    t0 = 0;
end
% if simTime is a vector, make sure that these specific timepoints are all included
timepoints = t0 + params.simTime;
epochT = diff([0 params.simTime]);
trajectory = [];
for ee = 1:length(epochT)
    if isempty(trajectory)
        [trajectory, parent] = runTrajectorySimulation_1epoch(epochT(ee), g0, a0, params, env);
    else
        par = parent;
        [tt, parent] = runTrajectorySimulation_1epoch(epochT(ee), [], trajectory(end), params, env);
        if abs(tt(1).time - timepoints(ee-1))>1e-5
            warning('Ignoring the fact that simulation stopped early in epoch %d.', ee)
            tt(1).time = timepoints(ee);
        end
        trajectory = [trajectory(1:end-1), tt];
        parent(1:length(par)) = par;
    end
end
end

function [trajectory, parent] = runTrajectorySimulation_1epoch(T, g0, a0, params, env)
% Runs a trajectory from the supplied starting point, for a specified
% length of time or until equilibration
if isfield(params, 'DEBUG')
    DEBUG = params.DEBUG;
else
    DEBUG = false;
end

trajectoryLengthLimit = params.trajectoryLengthLimit;
if ~isa(env,'function_handle')
    env = @(t)env;
end

% Initialize the state
currentState = initialState(g0, a0, env, params);
time0 = currentState.time; % normally 0, unless we're continuing a past trajectory
timeFinal = time0+T;

% Initialize the trajectory
trajectory(1) = currentState;
if isinf(params.trajectoryLengthLimit)
    trajectory(2) = currentState; % only two states to populate: initial and final
else
    trajectory(trajectoryLengthLimit) = currentState; % preallocate memory
end

% to track species descent:
parent = zeros(1,1e4,'uint32');
% all initial species are automatically labeled as having derived from
% "genotype 0" (i.e. were present initially)

nextT = 2; % next trajectory entry to be set
while currentState.time < timeFinal && nextT<=params.trajectoryLengthLimit
    if isempty(currentState.abd)
        error('Inconsistency. Simulation step too large?');
    end
    
    % Propagate ecology (solve CRM odes) until event occurs (stopR)
    %   Possible events:
    %       -- mutation is selected
    %       -- record state (stored in params.tickEvery)
    %   Note: rates for these events change during this deterministic step,
    %   but the their total is tracked until it reaches the random number
    %   'stopR', at which point of the random events takes place...
    stopR = exprnd(1);
    currentState = propagateEcologicalDynamics(currentState, timeFinal, stopR, env, params);
    
    % Either we're done, or it's time to apply an event
    if currentState.time < timeFinal
        [luckyMutation, df] = pickEvent(currentState, params);
        [currentState, parent] = applyMutation(luckyMutation, df, currentState, parent, params);
        if luckyMutation>0
            currentState = updateDependentFields(currentState, env, params.T2h);
        end
        [currentState, anyRemoved] = pruneLosers(currentState, params.minAbd);
        if anyRemoved
            currentState = updateDependentFields(currentState, env, params.T2h);
        end        
    end
    
    % Record the state
    trajectory(nextT) = currentState;
    if ~isinf(params.trajectoryLengthLimit)
        % If params.trajectoryLengthLimit==Inf, keep nextT at 2 and never
        % increase it.This will keep writing the current state into the
        % same slot each time.
        nextT = nextT+1;
    end
    
    if DEBUG
        printDebugMessage(currentState, params);
    end
end
if ~isinf(params.trajectoryLengthLimit)
    trajectory(nextT:end) = [];
end
parent(currentState.nextId:end)=[];
end

function v = dRdt(a,state, params)
    % total event rate dynamics, accounting for both record event and
    % mutation events
    state.abd = a;
    v = rTotal(state, params);
end

function [value,isterminal,direction] = stopEventFcn(x,stopR)
% standard event function to stop ode solver. In this case, occurs when
% integrated total rate R = stopR
% Returns the time at which the event occurs
value = x(1)-stopR;
isterminal = 1;
direction = 0; 
end

function state = propagateEcologicalDynamics(state, timeFinal, stopR, environment, params)
% integrate ODE for variable: x = [R, a] where R is integrated total rate
% of any event occurring:
%          R = integral dt r_tot(t) from state.t to state.t+tau,
%          where r_tot(t) = rate of mutations occuring + rate of recording
%          state and tau = the time at which R = stopR
T2h = params.T2h; % fcn handle converting total demand/exploitation to niche availability
a2h = @(t,a)T2h(a' * state.sigma, environment(t));
dxdt = @(t,x)[dRdt(x(2:end),state, params); % dynamics of total event rate
              x(2:end).*(state.sigma * a2h(t,x(2:end))' - state.chi)]; % abd dynamics

% event function that stops ode solver if R = stopR
options = odeset('Events', @(t,x)stopEventFcn(x,stopR));

if length(state.abd)>params.ode15s_threshold %params.L*3
    sol = ode45(dxdt, [state.time, timeFinal], [0; state.abd], options);
else
    % more accurate, but slow if too many types coexist
    sol = ode15s(dxdt, [state.time, timeFinal], [0; state.abd], options);
end
state.abd = sol.y(2:end,end);
state.abd(state.abd<0) = 0; % can happen within numerical error, ~1e-10  
state.time = sol.x(end);

state = updateDependentFields(state, environment, T2h);
end

function state = updateDependentFields(state, env, T2h)
    state.environment = env(state.time);
    totalDemand = state.abd' * state.sigma;
    state.h = T2h(totalDemand, state.environment);
    state.Delta = state.sigma*state.h' - state.chi;
end

function [state, parent] = applyMutation(luckyMutation, df, state, parent, params)
% Apply the event
if luckyMutation==0
    state.whoMutated = NaN;
    state.whereMutated = NaN;
    state.df = NaN;
    state.newTypeBorn = false;
    return
else
    % Add the new species next to its parent
    
    % Construct the new species
    [K,L] = size(state.sigma);
    [parentK, mutL] = ind2sub([K,L],luckyMutation);
    indexOfParent = state.idx(parentK);
    state.whoMutated = indexOfParent;
    state.whereMutated = mutL;
    state.df = df;
    
    mutG = state.sigma(parentK,:);
    mutG(mutL) = ~mutG(mutL); % flip bit mutL
    % is this mutant already in the population?
    same = all(state.sigma==mutG,2);
    if any(same)
        state.newTypeBorn = false;
        mutK = find(same); % used to transfer parent abd to extant mutant abd
    else
        state.newTypeBorn = true;
        % add a new species
        indexOfMutant = state.nextId;
        parent(indexOfMutant) = indexOfParent; 
        state.nextId = state.nextId+1;

        % free up the space next to parent, shifting everyone AFTER parent
        shift = (parentK+1):K;
        state.idx(shift+1,1) = state.idx(shift,1);
        state.abd(shift+1,1) = state.abd(shift,1);
        state.chi(shift+1,1) = state.chi(shift,1);
        state.sigma(shift+1,:) = state.sigma(shift,:);
        state.mutChi(shift+1,:) = state.mutChi(shift,:);
        
        % Initialize the mutant
        mutK = parentK+1;
        state.idx(mutK,1) = indexOfMutant;
        state.abd(mutK,1) = 0;
        state.sigma(mutK,:) = mutG;
        state.chi(mutK,1) = computeCost(mutG, params.costStructure);
        mutChi = computeMuts(mutG, params.costStructure);
        if any(mutChi<0)
            % a legacy check: with the new cost model, shouldn't be possible
            warning('my:PropCom:neg','Forbidding a negative-cost mutant.')
        end
        mutChi(mutChi<=0)=NaN;
        state.mutChi(mutK,:) = mutChi;
    end
    % Move 1/df individuals from parent to mutant
    abdToMove = min(state.abd(parentK), max(1/df,1));
    state.abd(mutK) = state.abd(mutK) + abdToMove;
    state.abd(parentK) = state.abd(parentK) - abdToMove;
end
end

function [state, anyRemoved] = pruneLosers(state, minAbd)
losers = state.abd<1 | (state.abd<minAbd & state.Delta<0);
anyRemoved = any(losers);
if anyRemoved
    winners = ~losers;
    state.sigma = state.sigma(winners,:);
    state.chi = state.chi(winners);
    state.mutChi = state.mutChi(winners,:);
    state.abd = state.abd(winners);
    state.idx = state.idx(winners);
end
end

function [state, costStructure] = initialState(g0, a0, env, params)
costStructure = params.costStructure;
if isempty(g0)
    state = a0;
    state.environment = env(state.time); % 1 by L, double
    assert(all(size(state.environment.K)==size(state.h)),'L mismatch');
    if isempty(state.mutChi)
        [K,L] = size(state.sigma);
        state.mutChi = NaN(K,L);
        for k=1:K
            state.mutChi(k,:) = computeMuts(state.sigma(k,:), costStructure); % K-by-L, double
        end
    end
else
    assert(size(g0,1)==length(a0),'Initial abundance and initial genotypes: size mismatch.')
    state.time = 0;
    state.sigma = logical(g0); % K by L
    [K,L] = size(g0);
    state.chi = computeCost(g0, costStructure); % K-by-1, double 
    state.mutChi = NaN(K,L);
    for k=1:K
        state.mutChi(k,:) = computeMuts(g0(k,:), costStructure); % K-by-L, double
    end
    state.abd = a0(:); % K-by-1, double
    state.environment = env(state.time); % 1 by L, double
    state.idx = uint32(1:K)';    % K by 1, uint32 (id of each extant genotype, to track phylogeny)
    state.nextId = uint32(K+1); % uint32 (next unassigned id)
    state.df = NaN;
    state.whoMutated = NaN;
    state.whereMutated = NaN;
    state.newTypeBorn = false;
    state.h = NaN(size(state.environment));
end
state = updateDependentFields(state,env,params.T2h);
end

function [r, fMut] = mutationEventRates(state, mutationRate)
    % For each species, L possible mutants
    % So K-by-L matrix of mutants. Their chi is precomputed. 
    % Their harvest differs from the harvest of their parents by the
    % corresponding h: 
    %   If mutant flips position i from 0 to 1 -> harvest is +h_i
    %   If mutant flips position i from 1 to 0 -> harvest is -h_i
    % So the sign is the opposite from sigma of parent -0.5
    mutantHarvestChange = state.h .* sign(0.5-state.sigma);
    fMut = state.Delta + state.chi + mutantHarvestChange - state.mutChi;
    % compare to mean fitness in population
%     fMean = (state.Delta' *state.abd)/sum(state.abd);
%     fMut = fMut-fMean;
    fMut(fMut<0 | isnan(fMut)) = 0;
    r = mutationRate * fMut .* state.abd;
end

function r = rTotal(state, params)
    % Sum of all rates of mutations occurring & "record" event rate
    recordEventRate = 1/params.tickEvery;
    mus = mutationEventRates(state, params.mutationRate);
    r = sum(mus(:)) + recordEventRate;
end

% returns a number (index of the lucky mutation to escape drift),
% or zero if it's just the record-event bell (not a mutation)
function [luckyEvent, df] = pickEvent(currentState, params)
    [mus, dfs] = mutationEventRates(currentState, params.mutationRate);
    mus = [1/params.tickEvery; mus(:)];
    dfs = [NaN; dfs(:)];
    luckyEvent = randsample(length(mus), 1, true, mus)-1;
    df = dfs(luckyEvent+1);
end

function printDebugMessage(currentState, params)
    totAbd = sum(currentState.abd);
    edges = 10.^[-Inf, -7:0];
    phenotypeBins = discretize(currentState.abd, totAbd*edges);
    abdHist = accumarray(phenotypeBins, 1);
    mus = mutationEventRates(currentState, params.mutationRate);
    mutPotential = sum(mus,2);
    mutPotHist = accumarray(phenotypeBins, mutPotential);
    str = [...
        sprintf('Time = %.0f, total abundance %g\n', currentState.time, totAbd),...
        sprintf('Phenotypes:  %3d alive (%3d growing).', length(currentState.Delta), sum(currentState.Delta>0)),...
        sprintf('\nAbd bin:'), sprintf('%5.2g     ', edges)               ,...
        sprintf('\nSpecies:     '), sprintf('%5d     ', abdHist)          ,...
        sprintf('\nMut pot:     '), sprintf('%5.2g     ', mutPotHist)];
    fprintf('\n\n%s\n',str);
end

