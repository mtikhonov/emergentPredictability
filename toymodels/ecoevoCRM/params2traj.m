function data = params2traj(params, loadOnly)
% Function to interpret and convert parameter structure into inputs for
% runTrajectorSimulation and also acts as load/save wrapper for
% trajectories, checking no parameter collisions occur
if nargin<2
    loadOnly = false;
end
fname = fullfile(params.dataFolder, ['sigmoid_', getHash(params),'.mat']);
data = [];
if exist(fname,'file') && (~isfield(params,'fromScratch') || ~params.fromScratch)
    fprintf('Loading from file...\n')
    dLoad = load(fname);
    
    if paramsCoincide(params, dLoad.data.params)
        data = dLoad.data;
        data.params.T2h = params.T2h; % function handles don't always load well
        fprintf('Done.\n')
    else
        warning('Hash collision! (2 parameter sets are pointing to same filename)')
        fprintf('Continuing on to rerun and override trajectory file...\n')
        pause(5)
    end
end
if isempty(data)
    if loadOnly
        warning('Requested file does not exist: %s', fname);
    else
        data = params2traj_internal(params);
        if  ~isfield(data.params, 'noSave') || ~data.params.noSave
            save(fname,'data');
        end
    end
end
end

function tf = paramsCoincide(p, pLoaded)
% strip irrelevant and dependent param fields
fields2rmv = {'dataFolder','DEBUG','fromScratch','costStructure','g0','a0',...
              'dE_epsilon','envN'};
for ff = 1:length(fields2rmv)
    if isfield(p,fields2rmv{ff})
        p = rmfield(p,fields2rmv{ff});
    end
    if isfield(pLoaded,fields2rmv{ff})
        pLoaded = rmfield(pLoaded,fields2rmv{ff});
    end
end
% convert function handles to string for comparision
if isa(pLoaded.T2h, 'function_handle')
    pLoaded.T2h = func2str(pLoaded.T2h);
end
if isa(p.T2h, 'function_handle')
    p.T2h = func2str(p.T2h);
end
% compare
tf = isequaln(p, pLoaded);
end

function hash = getHash(p)
% remove irrelevant fields
fields2rmv = {'dataFolder','DEBUG','fromScratch','costStructure','g0','a0',...
              'dE_epsilon','envN'};
for ff = 1:length(fields2rmv)
    if isfield(p,fields2rmv{ff})
        p = rmfield(p,fields2rmv{ff});
    end
end
% To get a unique hash string for the parameter set:
% convert function handles to strings
p.T2h = func2str(p.T2h);
% convert structs to arrays
if isstruct(p.environment)
    p.environmentb = p.environment.b;
    p.environmentK = p.environment.K;
    p = rmfield(p,'environment');
end
hash = struct2hash(p);
end

function data = params2traj_internal(params)
% Biochemistry
params.costStructure = getCostStructure(params);
[params.g0, params.a0] = getStartingPoint(params);
[data.trajectory, data.parent] = runTrajectorySimulation(params);
data.params = params;
end

function [g0, a0] = getStartingPoint(params)
rng(params.startSeed,'combRecursive');
% initialize random phenotype (g)
g = rand(100,length(params.costStructure.J))<0.5;
g(:,params.environment.K == 0) = 0; % zero out regions of large cost
% compute its cost
chi = computeCost(g,params.costStructure);
a0 = [];
for i=1:size(g, 1)
    % try possible starting points until the first viable one
    g0 = g(i,:);
    delta = @(a)params.T2h(a, params.environment)*g0'-chi(i);
    if delta(1)*delta(1e12)<0
        a0 = fzero(delta, [1 1e12]);
        break;
    end
end
assert(~isempty(a0),'Starting point not viable for seed %d', params.startSeed);
end

function costStructure = getCostStructure(params)
costStructure.J = params.J0*interactionMatrix_sigmoid(params.Linf, params.nstar, params.delta, params.biochemistrySeed);
costStructure.chii = params.chi0 .* ones(1, params.Linf);
costStructure.c = params.c;
end