function perturbData = runEnvPerturbation(traj, perturbationName, perturbation)
p = mergestructs(traj.params,perturbation);
p.pertb = cell2mat(arrayfun(@(x)x.b,p.envpool,'un',0));
p.pertK = cell2mat(arrayfun(@(x)x.K,p.envpool,'un',0));
fname = fullfile(traj.params.dataFolder, [perturbationName '_' getHash(p) '.mat']);

T0 = 1e7;
perturbFun = @(params, state) perturbEnvironment(params, state, perturbation.envpool);

if exist(fname,'file') && (~isfield(traj.params,'fromScratch') || ~traj.params.fromScratch)
    fprintf('Loading from file...\n')
    dLoad = load(fname);
    
    if paramsCoincide(traj.params, dLoad.data{1}.params)
        data = dLoad.data;
        data{1}.params.T2h = traj.params.T2h; % function handles don't always load well
        fprintf('Done.\n')
    else
        warning('Hash collision!');
    end
else
    data = recordPerturbationResponse({traj}, T0, perturbFun, perturbationName);
    save(fname,'data');
end
perturbData = data{1}.(perturbationName);
end

%% Helper functions
function hash = getHash(p)
% remove irrelevant/dependent param fields
fields2rmv = {'dataFolder','DEBUG','fromScratch','costStructure','g0','a0',...
              'envpool',};
for ff = 1:length(fields2rmv)
    if isfield(p,fields2rmv{ff})
        p = rmfield(p,fields2rmv{ff});
    end
end
% convert any fields that are structs to arrays
if isstruct(p.environment)
    p.environmentb = p.environment.b;
    p.environmentK = p.environment.K;
    p = rmfield(p,'environment');
end
if isfield(p,'strainpool')
    if isstruct(p.strainpool)
        p.sigmapool = p.strainpool.sigma;
        p.poolabd = p.strainpool.abd;
        p = rmfield(p,'strainpool');
    end
end
% convert fcn handles to strings
p.T2h = func2str(p.T2h);
hash = struct2hash(p);
end

function tf = paramsCoincide(p, pLoaded)
% remove any irrelevant/dependent param fields
fields2rmv = {'dataFolder','DEBUG','fromScratch','costStructure','g0','a0'};
for ff = 1:length(fields2rmv)
    if isfield(p,fields2rmv{ff})
        p = rmfield(p,fields2rmv{ff});
    end
    if isfield(pLoaded,fields2rmv{ff})
        pLoaded = rmfield(pLoaded,fields2rmv{ff});
    end
end
% convert fcn handles to strings for hash conversion
if isa(pLoaded.T2h, 'function_handle')
    pLoaded.T2h = func2str(pLoaded.T2h);
end
if isa(p.T2h, 'function_handle')
    p.T2h = func2str(p.T2h);
end
tf = isequaln(p, pLoaded);
end

function mergedStruct = mergestructs(x,y)
try
    mergedStruct = cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
catch
    % remove duplicate fields from one of the structs
    xfields = fieldnames(x);
    yfields = fieldnames(y);
    fields2rmv = intersect(xfields,yfields);
    for f = 1:length(fields2rmv)
        if isfield(x,fields2rmv{f})
            x = rmfield(x,fields2rmv{f});
        end
    end
    mergedStruct = cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
end
end