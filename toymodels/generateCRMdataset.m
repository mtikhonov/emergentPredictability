function generateCRMdataset(params)
for i = 1:params.invadersN
    missingdata_bool = ~exist(['figureData/CRM_invader_' num2str(i) '.mat'],'file');
end
invadersN = sum(missingdata_bool);
invaders = find(missingdata_bool);
if  ~isempty(invaders)
%% define a single random environment
state.environment.rho = params.rho * ones(1,params.L);
state.environment.K = params.K0 .*(1 + params.dK*randn(1,params.L));
%% Get pool to subsample; put in all strains to find coexisting community
rng(params.seed,'combRecursive');
pool = rand(params.M,params.L)<params.p;
while ~all(any(pool,2))
    nullRows = ~any(pool,2);
    pool(nullRows,:) = rand(sum(nullRows),params.L)<params.p;
end
[~,ord] = sort(sum(pool));
pool = pool(:,ord);
% define cost model
cost = params.c + sum(params.lambda*pool,2) + params.dchi*randn(params.M,1);
full_comm.sigma = pool;
full_comm.chi = cost;
full_comm.abd = params.K0*rand(size(full_comm.chi));
full_comm.environment = state.environment;
full_comm_finalabd = equilibrateEcology(full_comm, params.T2h, params.simTime);
full_comm_finalabd(full_comm_finalabd < 1e-6) = 0;
pool = pool(full_comm_finalabd > 0,:);
cost = cost(full_comm_finalabd > 0);
%% Draw random invader strains
inv_sigma = rand(invadersN,params.L)<params.p;
inpool_bool = ismember(pool,inv_sigma,'rows');
while any(inpool_bool)
    inv_sigma(inpool_bool,:) = rand(sum(inpool_bool),params.L)<params.p;
    inpool_bool = ismember(pool,inv_sigma,'rows');
end
inv_cost = params.c + sum(params.lambda*inv_sigma,2) + params.dchi*randn(invadersN,1);
%% Generate trial community subsets for pathogen experiments
N = size(pool,1);
loop_iter = 1;
for i = invaders
    focalstrain.sigma = inv_sigma(loop_iter,:);
    focalstrain.cost = inv_cost(loop_iter);

    lowdiv_pres_abs_bool = false(params.nTrials, N);
    for trial = 1:params.nTrials
        assert(N>=5,'Full community of strains not diverse enough to make low/mid/high richness subsets (up to 25 strains).')
        % how many strains in subset?
        m = randi([1 5]);
        sel = randsample(N, m, false);
        lowdiv_pres_abs_bool(trial,sel) = true;
    end

    middiv_pres_abs_bool = false(params.nTrials, N);
    for trial = 1:params.nTrials
        assert(N>=15,'Full community of strains not diverse enough to make mid/high richness subsets (up to 25 strains).')
        % how many strains in subset?
        m = randi([11 15]);
        sel = randsample(N, m, false);
        middiv_pres_abs_bool(trial,sel) = true;
    end

    highdiv_pres_abs_bool = false(params.nTrials, N);
    for trial = 1:params.nTrials
        assert(N>=25,'Full community of strains not diverse enough to make high richness subsets (up to 25 strains).')
        % how many strains in subset?
        m = randi([21 25]);
        sel = randsample(N, m, false);
        highdiv_pres_abs_bool(trial,sel) = true;
    end

    %% Run pathogen experiments
    % X_temp = zeros(size(lowdiv_pres_abs_bool));
    X_temp = cell(1,params.nTrials);
    invrate_temp = nan(size(lowdiv_pres_abs_bool,1),1);
    postinvAbd_temp = nan(size(lowdiv_pres_abs_bool,1),1);
    parfor trial = 1:size(lowdiv_pres_abs_bool,1)
        % fprintf('Running trial %d/%d...\n',trial,params.nTrials)
        sel = lowdiv_pres_abs_bool(trial,:);
        state_temp = struct();
        state_temp.sigma = pool(sel,:);
        state_temp.chi = cost(sel);
        state_temp.abd = params.K0*rand(size(state_temp.chi));
        state_temp.environment = state.environment;
        finalabd = equilibrateEcology(state_temp, params.T2h, params.simTime);

        % X_temp(sel,trial) = finalabd;
        X_temp{trial} = finalabd;
        T = finalabd' * state_temp.sigma;
        invrate_temp(trial) = focalstrain.sigma*params.T2h(T,state.environment)'-focalstrain.cost;
        % equilbrate invasion
        state_temp.sigma = [state_temp.sigma; focalstrain.sigma];
        state_temp.chi = [state_temp.chi; focalstrain.cost];
        state_temp.abd = [finalabd; params.K0*rand(1)];
        finalabd = equilibrateEcology(state_temp, params.T2h, params.simTime);
        postinvAbd_temp(trial) = finalabd(end);
    end
    dataCRM.lowdiv.invrate = invrate_temp;
    postinvAbd_temp(postinvAbd_temp<1e-6) = 0;
    dataCRM.lowdiv.postinv_abd = postinvAbd_temp;
    dataCRM.lowdiv.X = zeros(size(lowdiv_pres_abs_bool));
    for trial = 1:params.nTrials
        dataCRM.lowdiv.X(trial,lowdiv_pres_abs_bool(trial,:)) = X_temp{trial};
    end

    X_temp = cell(1,params.nTrials);
    invrate_temp = nan(size(middiv_pres_abs_bool,1),1);
    postinvAbd_temp = nan(size(middiv_pres_abs_bool,1),1);
    parfor trial = 1:size(middiv_pres_abs_bool,1)
        % fprintf('Running trial %d/%d...\n',trial,params.nTrials)
        sel = middiv_pres_abs_bool(trial,:);
        state_temp = struct();
        state_temp.sigma = pool(sel,:);
        state_temp.chi = cost(sel);
        state_temp.abd = params.K0*rand(size(state_temp.chi));
        state_temp.environment = state.environment;
        finalabd = equilibrateEcology(state_temp, params.T2h, params.simTime);

        % X_temp(sel,trial) = finalabd;
        X_temp{trial} = finalabd;
        T = finalabd' * state_temp.sigma;
        invrate_temp(trial) = focalstrain.sigma*params.T2h(T,state.environment)'-focalstrain.cost;
        % equilbrate invasion
        state_temp.sigma(end+1,:) = focalstrain.sigma;
        state_temp.chi(end+1) = focalstrain.cost;
        state_temp.abd = finalabd;
        state_temp.abd(end+1) = params.K0*rand(1);
        finalabd = equilibrateEcology(state_temp, params.T2h, params.simTime);
        postinvAbd_temp(trial) = finalabd(end);
    end
    dataCRM.middiv.invrate = invrate_temp;
    postinvAbd_temp(postinvAbd_temp<1e-6) = 0;
    dataCRM.middiv.postinv_abd = postinvAbd_temp;
    dataCRM.middiv.X = zeros(size(middiv_pres_abs_bool));
    for trial = 1:params.nTrials
        dataCRM.middiv.X(trial,middiv_pres_abs_bool(trial,:)) = X_temp{trial};
    end

    X_temp = cell(1,params.nTrials);
    invrate_temp = nan(size(highdiv_pres_abs_bool,1),1);
    postinvAbd_temp = nan(size(highdiv_pres_abs_bool,1),1);
    parfor trial = 1:size(highdiv_pres_abs_bool,1)
        % fprintf('Running trial %d/%d...\n',trial,params.nTrials)
        sel = highdiv_pres_abs_bool(trial,:);
        state_temp = struct();
        state_temp.sigma = pool(sel,:);
        state_temp.chi = cost(sel);
        state_temp.abd = params.K0*rand(size(state_temp.chi));
        state_temp.environment = state.environment;
        finalabd = equilibrateEcology(state_temp, params.T2h, params.simTime);

        % X_temp(sel,trial) = finalabd;
        X_temp{trial} = finalabd;
        T = finalabd' * state_temp.sigma;
        invrate_temp(trial) = focalstrain.sigma*params.T2h(T,state.environment)'-focalstrain.cost;
        % equilbrate invasion
        state_temp.sigma(end+1,:) = focalstrain.sigma;
        state_temp.chi(end+1) = focalstrain.cost;
        state_temp.abd = finalabd;
        state_temp.abd(end+1) = params.K0*rand(1);
        finalabd = equilibrateEcology(state_temp, params.T2h, params.simTime);
        postinvAbd_temp(trial) = finalabd(end);
    end
    dataCRM.highdiv.invrate = invrate_temp;
    postinvAbd_temp(postinvAbd_temp<1e-6) = 0;
    dataCRM.highdiv.postinv_abd = postinvAbd_temp;
    dataCRM.highdiv.X = zeros(size(highdiv_pres_abs_bool));
    for trial = 1:params.nTrials
        dataCRM.highdiv.X(trial,highdiv_pres_abs_bool(trial,:)) = X_temp{trial};
    end

    %% Clean-up and save
    dataCRM.lowdiv.pres_abs_bool = lowdiv_pres_abs_bool;
    dataCRM.middiv.pres_abs_bool = middiv_pres_abs_bool;
    dataCRM.highdiv.pres_abs_bool = highdiv_pres_abs_bool;
    save(['figureData/CRM_invader_' num2str(i) '.mat'],'dataCRM')
    loop_iter = loop_iter + 1;
end
end