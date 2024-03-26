function db = recordPerturbationResponse(db, T0, perturbFun, fName)
% Apply perturbations as specified by perturbFcn to each element of db,
% which is a cell array of evo trajectory structures.
% Record the perturbed trajectories as a new field in db cell, 'fName'
% perturbFcn is a fcn handle accepting arguments (params, state) and
% returns cell array of param structs, each with a different perturbation
% instance
NN = numel(db);
for r=1:NN
    d = db{r};
    state = d.trajectory(end);
    params = perturbFun(d.params, state);

    iterN = length(params);
    newData = cell(1,iterN);

    parfor i=1:iterN
        fprintf('%d/%d\n',i,iterN);
        p = params{i};
        p.simTime = T0;
        newData{i} = runTrajectorySimulation(p);
    end    
    db{r}.(fName) = newData;
end
end