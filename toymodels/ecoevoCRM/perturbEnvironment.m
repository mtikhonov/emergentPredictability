function  params = perturbEnvironment(params0, state, envList)
% if envList is a set of environments, use those as supplied.
% if envList is just a number, sequentially perturb each resource by that amount 
if numel(envList)==1
    epsilon = envList;
    envN = length(state.environment);
    envList = state.environment + epsilon*eye(envN);
else
    envN = size(envList, 1);
end

params0.g0 = [];
params0.a0 = state;
params = repmat({params0}, [1, envN]);
for i=1:envN
    params{i}.environment = envList(i,:);
    params{i}.evoSeed = i;
end
end