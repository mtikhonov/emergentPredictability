function permutedPareto = permutationTest(X,y,pres_bool,input_format,complexity_opt,obj2perm)
rng(1,'combRecursive')
if strcmp(obj2perm,'permute_y')
    y = y(randperm(length(y)));
elseif strcmp(obj2perm,'permute_Xclark')
    for species = 1:size(X,2)
        pres_inds = find(pres_bool(:,species));
        row_perm = randsample(pres_inds,length(pres_inds));
        X(pres_inds,species) = X(row_perm,species);
    end
elseif strcmp(obj2perm,'permute_Xkehe')
    for species = 1:size(X,2)
        X(:,species) = X(randperm(size(X,1)),species);
    end
end
% get coarse-graining Pareto Front of permuted dataset
permutedPareto = getCGParetoFront(X,y,pres_bool,input_format,[],complexity_opt);