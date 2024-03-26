function [rmse_median,rmse_std,partition_entropy] = mutateParetoFront(X,y,pres_bool,input_format,grp_assign)
% set up format of input data
if strcmp(input_format,'standardized')
    X = (X - mean(X))./std(X);
end
% mutate groupings
mut_grp_assign = [];
for i = 1:size(grp_assign,2)
    mut_grp_assign = [mut_grp_assign mutateassignment(grp_assign(:,i))];
end
% split into training-testing set
seed = 1; numfolds = 100; percent_test = 0.5;
testing_bool = split4TrainingTesting(X,pres_bool,percent_test,numfolds,seed);
% evaluate mutated coarse-grainings
partition_entropy = nan(1,size(mut_grp_assign,2));
rmse_median = nan(size(mut_grp_assign,2),size(y,2));
rmse_std = nan(size(mut_grp_assign,2),size(y,2));
parfor i = 1:size(mut_grp_assign,2)
    [cg_pa,partition_entropy(i)] = coarsegrainabd_v2(X,mut_grp_assign(:,i));
    [~,rmse_median(i,:),rmse_std(i,:)] = linreg(cg_pa,y,testing_bool);
end
end

%% functions
function muts = mutateassignment(groupassign)
grps = [unique(groupassign); max(groupassign)+1];
muts = repmat(groupassign,[1 length(groupassign)*(length(grps)-1)]);
% loop through each strain
iter = 1;
for strain = 1:length(groupassign)
    % loop through all possible alternative assignments
    altgrps = grps(grps~=groupassign(strain));
    for i = 1:length(altgrps)
        muts(strain,iter) = altgrps(i);
        % check for empty group
        for check = 1:max(groupassign)
            if ~ismember(muts(:,iter),check)
                muts(muts(:,iter)>check,iter) = muts(muts(:,iter)>check,iter) - 1;
            end
        end
        iter = iter + 1;
    end
end
% check if effectively same as original
bool = ismember(muts',groupassign','rows');
muts(:,bool) = [];
% check if same as biomass only description
bool = ismember(muts',ones(1,length(groupassign)),'rows');
muts(:,bool) = [];
% check if same as microscopic description
bool = ismember(muts',1:length(groupassign),'rows');
muts(:,bool) = [];
end


