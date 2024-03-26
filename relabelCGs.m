function dataOut = relabelCGs(X,y,pres_bool,input_format,grp_assign,complexity_opt)
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
desc_complexity = nan(1,size(mut_grp_assign,2));
FVUmedian = nan(size(mut_grp_assign,2),1); FVUlower = nan(size(mut_grp_assign,2),1); FVUupper = nan(size(mut_grp_assign,2),1);
AVUmedian = nan(size(mut_grp_assign,2),1); AVUlower = nan(size(mut_grp_assign,2),1); AVUupper = nan(size(mut_grp_assign,2),1);
ypredmed = nan(length(y),size(mut_grp_assign,2)); ypredstd = nan(length(y),size(mut_grp_assign,2));
parfor i = 1:size(mut_grp_assign,2)
    [cg_pa,partition_entropy] = coarsegrainabd(X,mut_grp_assign(:,i));
    if strcmp(complexity_opt,'num_classes')
        desc_complexity(i) = max(mut_grp_assign(:,i));
    else
        desc_complexity(i) = partition_entropy;
    end
    tempregdata = linreg_training_testing(cg_pa,y,testing_bool);
    FVUmedian(i) = tempregdata.FVUmed; FVUlower(i) = tempregdata.FVUlower; FVUupper(i) = tempregdata.FVUupper;
    AVUmedian(i) = tempregdata.AVUmed; AVUlower(i) = tempregdata.AVUlower; AVUupper(i) = tempregdata.AVUupper;
    ypredmed(:,i) = tempregdata.ypredmed; ypredstd(:,i) = tempregdata.ypredstd;
end
dataOut.group_assignment = mut_grp_assign;
dataOut.FVUmedian = FVUmedian; dataOut.FVUlower = FVUlower; dataOut.FVUupper = FVUupper;
dataOut.AVUmedian = AVUmedian; dataOut.AVUlower = AVUlower; dataOut.AVUupper = AVUupper;
dataOut.ypredmedian = ypredmed; dataOut.ypredstd = ypredstd;
dataOut.desc_complexity = desc_complexity;
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


