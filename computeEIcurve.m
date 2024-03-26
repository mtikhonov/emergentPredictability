function dataOut = computeEIcurve(X,y,pa_table,input_format,manual_info,complexity_opt)
%% Specify input format (eg., standardized)
if strcmp(input_format,'standardized')
    X = (X - mean(X))./std(X);
end
nStrains = size(X,2);
%% Fit microscopic descriptions
[coeffs,coeffs_SE] = getmircoscopicfit(X,y);
%% Get Pareto Front groupings
group_assignment = getParetoFrontgroupings(coeffs,coeffs_SE);
group_assignment = [ones(nStrains,1) group_assignment (1:nStrains)'];
% add manual groupings based on inspection
if ~isempty(manual_info)
    if manual_info.varID == 1 % butyrate
        if strcmp(manual_info.div,'low')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'AC'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'BV'),end) = 3;
        elseif strcmp(manual_info.div,'mid')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,{'AC','EL'}),end) = 2;
            group_assignment(ismember(manual_info.varnames,'RI'),end) = 3;
        elseif strcmp(manual_info.div,'high')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'AC'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'DP'),end) = 3;
        end
        
    elseif manual_info.varID == 2 % acetate
        if strcmp(manual_info.div,'low')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'AC'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'BP'),end) = 3;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'BP'),end) = 2;

        elseif strcmp(manual_info.div,'mid')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,{'AC','EL'}),end) = 2;
            group_assignment(ismember(manual_info.varnames,{'BA','BP'}),end) = 3;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'EL'),end) = 2;

        elseif strcmp(manual_info.div,'high')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'ER'),end) = 2;
        end

    elseif manual_info.varID == 3 % lactate
        if strcmp(manual_info.div,'low')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'AC'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'CC'),end) = 3;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'CC'),end) = 2;

        elseif strcmp(manual_info.div,'mid')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,{'AC','EL'}),end) = 2;
            group_assignment(ismember(manual_info.varnames,{'RI','BP'}),end) = 3;

        elseif strcmp(manual_info.div,'high')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'BA'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'CA'),end) = 3;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'BA'),end) = 2;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'CA'),end) = 2;
        end

    elseif manual_info.varID == 4 % succinate
        if strcmp(manual_info.div,'low')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'PC'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'CC'),end) = 3;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'PC'),end) = 2;

        elseif strcmp(manual_info.div,'mid')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,{'DL','BP','BL','BA','DP','ER'}),end) = 2;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'DL'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'BO'),end) = 3;

        elseif strcmp(manual_info.div,'high')
            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,{'BF','RI','CA','BT','BL','PC','BY','ER'}),end) = 2;

            group_assignment(:,end+1) = ones(1,25);
            group_assignment(ismember(manual_info.varnames,'ER'),end) = 2;
            group_assignment(ismember(manual_info.varnames,'BA'),end) = 3;
        end
    end
end
%% Split into training/testing
seed = 1; numfolds = 100; percent_test = 0.5;
testing_bool = split4TrainingTesting(X,pa_table,percent_test,numfolds,seed);
%% Evaluate coarse-graining
desc_complexity = nan(1,size(group_assignment,2));
FVUmedian = nan(size(group_assignment,2),1); FVUlower = nan(size(group_assignment,2),1); FVUupper = nan(size(group_assignment,2),1);
AVUmedian = nan(size(group_assignment,2),1); AVUlower = nan(size(group_assignment,2),1); AVUupper = nan(size(group_assignment,2),1);
ypredmed = nan(length(y),size(group_assignment,2)); ypredstd = nan(length(y),size(group_assignment,2));
parfor i = 1:size(group_assignment,2)
    [cg_pa,partition_entropy] = coarsegrainabd(X,group_assignment(:,i));
    if strcmp(complexity_opt,'num_classes')
        desc_complexity(i) = max(group_assignment(:,i));
    else
        desc_complexity(i) = partition_entropy;
    end
    tempregdata = linreg_training_testing(cg_pa,y,testing_bool);
    FVUmedian(i) = tempregdata.FVUmed; FVUlower(i) = tempregdata.FVUlower; FVUupper(i) = tempregdata.FVUupper;
    AVUmedian(i) = tempregdata.AVUmed; AVUlower(i) = tempregdata.AVUlower; AVUupper(i) = tempregdata.AVUupper;
    ypredmed(:,i) = tempregdata.ypredmed; ypredstd(:,i) = tempregdata.ypredstd;
end
dataOut.group_assignment = group_assignment;
dataOut.FVUmedian = FVUmedian; dataOut.FVUlower = FVUlower; dataOut.FVUupper = FVUupper;
dataOut.AVUmedian = AVUmedian; dataOut.AVUlower = AVUlower; dataOut.AVUupper = AVUupper;
dataOut.ypredmedian = ypredmed; dataOut.ypredstd = ypredstd;
dataOut.desc_complexity = desc_complexity;