function subsampPareto = downsampleTest(X,y,pres_bool1,pres_bool2,input_format,complexity_opt)
% downsample low-diversity dataset to match high-diversity dataset size
totalobs = sum(unique(pres_bool2,'rows'),'all'); % pres_bool2 assumed to be inputted as the one to match
[pa_table_uniq,~,pa_ic] = unique(pres_bool1,'rows','stable');

nsamps = 5;
rng(1,'combRecursive');
subsampPareto = cell(1,nsamps);
for samp = 1:nsamps
    clearvars comm_ind_samp y_sample X_sample pa_sample
    rowperminds = randperm(size(pa_table_uniq,1));
    pa_table_perm = pa_table_uniq(rowperminds,:);
    sumobs = 0; rowcntr = 1;
    while sumobs < totalobs
        comm_ind_samp(rowcntr) = rowperminds(rowcntr);
        sumobs = sumobs + sum(pa_table_perm(rowcntr,:));
        rowcntr = rowcntr + 1;
    end
    ic_outer = pa_ic(ismember(pa_ic,comm_ind_samp','rows'));
    y_sample = y(ic_outer); % only 1 variable should be supplied
    X_sample = X(ic_outer,:);
    pa_sample = pres_bool1(ic_outer,:);
    % get coarse-graining Pareto Front for each subsample dataset
    subsampPareto{samp} = getCGParetoFront(X_sample,y_sample,pa_sample,input_format,[],complexity_opt);
end
