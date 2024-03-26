function testing_bool = split4TrainingTesting(X,pa_table,percent_test,numfolds,seed)
rng(seed,'combRecursive')

[~,~,ic] = unique(pa_table,'rows');
testing_bool = nan(size(X,1),numfolds);
temp = rand(numfolds, max(ic));
[~, order] = sort(temp, 2);
randsamples = order(:,1:floor(max(ic)*percent_test));
for k = 1:numfolds
    testing_bool(:,k) = ismember(ic,randsamples(k,:));
end
testing_bool = logical(testing_bool);