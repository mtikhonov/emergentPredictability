function group_assignment = getParetoFrontgroupings(coeffs,coeffs_SE)
nVars = size(coeffs,2);
% kMin = 2; kMax = size(coeffs,1) - 2;
kMax = size(coeffs,1) - 2;
k_arr = [2 3 4:4:kMax];
group_assignment = nan(size(coeffs,1)-1,length(k_arr)*nVars);
iter = 1;
for i = 1:nVars
    for k = k_arr
        % opts = statset('Display','final','MaxIter',5000);
        group_assignment(:,iter) = ukmeans(coeffs(2:end,i),coeffs_SE(2:end,i),k,'Replicates',50,'OnlinePhase','on');
        % group_assignment(:,iter) = kmeans(coeffs(2:end,i),k,'Replicates',500);
        iter = iter + 1;
    end
end