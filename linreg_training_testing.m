function dataOut = linreg_training_testing(X,y,leaveOut_bool)
% train and test
for i = size(leaveOut_bool,2):-1:1
    % remove replicates
    ytest = y(leaveOut_bool(:,i),:);
    Xtest = X(leaveOut_bool(:,i),:);
    ytrain = y; ytrain(leaveOut_bool(:,i)) = [];
    Xtrain = X; Xtrain(leaveOut_bool(:,i),:) = [];
    % train linear regression
    mdl = fitlm(Xtrain,ytrain,'RobustOpts','off');
    betas = mdl.Coefficients.Estimate;
    ypred{i} = betas(1) + Xtest*betas(2:end);
    mse(i,:) = mean((ytest - ypred{i}).^2);
    fvu(i,:) = mse(i,:) ./ var(ytest);
end
% collect output
quants = quantile(fvu,3);
dataOut.FVUmed = quants(2);
dataOut.FVUupper = quants(3);
dataOut.FVUlower = quants(1);

quants = quantile(mse,3);
dataOut.AVUmed = quants(2);
dataOut.AVUupper = quants(3);
dataOut.AVUlower = quants(1);
for i = 1:size(leaveOut_bool,1)
    % find community index in testing subset
    if i > 1
        ypred_single = cellfun(@(x,y)x(y),ypred(leaveOut_bool(i,:)),num2cell(sum(leaveOut_bool(1:i,leaveOut_bool(i,:)))));
    else
        ypred_single = cellfun(@(x)x(i),ypred(leaveOut_bool(i,:)));
    end
    dataOut.ypredmed(i) = median(ypred_single);
    dataOut.ypredstd(i) = std(ypred_single);
end
end