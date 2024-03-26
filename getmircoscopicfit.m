function [coeffs,coeffs_SE] = getmircoscopicfit(X,y)
nVars = size(y,2);
for j = nVars:-1:1
    mdl = fitlm(X,y(:,j),'RobustOpts','off');
    coeffs(:,j) = mdl.Coefficients.Estimate;
    coeffs_SE(:,j) = mdl.Coefficients.SE;
end