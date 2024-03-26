function J = interactionMatrix_sigmoid(L, nstar, delta, seed)
% This is for pairwise INTERACTIONS (no diagonal terms)
%
% generates a triangular L by L matrix with zero diagonal
% whose entries fall off in mean as mu*N^-alpha, 
% and in standard deviation as sigma*N^-beta
n = max(1:L, (1:L)');
sd = 1./(1+exp((n-nstar)/delta));

half = triu(true(L),1);
% Using triu instead of tril is cruicial!
% Using triu set the matrix values such that, at the same seed, increasing 
% L would keep the exact same values and add to them (rather than resetting
% the entire matrix). That way increasing L would be like a more refined
% view of SAME biochemistry.

sd = sd(half);
rng(seed,'combRecursive');
J = double(half);
J(half) = sd.*randn(size(sd));
end