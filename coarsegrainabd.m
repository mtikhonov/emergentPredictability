function [cgabd,entropy] = coarsegrainabd(abd,groupassign)
% Input:
%   - abd: 2-D array of strain abundances with columns corresponding to
%          different strains and rows corresponding to different
%          observations
%   - groupassign: 1-D array assigning strains to groups
% Output:
%   - cgabd: combined abundance of strains in coarse-grained group
%   - entropy: information content of coarse-graining
Nstrains = size(abd,2);
p = ones(1,Nstrains)/Nstrains; % assuming strains have equal probability/weight
for i = max(groupassign):-1:1
    cgabd(:,i) = sum(abd(:,groupassign == i),2);
    cgp(i) = sum(p(groupassign == i));
end
entropy = -sum(cgp.*log2(cgp));
end