function chis = computeMuts(mutG, costStructure)
L = size(mutG,2);
muts = xor(repmat(mutG,[L,1]), eye(L));
chis = computeCost(muts, costStructure)';
end
