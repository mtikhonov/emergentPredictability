function chi = computeCost(g, costStructure)
% costStructure.J is interaction matrix (normally, zero diagonal)
% costStructure.chii is cost per pathway (must be a vector)
% costStructure.c is baseline cost for everyone
chi = costStructure.c + g*costStructure.chii(:) + sum(g.*(g*costStructure.J), 2);
end
