function [paretoData,conv_info] = getCGParetoFront(X,y,pa_table,input_format,manual_info,complexity_opt)
% [group_assignment,FVUmedian,FVUlower,FVUupper,desc_complexity,conv_info] = getCGParetoFront(X,y,pa_table,input_format,manual_info,complexity_opt)
% get initial guess based on coefficients from fitting microscopic desc.
% X(X < 0) = 0;
paretoData = computeEIcurve(X,y,pa_table,input_format,manual_info,complexity_opt);
paretoData = sort_by_complexity(paretoData);
% evaluate all possible single-point relabelings
relabelData = relabelCGs(X,y,pa_table,input_format,paretoData.group_assignment,complexity_opt);
relabelData = sort_by_complexity(relabelData);
% check if any single-point relabelings are outside convex hull of initial
% guess
x_lower = paretoData.desc_complexity;
y_lower = paretoData.FVUlower;
[~,~,ic] = unique(x_lower); % find distinct x-axis values for interpolation
counts = accumarray(ic,1);
nonuniqgrp = find(counts > 1);
if ~isempty(nonuniqgrp)
    for i = 1:length(nonuniqgrp)
        y_lower(ic == nonuniqgrp(i)) = min(y_lower(ic == nonuniqgrp(i)));
    end
    [x_lower,ia,~] = unique(x_lower,'stable');
    y_lower = y_lower(ia);
end
fvu_interp = interp1(x_lower,y_lower,relabelData.desc_complexity); % lower bound = median - std
below_bool = fvu_interp > relabelData.FVUmedian'; % rmse of relabelled group assignment falls 1 std below line connecting current guess?
conv_info.numiterations = 0; % number of iterations took to converge
conv_info.num_added = 0; % number of new groupings added to Pareto front from relabeling step
conv_info.iter_thres = 50; % number of iterations until stopping loop and throwing warning did not converge (either too slow or not converging)
conv_info.converged_bool = false;
while any(below_bool) % repeat until convex hull doesn't change
    conv_info.numiterations = conv_info.numiterations + 1;
    if conv_info.numiterations > conv_info.iter_thres
        conv_info.converged_bool = false;
        break
    end
    % if yes, update current guess of Pareto front (add points that fall
    % below/outside convex hull)
    paretoData = updatePareto(paretoData,relabelData,below_bool);
    if length(unique(paretoData.desc_complexity)) > length(x_lower)
        num_pts_added = length(setdiff(unique(paretoData.desc_complexity),x_lower));
    else
        num_pts_added = length(setdiff(x_lower,unique(paretoData.desc_complexity)));
    end
    conv_info.num_added = conv_info.num_added + num_pts_added;
    % check again by perturbing updated Pareto front
    relabelData = relabelCGs(X,y,pa_table,input_format,paretoData.group_assignment,complexity_opt);
    relabelData = sort_by_complexity(relabelData);
    % find distinct x-axis values for interpolation
    x_lower = paretoData.desc_complexity;
    y_lower = paretoData.FVUlower;
    [~,~,ic] = unique(x_lower);
    counts = accumarray(ic,1);
    nonuniqgrp = find(counts > 1);
    if ~isempty(nonuniqgrp)
        for i = 1:length(nonuniqgrp)
            y_lower(ic == nonuniqgrp(i)) = min(y_lower(ic == nonuniqgrp(i)));
        end
        [x_lower,ia,~] = unique(x_lower,'stable');
        y_lower = y_lower(ia);
    end
    fvu_interp = interp1(x_lower,y_lower,relabelData.desc_complexity);
    below_bool = fvu_interp > relabelData.FVUmedian'; % rmse of relabelled group assignment falls below line connecting current guess?
    fprintf('Iter %d done\n',conv_info.numiterations)
end
if conv_info.numiterations < conv_info.iter_thres
    conv_info.converged_bool = true;
end
if ~conv_info.converged_bool
    warning(['Search algorithm did not converge within ' num2str(conv_info.iter_thres) ' iterations.'])
end
end

%% Helpers
function paretoStruct = sort_by_complexity(paretoStruct)
[paretoStruct.desc_complexity,ord] = sort(paretoStruct.desc_complexity);
paretoStruct.group_assignment = paretoStruct.group_assignment(:,ord);
paretoStruct.FVUmedian = paretoStruct.FVUmedian(ord); paretoStruct.FVUlower = paretoStruct.FVUlower(ord); paretoStruct.FVUupper = paretoStruct.FVUupper(ord);
paretoStruct.AVUmedian = paretoStruct.AVUmedian(ord); paretoStruct.AVUlower = paretoStruct.AVUlower(ord); paretoStruct.AVUupper = paretoStruct.AVUupper(ord);
paretoStruct.ypredmedian = paretoStruct.ypredmedian(:,ord); paretoStruct.ypredstd = paretoStruct.ypredstd(:,ord);
end

function pareto = updatePareto(pareto,relabelled,pts2add)
pareto.group_assignment = [pareto.group_assignment relabelled.group_assignment(:,pts2add)];
pareto.desc_complexity = [pareto.desc_complexity relabelled.desc_complexity(pts2add)];
pareto.FVUmedian = [pareto.FVUmedian;  relabelled.FVUmedian(pts2add)]; 
pareto.FVUlower= [pareto.FVUlower;  relabelled.FVUlower(pts2add)];
pareto.FVUupper = [pareto.FVUupper;  relabelled.FVUupper(pts2add)];
pareto.AVUmedian = [pareto.AVUmedian;  relabelled.AVUmedian(pts2add)]; 
pareto.AVUlower= [pareto.AVUlower;  relabelled.AVUlower(pts2add)];
pareto.AVUupper = [pareto.AVUupper;  relabelled.AVUupper(pts2add)];
pareto.ypredmedian = [pareto.ypredmedian relabelled.ypredmedian(:,pts2add)];
pareto.ypredstd = [pareto.ypredstd relabelled.ypredstd(:,pts2add)];
pareto = sort_by_complexity(pareto);
k = boundary(pareto.desc_complexity',pareto.FVUmedian,1);
k = k(1:find(pareto.desc_complexity(k) == max(pareto.desc_complexity),1)); % go until index hits micro entropy; lower boundary
pareto.group_assignment = pareto.group_assignment(:,k);
pareto.desc_complexity = pareto.desc_complexity(k);
pareto.FVUmedian = pareto.FVUmedian(k); pareto.FVUlower = pareto.FVUlower(k); pareto.FVUupper = pareto.FVUupper(k);
pareto.AVUmedian = pareto.AVUmedian(k); pareto.AVUlower = pareto.AVUlower(k); pareto.AVUupper = pareto.AVUupper(k);
pareto.ypredmedian = pareto.ypredmedian(:,k); pareto.ypredstd = pareto.ypredstd(:,k);
end



