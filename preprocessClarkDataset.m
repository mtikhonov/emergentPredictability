function [X,y,varnames,pHdata,pptynames] = preprocessClarkDataset(REPROCESS)
%% Clean dataset from Clark etal 2021 for Moran & Tikhonov 2023
% this script loads the raw data published with Clark etal 2021, cleans it
% for the present purposes and bins the full dataset into low-, mid- and
% high-diversity datasets
if exist('datasets/processed_clark_dataset.mat','file') && ~REPROCESS
    fprintf('Dataset has already been preprocessed. Loading from disk...\n')
    load('datasets/processed_clark_dataset.mat','X','y','varnames','pHdata','pptynames')
else
    %% Load raw dataset
    abdData_master = readtable('datasets/2020_02_28_MasterDF.csv');
    %% Use only non-contaminated trials
    not_contam_bool = ismember(abdData_master.Contamination_,'No');
    abdData_nocontam = abdData_master(not_contam_bool,:);
    %% Use only trials that do not possess strain "HB"
    abdData_nocontam.HB(isnan(abdData_nocontam.HB)) = 0;
    HB_pa_bool = logical(abdData_nocontam.HB);
    abdData_nocontam_noHB = abdData_nocontam(~HB_pa_bool,:);
    %% Find communities for which pH was measured
    pH_measured_bool = ~isnan(abdData_nocontam_noHB.pH);
    %% Bin data into low-, mid-, high-diversity
    strainsN_pertrial = sum(abdData_nocontam_noHB{:,9:34},2);

    low_div = (strainsN_pertrial <= 5) & (strainsN_pertrial >= 1);
    low_div_pH_bool = low_div & pH_measured_bool;
    low_div_pH = abdData_nocontam_noHB.pH(low_div_pH_bool);
    % low_div_obs = [abdData_nocontam_noHB.Butyrate(low_div_pH_bool) abdData_nocontam_noHB.Acetate(low_div_pH_bool) ...
                   % abdData_nocontam_noHB.Lactate(low_div_pH_bool) abdData_nocontam_noHB.Succinate(low_div_pH_bool)];
    low_div_pH_bool = ismember(find(low_div),find(pH_measured_bool));
    low_div = abdData_nocontam_noHB(low_div,:);
    low_div(:,[49 55 58 59 61 65:end]) = [];
    low_div.Var1 = [];
    low_div.Rep = [];
    low_div.HB =[];

    mid_div = (strainsN_pertrial <= 15) & (strainsN_pertrial > 10);
    mid_div_pH_bool = mid_div & pH_measured_bool;
    mid_div_pH = abdData_nocontam_noHB.pH(mid_div_pH_bool);
    % mid_div_obs = [abdData_nocontam_noHB.Butyrate(mid_div_pH_bool) abdData_nocontam_noHB.Acetate(mid_div_pH_bool) ...
                   % abdData_nocontam_noHB.Lactate(mid_div_pH_bool) abdData_nocontam_noHB.Succinate(mid_div_pH_bool)];
    mid_div_pH_bool = ismember(find(mid_div),find(pH_measured_bool));
    mid_div = abdData_nocontam_noHB(mid_div,:);
    mid_div(:,[49 55 58 59 61 65:end]) = [];
    mid_div.Var1 = [];
    mid_div.Rep = [];
    mid_div.HB =[];

    high_div = strainsN_pertrial >= 20;
    high_div_pH_bool = high_div & pH_measured_bool;
    high_div_pH = abdData_nocontam_noHB.pH(high_div_pH_bool);
    % high_div_obs = [abdData_nocontam_noHB.Butyrate(high_div_pH_bool) abdData_nocontam_noHB.Acetate(high_div_pH_bool) ...
                   % abdData_nocontam_noHB.Lactate(high_div_pH_bool) abdData_nocontam_noHB.Succinate(high_div_pH_bool)];
    high_div_pH_bool = ismember(find(high_div),find(pH_measured_bool));
    high_div = abdData_nocontam_noHB(high_div,:);
    high_div(:,[49 55 58 59 61 65:end]) = [];
    high_div.Var1 = [];
    high_div.Rep = [];
    high_div.HB =[];

    % abd variables and pres/abs variables have different ordering!
    varnames = cellfun(@(x)x(1:2),low_div.Properties.VariableNames(32:56),'un',0);
    pptynames = {'Butyrate','Acetate','Lactate','Succinate','OD'};
    spsnames = low_div.Properties.VariableNames(7:31);
    %% Format data for inputs
    % reorder for pres/abs multiplication
    pa_2_odfrac = nan(1,25);
    for i = 1:25
        pa_2_odfrac(find(ismember(varnames,spsnames{i}))) = i;
    end
    pa_low = low_div{:,7:31};
    pa_low = pa_low(:,pa_2_odfrac);
    pa_mid = mid_div{:,7:31};
    pa_mid = pa_mid(:,pa_2_odfrac);
    pa_high = high_div{:,7:31};
    pa_high = pa_high(:,pa_2_odfrac);
    % get relative, absolute, and total biomass
    X.lowdiversity.pres_bool = pa_low;
    X.lowdiversity.ODtot = low_div.OD;
    X.lowdiversity.relabd = low_div{:,32:end} .* pa_low;
    X.lowdiversity.absabd = X.lowdiversity.ODtot .* X.lowdiversity.relabd;
    y.lowdiversity = [low_div.(pptynames{1}) low_div.(pptynames{2}) low_div.(pptynames{3}) low_div.(pptynames{4}) low_div.(pptynames{5})];
    pHdata.lowdiversity.pH = low_div_pH;
    pHdata.lowdiversity.pH_bool = low_div_pH_bool;
    % pHdata.lowdiversity.y = low_div_obs;

    X.middiversity.pres_bool = pa_mid;
    X.middiversity.ODtot = mid_div.OD;
    X.middiversity.relabd = mid_div{:,32:end} .* pa_mid;
    X.middiversity.absabd = X.middiversity.ODtot .* X.middiversity.relabd;
    y.middiversity = [mid_div.(pptynames{1}) mid_div.(pptynames{2}) mid_div.(pptynames{3}) mid_div.(pptynames{4}) mid_div.(pptynames{5})];
    pHdata.middiversity.pH = mid_div_pH;
    pHdata.middiversity.pH_bool = mid_div_pH_bool;
    % pHdata.middiversity.y = mid_div_obs;

    X.highdiversity.pres_bool = pa_high;
    X.highdiversity.ODtot = high_div.OD;
    X.highdiversity.relabd = high_div{:,32:end} .* pa_high;
    X.highdiversity.absabd = X.highdiversity.ODtot .* X.highdiversity.relabd;
    y.highdiversity = [high_div.(pptynames{1}) high_div.(pptynames{2}) high_div.(pptynames{3}) high_div.(pptynames{4}) high_div.(pptynames{5})];
    pHdata.highdiversity.pH = high_div_pH;
    pHdata.highdiversity.pH_bool = high_div_pH_bool;
    % pHdata.highdiversity.y = high_div_obs;
    %% Save preprocessed dataset
    save('datasets/processed_clark_dataset.mat','X','y','pHdata','varnames','pptynames')
end