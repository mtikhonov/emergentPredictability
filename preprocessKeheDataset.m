function [X,y,varnames] = preprocessKeheDataset(REPROCESS)
%% Clean dataset from Kehe etal 2019 for Moran & Tikhonov 2023
% this script loads the raw data published with Kehe etal 2019, cleans it
% for the present purposes and bins the full dataset into low- and
% high-diversity datasets
if exist('datasets/processed_kehe_dataset.mat','file') && ~REPROCESS
    fprintf('Dataset has already been preprocessed. Loading from disk...\n')
    load('datasets/processed_kehe_dataset.mat','X','y','varnames')
else
    %% Load raw data
    data = readtable('datasets\jared_data.csv');
    varnames = data(:,1:end-1).Properties.VariableNames;
    %% Split into diversity bins
    pres_abs_table = data{:,1:end-1};
    Y = data.fitness;
    trial_div = sum(pres_abs_table,2);
    low_div_bool = trial_div <= 5;
    high_div_bool = trial_div >= 10;

    y.lowdiversity = Y(low_div_bool);
    X.lowdiversity.pres_bool = pres_abs_table(low_div_bool,:);

    y.highdiversity = Y(high_div_bool);
    X.highdiversity.pres_bool = pres_abs_table(high_div_bool,:);
    %% Save processed dataset
    save('datasets/processed_kehe_dataset.mat','X','y','varnames')
end