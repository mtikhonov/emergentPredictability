function generateClarkKehePareto(options)
if strcmp(options.complexity_opt,'entropy') && strcmp(options.input_frmt,'standardized')
    fname = 'figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsEntropy.mat';
elseif strcmp(options.complexity_opt,'num_classes') && strcmp(options.input_frmt,'standardized')
    fname = 'figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsNumClasses.mat';
elseif strcmp(options.complexity_opt,'entropy') && strcmp(options.input_frmt,'raw')
    fname = 'figureData/clark_kehe_paretofronts_inputsRaw_complexityAsEntropy.mat';
end
if exist(fname,'file')
    fprintf('Dataset has already been preprocessed. Loading from disk...\n')
else
    %% Load processed data
    [X_clark,y_clark,varnames_clark,~] = preprocessClarkDataset(false);
    [X_kehe,y_kehe,~] = preprocessKeheDataset(false);
    %% Find Pareto Front Coarse-grainings and Evaluate Their Performance
    % NOTE: a design matrix rank deficiency warning is thrown by fitlm in cases
    % where a strain has a consistent input value (i.e., X(:,i) all same;
    % happens in low-diversity communities) and therefore results in a column
    % in design matrix that is linearly dependent with column representing the
    % intercept term (all 1s). Matlab handles these cases by setting corresponding
    % coefficient to 0, so we suppress this warning.
    parfevalOnAll(@warning,0,'off','stats:LinearModel:RankDefDesignMat');

    manual_grouping_info.varnames = varnames_clark;
    complexity_opt = options.complexity_opt;
    input_format = options.input_frmt;
    for var = 1:size(y_clark.lowdiversity,2)
        % low-diversity
        manual_grouping_info.varID = var; 
        manual_grouping_info.div = 'low';
        [paretoData,conv_info] = ...
            getCGParetoFront(X_clark.lowdiversity.absabd,y_clark.lowdiversity(:,var),X_clark.lowdiversity.pres_bool,input_format,manual_grouping_info,complexity_opt);
        clark.lowdiv{var} = paretoData;
        clark.lowdiv{var}.conv_info = conv_info;
        % mid-diversity
        manual_grouping_info.div = 'mid';
        [paretoData, conv_info] = ...
            getCGParetoFront(X_clark.middiversity.absabd,y_clark.middiversity(:,var),X_clark.middiversity.pres_bool,input_format,manual_grouping_info,complexity_opt);
        clark.middiv{var} = paretoData;
        clark.middiv{var}.conv_info = conv_info;
        % high-diversity
        manual_grouping_info.div = 'high';
        [paretoData,conv_info] = ...
            getCGParetoFront(X_clark.highdiversity.absabd,y_clark.highdiversity(:,var),X_clark.highdiversity.pres_bool,input_format,manual_grouping_info,complexity_opt);
        clark.highdiv{var} = paretoData;
        clark.highdiv{var}.conv_info = conv_info;
    end
    [paretoData,conv_info] = ...
        getCGParetoFront(X_kehe.lowdiversity.pres_bool,y_kehe.lowdiversity,X_kehe.lowdiversity.pres_bool,input_format,[],complexity_opt);
    kehe.lowdiv = paretoData;
    kehe.lowdiv.conv_info = conv_info;

    [paretoData,conv_info] = ...
        getCGParetoFront(X_kehe.highdiversity.pres_bool,y_kehe.highdiversity,X_kehe.highdiversity.pres_bool,input_format,[],complexity_opt);
    kehe.middiv = paretoData;
    kehe.middiv.conv_info = conv_info;
    if strcmp(options.complexity_opt,'num_classes') && strcmp(options.input_frmt,'standardized')
        %% Run tests on low-diveristy
        for var = 1:size(y_clark.lowdiversity,2)
            % "relabel" test
            [~,ord] = sort(clark.lowdiv{var}.desc_complexity);
            sorted_grp_assign = clark.lowdiv{var}.group_assignment(:,ord);
            relabelledPareto = ...
                relabelCGs(X_clark.lowdiversity.absabd,y_clark.lowdiversity(:,var),X_clark.lowdiversity.pres_bool,input_format,sorted_grp_assign(:,2:end-1),complexity_opt);
            clark.lowdiv{var}.relabeltest = relabelledPareto;
            % down-sampling test
            subsampledPareto = ...
                downsampleTest(X_clark.lowdiversity.absabd,y_clark.lowdiversity(:,var),X_clark.lowdiversity.pres_bool,X_clark.highdiversity.pres_bool,input_format,complexity_opt);
            clark.lowdiv{var}.subsampletest = subsampledPareto;
        end
        % "relabeling" test
        [~,ord] = sort(kehe.lowdiv.desc_complexity);
        sorted_grp_assign = kehe.lowdiv.group_assignment(:,ord);
        relabelledPareto = ...
            relabelCGs(X_kehe.lowdiversity.pres_bool,y_kehe.lowdiversity,X_kehe.lowdiversity.pres_bool,input_format,sorted_grp_assign(:,2:end-1),complexity_opt);
        kehe.lowdiv.relabeltest = relabelledPareto;
        % down-sampling test
        subsampledPareto = ...
            downsampleTest(X_kehe.lowdiversity.pres_bool,y_kehe.lowdiversity,X_kehe.lowdiversity.pres_bool,X_kehe.highdiversity.pres_bool,input_format,complexity_opt);
        kehe.lowdiv.subsampletest = subsampledPareto;
        %% Run tests on high-diversity
        for var = 1:size(y_clark.lowdiversity,2)
            % shuffle observable
            clark.highdiv{var}.permuteY = ...
                permutationTest(X_clark.highdiversity.absabd,y_clark.highdiversity(:,var),X_clark.highdiversity.pres_bool,input_format,complexity_opt,'permute_y');
            % shuffle abundance table
            clark.highdiv{var}.permuteX = ...
                permutationTest(X_clark.highdiversity.absabd,y_clark.highdiversity(:,var),X_clark.highdiversity.pres_bool,input_format,complexity_opt,'permute_Xclark');
        end
        % shuffle observable values
        kehe.middiv.permuteY = ...
            permutationTest(X_kehe.highdiversity.pres_bool,y_kehe.highdiversity,X_kehe.highdiversity.pres_bool,input_format,complexity_opt,'permute_y');
        % shuffle abundance table
        kehe.middiv.permuteX = ...
            permutationTest(X_kehe.highdiversity.pres_bool,y_kehe.highdiversity,X_kehe.highdiversity.pres_bool,input_format,complexity_opt,'permute_Xkehe');
    end
    %% Clean-up and save
    clark.X = X_clark; clark.y = y_clark;
    kehe.X = X_kehe; kehe.y = y_kehe;
    save(fname,'clark','kehe')
end
