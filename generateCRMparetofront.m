function generateCRMparetofront(focal_id)
if exist(['figureData/CRM_invader_' num2str(focal_id) '_pareto_data.mat'],'file')
    fprintf('Dataset has already been preprocessed. Loading from disk...\n')
else
    load(['/data/jmoran/figureData/CRM_invader_' num2str(focal_id) '.mat'],'dataCRM')
    % run Pareto front search algorithm
    parfevalOnAll(@warning,0,'off','stats:LinearModel:RankDefDesignMat');
    
    % format data for proper inputs
    manual_grouping_info = [];
    X = dataCRM.lowdiv.X; X(X<1e-6) = 0;
    y = dataCRM.lowdiv.postinv_abd; y(y<1e-6) = 0;
    pres_abs_bool = dataCRM.lowdiv.pres_abs_bool;
    [CRM.lowdiv, CRM.lowdiv.conv_info] = ...
        getCGParetoFront(X,y,pres_abs_bool,'standardized',manual_grouping_info,'entropy');
    fprintf('Low Div done.\n')
    X = dataCRM.middiv.X; X(X<1e-6) = 0;
    y = dataCRM.middiv.postinv_abd;
    pres_abs_bool = dataCRM.middiv.pres_abs_bool;
    [CRM.middiv, CRM.middiv.conv_info] = ...
        getCGParetoFront(X,y,pres_abs_bool,'standardized',manual_grouping_info,'entropy');
    fprintf('Mid Div done.\n')
    X = dataCRM.highdiv.X; X(X<1e-6) = 0;
    y = dataCRM.highdiv.postinv_abd; y(y<1e-6) = 0;
    pres_abs_bool = dataCRM.highdiv.pres_abs_bool;
    [CRM.highdiv, CRM.highdiv.conv_info] = ...
        getCGParetoFront(X,y,pres_abs_bool,'standardized',manual_grouping_info,'entropy');Hi
    fprintf('High Div done.\n')

    % save
    save(['figureData/CRM_invader_' num2str(focal_id) '_pareto_data.mat'],'CRM')
end
end