function generateGLVparetofront(invaderID)
if exist(['figureData/GLV_invader_' num2str(invaderID) '_pareto_data.mat'],'file')
    fprintf('Dataset has already been preprocessed. Loading from disk...\n')
else
    load(['figureData/GLV_invader_' num2str(invaderID) '.mat'],'data','params')
    % run Pareto front search algorithm
    % parpool(10);
    parfevalOnAll(@warning,0,'off','stats:LinearModel:RankDefDesignMat');
    
    % format data for proper inputs
    communitybool = false(1,params.S+1);
    communitybool(data.survivors_inds_full) = true;
    manual_grouping_info = [];
    X = data.lowdiv.abd_preinv(:,communitybool);
    y = data.lowdiv.abd_postinv(:,end);
    pres_abs_bool = data.lowdiv.present(:,communitybool);
    [GLV.lowdiv, GLV.lowdiv.conv_info] = ...
        getCGParetoFront(X,y,pres_abs_bool,'standardized',manual_grouping_info,'entropy');
    fprintf('Low Div done.\n')
    X = data.middiv.abd_preinv(:,communitybool);
    y = data.middiv.abd_postinv(:,end);
    pres_abs_bool = data.middiv.present(:,communitybool);
    [GLV.middiv, GLV.middiv.conv_info] = ...
        getCGParetoFront(X,y,pres_abs_bool,'standardized',manual_grouping_info,'entropy');
    fprintf('Mid Div done.\n')
    X = data.highdiv.abd_preinv(:,communitybool);
    y = data.highdiv.abd_postinv(:,end);
    pres_abs_bool = data.highdiv.present(:,communitybool);
    [GLV.highdiv, GLV.highdiv.conv_info] = ...
        getCGParetoFront(X,y,pres_abs_bool,'standardized',manual_grouping_info,'entropy');

    % save
    save(['figureData/GLV_invader_' num2str(invaderID) '_pareto_data.mat'], 'GLV')
end