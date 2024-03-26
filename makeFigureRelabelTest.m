function makeFigureRelabelTest
%% load clark etal and kehe etal data
load('figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsEntropy.mat')
%%
figure
clf
W = 21.55; % two column figure
H = 5.5;
w = 3.65; h = 3.65;
dw = 0.5;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);
clrs = [197/256, 90/256, 17/256;
        175/256, 170/256 170/256;
        143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));

observables = {'Butyrate', 'Acetate','Lactate','Succinate','Focal Strain Abundance'};
w0 = 1; h0 = 1;
for obs = 1:length(observables)
    ax = axes('Units','centimeters','Position',[w0 h0 w h]);
    box on
    hold on

    if obs < length(observables)
        thresh4vis = min(clark.lowdiv{obs}.FVUmedian) + 0.2;
        x2plt = clark.lowdiv{obs}.relabeltest.desc_complexity;
        y2plt = clark.lowdiv{obs}.relabeltest.FVUmedian;
        % negerr2plt = y2plt - clark.lowdiv{obs}.relabeltest.FVUlower;
        % poserr2plt = clark.lowdiv{obs}.relabeltest.FVUupper - y2plt;
        keep = y2plt < thresh4vis;
        y2plt(~keep) = []; x2plt(~keep) = []; %negerr2plt(~keep) = []; poserr2plt(~keep) = [];
        scatter(x2plt,y2plt,20,cmap(color_inds(1),:),'filled','Marker','s','MarkerFaceAlpha',0.35);
        % low-div Pareto front for reference
        y2plt = clark.lowdiv{obs}.FVUmedian;
        x2plt = clark.lowdiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        negerr2plt = y2plt - clark.lowdiv{obs}.FVUlower;
        poserr2plt = clark.lowdiv{obs}.FVUupper - y2plt;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
        ph1 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
            'linewidth',1.5,'MarkerSize',12,'CapSize',2);
        ylim([0 1.1])

        % high-div Pareto front for reference (no error-bars)
        y2plt = clark.highdiv{obs}.FVUmedian;
        x2plt = clark.highdiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord);
        ph3 = plot(ax,x2plt,y2plt,':', 'Color',cmap(color_inds(3),:),'MarkerFaceColor',cmap(color_inds(3),:),...
            'linewidth',1.5);

        yticks(ax,[0 0.5 1])
        xticks(ax,[0 2 4])
        if obs == 3
            xlabel(ax,'Description Complexity')
        end
        if obs < 2
            ylabel(ax,'FVU')
        end
        title(ax,observables{obs})
    else
        thresh4vis = min(kehe.lowdiv.FVUmedian) + 0.2;
        x2plt = kehe.lowdiv.relabeltest.desc_complexity;
        y2plt = kehe.lowdiv.relabeltest.FVUmedian;
        % negerr2plt = y2plt - kehe.lowdiv.relabeltest.FVUlower;
        % poserr2plt = kehe.lowdiv.relabeltest.FVUupper - y2plt;
        keep = y2plt < thresh4vis;
        y2plt(~keep) = []; x2plt(~keep) = []; %err2plt(~keep) = [];
        scatter(x2plt,y2plt,20,cmap(color_inds(1),:),'filled','Marker','s','MarkerFaceAlpha',0.35);
        % add low-div Pareto front
        y2plt = kehe.lowdiv.FVUmedian;
        negerr2plt = y2plt - kehe.lowdiv.FVUlower;
        poserr2plt = kehe.lowdiv.FVUupper - y2plt;
        errorbar(kehe.lowdiv.desc_complexity,y2plt,negerr2plt,poserr2plt,...
            '.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
            'linewidth',1.5,'MarkerSize',12,'CapSize',2)
        ylim([0 1.1])

        y2plt = kehe.middiv.FVUmedian;
        plot(ax,kehe.middiv.desc_complexity,y2plt,...
            ':', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
            'linewidth',1.5)

        yticks(ax,[0 0.5 1])
        xticks(ax,[0 2 4])
        % xlabel(ax,'Description Complexity')
        title(ax,observables{obs})
    end
    if obs > 1
        yticklabels({})
    end

    [w0,h0] = getPanelCoords(w0,h0,w,dw,0);
end
end

%%
% get panel dimensions
function [x0,y0] = getPanelCoords(xc,yc,w,dw,dh)
if nargin < 5
    x0 = xc + w + dw;
    y0 = yc;
else
    x0 = xc + w + dw;
    y0 = yc - dh;
end
end