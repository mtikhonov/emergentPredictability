function makeFigureNclasses
options.complexity_opt = 'num_classes';
options.input_frmt = 'standardized';
generateClarkKehePareto(options);
plotFigure;
end

function plotFigure
load('figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsNumClasses.mat')
%%
figure
clf
W = 21.55; % two column figure
H = 5.75;
w = 3.65; h = 3.65;
dw = 0.5;
font_sz = 11; lw = 2; mkrsz = 16;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);
clrs = [197/256, 90/256, 17/256;
        175/256, 170/256 170/256;
        143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));

observables = {'Butyrate', 'Acetate','Lactate','Succinate','Focal Strain Abundance'};
w0 = 1; h0 = 1.5;
for obs = 1:length(observables)
    ax = axes('Units','centimeters','Position',[w0 h0 w h]);
    box on
    hold on
    
    if obs < length(observables)
        y2plt = clark.middiv{obs}.FVUmedian;
        x2plt = clark.middiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        negerr2plt = y2plt - clark.middiv{obs}.FVUlower;
        poserr2plt = clark.middiv{obs}.FVUupper - y2plt;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
        phmid = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
            'linewidth',lw-0.5,'MarkerSize',mkrsz-2,'CapSize',2);
        set([phmid.Bar, phmid.Line], 'ColorType', 'truecoloralpha', 'ColorData', [phmid.Line.ColorData(1:3); 255*0.75])

        % low-div errorbar for reference
        y2plt = clark.lowdiv{obs}.FVUmedian;
        x2plt = clark.lowdiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        negerr2plt = y2plt - clark.lowdiv{obs}.FVUlower;
        poserr2plt = clark.lowdiv{obs}.FVUupper - y2plt;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
        ph1 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
            'linewidth',lw,'MarkerSize',12,'CapSize',2);
        ylim([0 1.1])

        % high-div line for reference (no error-bars)
        y2plt = clark.highdiv{obs}.FVUmedian;
        x2plt = clark.highdiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        negerr2plt = y2plt - clark.highdiv{obs}.FVUlower;
        poserr2plt = clark.highdiv{obs}.FVUupper - y2plt;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
        ph2 = errorbar(ax,x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(3),:),'MarkerFaceColor',cmap(color_inds(3),:),...
            'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);

        yticks(ax,[0 0.5 1])
        xlim([0 26])
        xticks(ax,[0 5 10 15 20 25])
        % xlabel(ax,'Description Complexity')
        if obs < 2
            ylabel(ax,'FVU')
            lhand = legend(ax,[ph1 phmid ph2],{'Low','Mid','High'},'Location','northeast','FontSize',font_sz-2);
            lhand.Position(2) = lhand.Position(2) - 0.04;
            text(lhand.Position(1)+13,lhand.Position(2)+0.45,lhand.Position(3),'Richness','FontSize',font_sz-2)
        end
        if obs == 3
            xlabel(ax,{'Description Complexity', '[Number of Classes]'})
        end
        title(ax,observables{obs})
    else

        y2plt = kehe.lowdiv.FVUmedian;
        negerr2plt = kehe.lowdiv.FVUmedian - kehe.lowdiv.FVUlower;
        poserr2plt = kehe.lowdiv.FVUupper - kehe.lowdiv.FVUmedian;
        errorbar(kehe.lowdiv.desc_complexity,y2plt,negerr2plt,poserr2plt,...
            '.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
            'linewidth',lw,'MarkerSize',12,'CapSize',2)
        ylim([0 1.1])

        x2plt = kehe.middiv.desc_complexity;
        y2plt = kehe.middiv.FVUmedian;
        negerr2plt = y2plt - kehe.middiv.FVUlower;
        poserr2plt = kehe.middiv.FVUupper - y2plt;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
        ph2 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,...
            '.-', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
            'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);

        yticks(ax,[0 0.5 1])
        xticks(ax,[0 5 10 15 20])
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