function makefigureDownSample
%% load clark etal data
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
color_inds = floor(linspace(1,150,20));

observables = {'Butyrate', 'Acetate','Lactate','Succinate','Focal Strain Abundance'};
w0 = 1; h0 = 1;
for obs = 1:length(observables)
    ax = axes('Units','centimeters','Position',[w0 h0 w h]);
    box on
    hold on
    
    if obs < length(observables)
        xs2plt = cellfun(@(x)x.desc_complexity,clark.lowdiv{obs}.subsampletest,'un',0);
        ys2plt = cellfun(@(x)x.FVUmedian,clark.lowdiv{obs}.subsampletest,'un',0);
        negerrs2plt = cellfun(@(x)x.FVUmedian - x.FVUlower,clark.lowdiv{obs}.subsampletest,'un',0);
        poserrs2plt = cellfun(@(x)x.FVUupper - x.FVUmedian,clark.lowdiv{obs}.subsampletest,'un',0);
        for sample = 1:length(xs2plt)
            x2plt = xs2plt{sample};
            y2plt = ys2plt{sample};
            negerr2plt = negerrs2plt{sample};
            poserr2plt = poserrs2plt{sample};
            [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
            erh = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'<-',...
                'linewidth',1,'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),'MarkerSize',3,'CapSize',2);
            % Set transparency (undocumented)
            alpha = 0.25;
            set([erh.Bar, erh.Line], 'ColorType', 'truecoloralpha', 'ColorData', [erh.Line.ColorData(1:3); 255*alpha])
        end
        % low-div errorbar for reference
        y2plt = clark.lowdiv{obs}.FVUmedian;
        x2plt = clark.lowdiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        negerr2plt = y2plt - clark.lowdiv{obs}.FVUlower;
        poserr2plt = clark.lowdiv{obs}.FVUupper - y2plt;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
        ph1 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
            'linewidth',1.5,'MarkerSize',12,'CapSize',2);
        ylim([0 1.1])

        % high-div line for reference (no error-bars)
        y2plt = clark.highdiv{obs}.FVUmedian;
        x2plt = clark.highdiv{obs}.desc_complexity;
        x2plt(x2plt<0) = 0;
        [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord);
        ph3 = plot(ax,x2plt,y2plt,':', 'Color',cmap(color_inds(end),:),'MarkerFaceColor',cmap(color_inds(end),:),...
            'linewidth',1.5);

        yticks(ax,[0 0.5 1])
        xticks(ax,[0 2 4])
        % xlabel(ax,'Description Complexity')
        if obs < 2
            ylabel(ax,'FVU')
        end
        if obs == 3
            xlabel(ax,'Description Complexity')
        end
        title(ax,observables{obs})
    else
        xs2plt = cellfun(@(x)x.desc_complexity,kehe.lowdiv.subsampletest,'un',0);
        ys2plt = cellfun(@(x)x.FVUmedian,kehe.lowdiv.subsampletest,'un',0);
        negerrs2plt = cellfun(@(x)x.FVUmedian - x.FVUlower,kehe.lowdiv.subsampletest,'un',0);
        poserrs2plt = cellfun(@(x)x.FVUupper - x.FVUmedian,kehe.lowdiv.subsampletest,'un',0);
        for sample = 1:length(xs2plt)
            x2plt = xs2plt{sample};
            y2plt = ys2plt{sample};
            negerr2plt = negerrs2plt{sample}; poserr2plt = poserrs2plt{sample};
            erh = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'<-',...
                'linewidth',1,'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),'MarkerSize',3,'CapSize',2);
            % Set transparency (undocumented)
            alpha = 0.25;
            set([erh.Bar, erh.Line], 'ColorType', 'truecoloralpha', 'ColorData', [erh.Line.ColorData(1:3); 255*alpha])
        end
        y2plt = kehe.lowdiv.FVUmedian;
        negerr2plt = kehe.lowdiv.FVUmedian - kehe.lowdiv.FVUlower;
        poserr2plt = kehe.lowdiv.FVUupper - kehe.lowdiv.FVUmedian;
        errorbar(kehe.lowdiv.desc_complexity,y2plt,negerr2plt,poserr2plt,...
            '.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
            'linewidth',1.5,'MarkerSize',12,'CapSize',2)
        ylim([0 1.1])

        y2plt = kehe.middiv.FVUmedian;
        plot(ax,kehe.middiv.desc_complexity,y2plt,...
            ':', 'Color',cmap(color_inds(11),:),'MarkerFaceColor',cmap(color_inds(11),:),...
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