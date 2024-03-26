function makeFigureCRMinvaders(num_invaders)
figure
clrs = [197/256, 90/256, 17/256;
        175/256, 170/256 170/256;
        143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));
W = 21.55; % two column figure
H = 3.5;
w = 1.85; h = 1.85;
dw = 0.2;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);

w0 = 1; h0 = 1;
invader_inds = 1:num_invaders;
for i = invader_inds
    load(['figureData/CRM_invader_' num2str(i) '_pareto_data.mat'],'CRM')
    ax = axes('Units','centimeters','Position',[w0 h0 w h]);
    box on
    hold on

    y2plt = CRM.middiv.FVUmedian;
    x2plt = CRM.middiv.desc_complexity;
    x2plt(x2plt<0) = 0;
    negerr2plt = y2plt - CRM.middiv.FVUlower;
    poserr2plt = CRM.middiv.FVUupper - y2plt;
    [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
    phmid = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
        'linewidth',1.5,'MarkerSize',14,'CapSize',2);
    set([phmid.Bar, phmid.Line], 'ColorType', 'truecoloralpha', 'ColorData', [phmid.Line.ColorData(1:3); 255*0.75])
    y2plt = CRM.lowdiv.FVUmedian;
    x2plt = CRM.lowdiv.desc_complexity;
    x2plt(x2plt<0) = 0;
    negerr2plt = y2plt - CRM.lowdiv.FVUlower;
    poserr2plt = CRM.lowdiv.FVUupper - y2plt;
    [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
    phmid = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
        'linewidth',2,'MarkerSize',16,'CapSize',2);
    set([phmid.Bar, phmid.Line], 'ColorType', 'truecoloralpha', 'ColorData', [phmid.Line.ColorData(1:3); 255*0.75])
    y2plt = CRM.highdiv.FVUmedian;
    x2plt = CRM.highdiv.desc_complexity;
    x2plt(x2plt<0) = 0;
    negerr2plt = y2plt - CRM.highdiv.FVUlower;
    poserr2plt = CRM.highdiv.FVUupper - y2plt;
    [x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
    phmid = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(3),:),'MarkerFaceColor',cmap(color_inds(3),:),...
        'linewidth',2,'MarkerSize',16,'CapSize',2);
    set([phmid.Bar, phmid.Line], 'ColorType', 'truecoloralpha', 'ColorData', [phmid.Line.ColorData(1:3); 255*0.75])
    ylim([0 1.1])
    xlim([0 5])
    yticks([0 0.5 1])
    xticks(ax,[0 2 4])

    if i < 2 %|| i == 6
        ylabel(ax,'FVU')
    end
    if i == 5 %|| i == 8
        xlabel(ax,'Description Complexity')
    end
    if i > 1 %&& i~=6
        yticklabels({})
    end
    title(ax,['Invader ' num2str(i)])
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