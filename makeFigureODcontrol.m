function makeFigureODcontrol
%% Load clark etal data
load('figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsEntropy.mat')
%%
figure
clf
font_sz = 11; lw = 2; mkrsz = 16;
W = 8.5;
H = 5.5;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);
clrs = [197/256, 90/256, 17/256;
        175/256, 170/256 170/256;
        143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));
x0 = 1.1; y0 = 1.5; w = 3.1; h = 3.25; dw = 0.75;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
obsID = 5; % total OD
ax.FontSize = font_sz-2;
box on
hold on

y2plt = clark.middiv{obsID}.FVUmedian;
x2plt = clark.middiv{obsID}.desc_complexity;
x2plt(x2plt<0) = 0;
negerr2plt = y2plt - clark.middiv{obsID}.FVUlower;
poserr2plt = clark.middiv{obsID}.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
phmid = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
        'linewidth',lw-0.5,'MarkerSize',mkrsz-2,'CapSize',2);
set([phmid.Bar, phmid.Line], 'ColorType', 'truecoloralpha', 'ColorData', [phmid.Line.ColorData(1:3); 255*0.75])

y2plt = clark.lowdiv{obsID}.FVUmedian;
x2plt = clark.lowdiv{obsID}.desc_complexity;
x2plt(x2plt<0) = 0;
negerr2plt = y2plt - clark.lowdiv{obsID}.FVUlower;
poserr2plt = clark.lowdiv{obsID}.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
ph1 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
        'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);
ylim([0 1.1])
% xlim([0 26])

y2plt = clark.highdiv{obsID}.FVUmedian;
x2plt = clark.highdiv{obsID}.desc_complexity;
x2plt(x2plt<0) = 0;
negerr2plt = y2plt - clark.highdiv{obsID}.FVUlower;
poserr2plt = clark.highdiv{obsID}.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
ph2 = errorbar(ax,x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(3),:),'MarkerFaceColor',cmap(color_inds(3),:),...
        'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);

yticks(ax,[0 0.5 1])
ylabel(ax,'FVU (MSE/var)','FontSize',font_sz-1)
title(ax,'Total OD','FontSize',font_sz)
% add common legend for panels A and B
lhand = legend(ax,[ph1 phmid ph2],{'Low','Mid','High'},'Location','northeast','FontSize',font_sz-2);
lhand.Position(2) = lhand.Position(2) - 0.04;
text(lhand.Position(1)+1.8,lhand.Position(2)+0.525,lhand.Position(3),'Richness','FontSize',font_sz-2)
xlabel({'Description Complexity', '[Partition Entropy]'})
%%
load('figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsNumClasses.mat')
%%
x0 = x0 + w + dw;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
ax.FontSize = font_sz-2;
cla
box on
hold on
y2plt = clark.middiv{obsID}.FVUmedian;
x2plt = clark.middiv{obsID}.desc_complexity;
x2plt(x2plt<0) = 0;
negerr2plt = y2plt - clark.middiv{obsID}.FVUlower;
poserr2plt = clark.middiv{obsID}.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
phmid = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
        'linewidth',lw-0.5,'MarkerSize',mkrsz-2,'CapSize',2);
set([phmid.Bar, phmid.Line], 'ColorType', 'truecoloralpha', 'ColorData', [phmid.Line.ColorData(1:3); 255*0.75])

y2plt = clark.lowdiv{obsID}.FVUmedian;
x2plt = clark.lowdiv{obsID}.desc_complexity;
x2plt(x2plt<0) = 0;
negerr2plt = y2plt - clark.lowdiv{obsID}.FVUlower;
poserr2plt = clark.lowdiv{obsID}.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
ph1 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
        'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);
ylim([0 1.1])
xlim([0 26])

y2plt = clark.highdiv{obsID}.FVUmedian;
x2plt = clark.highdiv{obsID}.desc_complexity;
x2plt(x2plt<0) = 0;
negerr2plt = y2plt - clark.highdiv{obsID}.FVUlower;
poserr2plt = clark.highdiv{obsID}.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
ph2 = errorbar(ax,x2plt,y2plt,negerr2plt,poserr2plt,'.-', 'Color',cmap(color_inds(3),:),'MarkerFaceColor',cmap(color_inds(3),:),...
        'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);

yticks(ax,[0 0.5 1])
yticklabels({})
title(ax,'Total OD','FontSize',font_sz)
xlabel({'Description Complexity', '[Number of Classes]'})
end