function makeFigure3
options.complexity_opt = 'entropy';
options.input_frmt = 'standardized';
generateClarkKehePareto(options);
plotFigure;
end

function plotFigure
load('figureData/clark_kehe_paretofronts_inputsStandardized_complexityAsEntropy.mat')
%% butyrate
figure
clf
font_sz = 11; lw = 2; mkrsz = 16;
W = 21.55; % two column figure
H = 5;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);
% cmap = colormap(gray);
% color_inds = floor(linspace(1,150,3));
clrs = [197/256, 90/256, 17/256;
        175/256, 170/256 170/256;
        143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));
x0 = 1.1; y0 = 1.1; w = 3.1; h = 3.25; dw = 0.22;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
obsID = 1;
% ax = subplot(2,3,obsID+1);
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
% xticks(ax,[0 2 4])
% xlabel(ax,{'Description', 'Complexity'},'FontSize',font_sz-1)
% ylabel(ax,'Unexplained Variance','FontSize',font_sz-1)
ylabel(ax,'FVU (MSE/var)','FontSize',font_sz-1)
title(ax,'Butyrate','FontSize',font_sz)
% add common legend for panels A-E
lhand = legend(ax,[ph1 phmid ph2],{'Low','Mid','High'},'Location','northeast','FontSize',font_sz-2);
lhand.Position(2) = lhand.Position(2) - 0.04;
text(lhand.Position(1)+1.8,lhand.Position(2)+0.525,lhand.Position(3),'Richness','FontSize',font_sz-2)
% lhand = legend(ax,[ph1 ph2],{'Low Diversity','High Diversity'},'Location','southeast','Orientation','Horizontal','FontSize',font_sz);
% lhand.Position(1) = lhand.Position(1) + 0.4; lhand.Position(2) = lhand.Position(2) - 0.4;
%% acetate
obsID = 2;
% ax = subplot(2,3,obsID+1);
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
yticklabels({})
% xticks(ax,[0 2 4])
% xlabel(ax,{'Description', 'Complexity'},'FontSize',font_sz-1)
% ylabel(ax,'Normalized Prediction Error','FontSize',font_sz)
title(ax,'Acetate','FontSize',font_sz)

%% lactate
obsID = 3;
% ax = subplot(2,3,obsID+1);
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
yticks({})
% xticks(ax,[0 2 4])
xlabel(ax,'Description Complexity','FontSize',font_sz-1)
% ylabel(ax,'Normalized Prediction Error','FontSize',font_sz)
title(ax,'Lactate','FontSize',font_sz)

%% succinate
obsID = 4;
% ax = subplot(2,3,obsID+1);
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
yticks({})
% xticks(ax,[0 2 4])
% xlabel(ax,{'Description', 'Complexity'},'FontSize',font_sz-1)
% ylabel(ax,'Normalized Prediction Error','FontSize',font_sz)
title(ax,'Succinate','FontSize',font_sz)

%% kehe
% ax = subplot(2,3,6);
x0 = x0 + w + dw;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
ax.FontSize = font_sz-2;
cla
box on
hold on
x2plt = kehe.lowdiv.desc_complexity;
y2plt = kehe.lowdiv.FVUmedian;
negerr2plt = y2plt - kehe.lowdiv.FVUlower;
poserr2plt = kehe.lowdiv.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
p1h = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,...
    '.-', 'Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),...
        'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);
ylim([0 1.1])
% xlim([0 15])

x2plt = kehe.middiv.desc_complexity;
y2plt = kehe.middiv.FVUmedian;
negerr2plt = y2plt - kehe.middiv.FVUlower;
poserr2plt = kehe.middiv.FVUupper - y2plt;
[x2plt,ord] = sort(x2plt); y2plt = y2plt(ord); negerr2plt = negerr2plt(ord); poserr2plt = poserr2plt(ord);
ph2 = errorbar(x2plt,y2plt,negerr2plt,poserr2plt,...
    '.-', 'Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),...
        'linewidth',lw,'MarkerSize',mkrsz,'CapSize',2);

yticks(ax,[0 0.5 1])
yticks({})
% xticks(ax,[0 2 4])
% xlabel(ax,{'Description', 'Complexity'},'FontSize',font_sz-1)
% ylabel(ax,'Normalized Prediction Error','FontSize',font_sz)
title(ax,'Strain Abundance','FontSize',font_sz)
% lhand = legend(ax,[ph1 ph2],{'Low Ricness','High Richness'},'Location','southeast','Orientation','Horizontal','FontSize',font_sz-2);
% lhand.Position(1) = lhand.Position(1) - 0.22; lhand.Position(2) = lhand.Position(2) - 0.02;

%% Panels G & H
obsID = 1; % butyrate
x0 = x0 + w + dw + 1.5;
y0 = 3.1;
h = 1.6; w = 1.6;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
ax.FontSize = font_sz-2;
box on
hold on
errorbar(clark.y.lowdiversity(:,obsID),clark.lowdiv{obsID}.ypredmedian(:,2),clark.lowdiv{obsID}.ypredstd(:,2),'.',...
    'Color',cmap(color_inds(1),:),'MarkerSize',10,'CapSize',2)
maxVal = max([clark.y.lowdiversity(:,obsID); clark.lowdiv{obsID}.ypredmedian(:,2)]);
plot([0, maxVal], [0, maxVal], 'k--')
xticks([0 25 50])
yticks([0 25 50])
xticklabels({})
ylabel({'Predicted', '[Butyrate]'},'FontSize',font_sz-1)
% xlabel('Standardized Group Abundance')
thandle = title('Low Richness','FontSize',font_sz-2);
thandle.Position(2) = thandle.Position(2) - 2.75;
% thandle.Position(2) = thandle.Position(2)-0.1;
% set(gca,'FontSize',font_sz)

%%
y0 = 1; 
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
ax.FontSize = font_sz-2;
box on
hold on
errorbar(clark.y.highdiversity(:,obsID),clark.highdiv{obsID}.ypredmedian(:,2),clark.highdiv{obsID}.ypredstd(:,2),'.',...
    'Color',cmap(color_inds(3),:),'MarkerSize',10,'CapSize',2)
plot([0, maxVal], [0, maxVal], 'k--')
xlim([0 maxVal])
ylim([0 maxVal])
xticks([0 25 50])
yticks([0 25 50])
xticklabels({})
ylabel({'Predicted', '[Butyrate]'},'FontSize',font_sz-1)
% xhand = xlabel({'Standardized', 'Group Abundance'},'FontSize',font_sz-1);
xhand = xlabel('Measured [Butyrate]','FontSize',font_sz-1);
xhand.Position(2) = xhand.Position(2) - 10;
xhand.Position(1) = xhand.Position(1) - 5;
thandle = title('High Richness','FontSize',font_sz-2);
thandle.Position(2) = thandle.Position(2) - 2.25;
end