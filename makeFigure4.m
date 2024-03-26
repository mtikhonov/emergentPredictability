function makeFigure4
%% load Clark etal dataset
load('datasets\processed_clark_dataset.mat')
%% standardize all
% only look at functional properties (fermentation end productions)
ppty2incl = ismember(pptynames,{'Butyrate','Acetate','Lactate','Succinate'});
y_all = [y.lowdiversity(:,ppty2incl); y.middiversity(:,ppty2incl); y.highdiversity(:,ppty2incl)];
lowdiv_bool = ismember(y_all,y.lowdiversity(:,ppty2incl),'rows');
middiv_bool = ismember(y_all,y.middiversity(:,ppty2incl),'rows');
highdiv_bool = ismember(y_all,y.highdiversity(:,ppty2incl),'rows');
mu_all = mean(y_all);
sig_all = std(y_all);
y_all_strdrd = (y_all - mu_all) ./ sig_all;
%% do PCA on high div
[~,projs,~,~,~] = pca(y_all_strdrd,'Centered',false);
%%
figure
clf
font_sz = 11; lw = 2; mkrsz = 27;
W = 16.51; % two column figure
H = 10.8;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);
x0 = 1.1; y0 = 1.1; w = 8.25; h = 8.55; dw = 2.25;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
hold on
box on
clrs = [197/256, 90/256, 17/256;
        175/256, 170/256 170/256;
        143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));
scatter(ax,projs(lowdiv_bool,1),projs(lowdiv_bool,2),mkrsz-5,cmap(color_inds(1),:),'filled')
% scatter(ax,projs(middiv_bool,1),projs(middiv_bool,2),mkrsz-5,cmap(color_inds(2),:),'filled')
scatter(ax,projs(highdiv_bool,1),projs(highdiv_bool,2),mkrsz+10,cmap(color_inds(3),:),'filled')
% ax.XLim(2) = 7;
ylabel('PC2','FontSize',font_sz)
xlabel('PC1','FontSize',font_sz)
thand = title('PCA on Functional Space of All Communities','FontSize',font_sz);
thand.Position(2) = thand.Position(2) + 0.4;
lhand = legend('Low Richness','High Richness','FontSize',font_sz,'Location','northeast');
% lhand.Position(1) = lhand.Position(1); lhand.Position(2) = lhand.Position(2) - 0.775;
%% Variance Explained by Respective PCAs
yhigh_stdrd = (y.highdiversity(:,ppty2incl) - mean(y.highdiversity(:,ppty2incl))) ./ std(y.highdiversity(:,ppty2incl));
ylow_stdrd = (y.lowdiversity(:,ppty2incl) - mean(y.lowdiversity(:,ppty2incl))) ./ std(y.lowdiversity(:,ppty2incl));
[~,~,~,~,explained_high] = pca(yhigh_stdrd,'Centered',false);
[~,~,~,~,explained_low] = pca(ylow_stdrd,'Centered',false);
%%
x0 = x0 + w + dw; y0 = 6.2; w = 3.5; h = 3.5;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
hold on
box on

pl_low = plot(ax,explained_low,'s-','Color',cmap(color_inds(1),:),'MarkerFaceColor',cmap(color_inds(1),:),'LineWidth',lw,'MarkerSize',5);
% pl_mid = plot(explained_mid,'s-','Color',cmap(color_inds(2),:),'MarkerFaceColor',cmap(color_inds(2),:),'LineWidth',2,'MarkerSize',10);
pl_high = plot(ax,explained_high,'s-','Color',cmap(color_inds(3),:),'MarkerFaceColor',cmap(color_inds(3),:),'LineWidth',lw,'MarkerSize',5);
xlim([0.5 4.5])
xticks(1:4)
ylim([0 100])
yticks([0 25 50 75 100])
ylabel('% Variance Explained','FontSize',font_sz)
xlabel('Principal Component','FontSize',font_sz)
title({'High-Richness Subset', 'is Lower Dimensional'},'FontSize',font_sz)
% legend('Low Diversity','High Diversity','FontSize',font_sz-3,'Location','northeast');
%% Paint by pH
y0 = 1.1; w = 4.75; h = 3.5;
ax = axes('Units','centimeters','Position',[x0 y0 w h]);
hold on
box on

projslow = projs(lowdiv_bool,:);
projsmid = projs(middiv_bool,:);
projshigh = projs(highdiv_bool,:);

scatter(ax,projslow(~pHdata.lowdiversity.pH_bool,1),projslow(~pHdata.lowdiversity.pH_bool,2),...
    floor(mkrsz/2),0.8*[1 1 1],'filled');
scatter(ax,projsmid(~pHdata.middiversity.pH_bool,1),projsmid(~pHdata.middiversity.pH_bool,2),...
    floor(mkrsz/2),0.8*[1 1 1],'filled');
scatter(ax,projshigh(~pHdata.highdiversity.pH_bool,1),projshigh(~pHdata.highdiversity.pH_bool,2),...
    floor(mkrsz/2),0.8*[1 1 1],'filled');

scatter(ax,projslow(pHdata.lowdiversity.pH_bool,1),projslow(pHdata.lowdiversity.pH_bool,2),...
    mkrsz,pHdata.lowdiversity.pH,'filled');
scatter(ax,projsmid(pHdata.middiversity.pH_bool,1),projsmid(pHdata.middiversity.pH_bool,2),...
    mkrsz,pHdata.middiversity.pH,'filled');
scatter(ax,projshigh(pHdata.highdiversity.pH_bool,1),projshigh(pHdata.highdiversity.pH_bool,2),...
    mkrsz,pHdata.highdiversity.pH,'filled');
% scatter(ax,proj_on_highdiv_pcs.middiv(pHdata.middiversity.pH_bool,1),proj_on_highdiv_pcs.middiv(pHdata.middiversity.pH_bool,2),...
%     mkrsz,pHdata.middiversity.pH,'filled');
colormap(ax,'parula')
cb = colorbar;
set(cb.Title,'String','pH')
cb.Title.Position = [5 -10 0];
ylabel('PC2','FontSize',font_sz)
xlabel('PC1','FontSize',font_sz)
thand = title('pH Correlates with Dominant PC','FontSize',font_sz);
end