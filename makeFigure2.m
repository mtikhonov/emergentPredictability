function makeFigure2(num_invaders)
generateFigure2Data(num_invaders);
plotFigure(num_invaders);
end

function generateFigure2Data(invadersN)
addpath toymodels/
addpath toymodels/GLV/
addpath toymodels/ecoevoCRM/
%% Run GLV Model Experiments
params.S = 45; % number of strains
params.Teq = 1e4;

% Match params from Fig. S2 of Barbier et al, PNAS
std_K = 0.3;
abd_std = 0.1;
N = params.S;

rng(0,'combRecursive');
params.growth = ones([N,1]);
params.abd0   = exprnd(abd_std, [N,1]);
% params.abd0(end) = 0; % pathogen abd is 0
params.propertyOfInterest = @(abd) abd; % will get pathogen abd when finished

mu = 1;
sigma = 0.8;
A = normrnd(mu/N, sigma/sqrt(N), [N,N]);
A = tril(A,-1);
params.inters = A+A';

while true % get random K but ensure they are all positive
    params.K = 1 + std_K*randn([N,1]);
    if all(params.K>0), break, end
end

GLV_full = V2cflLotkaVolterra(params);
abd_full = GLV_full.getAbd(true(params.S,1));
abd_full(abd_full < 1e-6) = 0;
survivors_inds_full = find(abd_full > 0);
%% Invade communities of varying richness
for i = 1:invadersN
    if ~exist(['figureData/GLV_invader_' num2str(i) '.mat'],'file')
        % draw random invader
        invParams = params;
        invParams.growth = [invParams.growth; 1];
        invAbd0 = exprnd(abd_std,1);
        invParams.abd0 = [invParams.abd0; invAbd0];
        invA = normrnd(mu/N, sigma/sqrt(N), [1,N+1]);
        invA(end) = 0;
        invParams.inters = [invParams.inters; invA(1:end-1)];
        invParams.inters = [invParams.inters invA'];
        invK = 1 + std_K*randn(1);
        while invK <= 0
            invK = 1 + std_K*randn(1);
        end

        invParams.K = [invParams.K; invK];
        GLV_full_plusinv = V2cflLotkaVolterra(invParams);
        % invade
        fprintf('iteration %d\n',i)
        focal_ID = N+1;
        fprintf('Low-diversity samples...\n')
        divRange = [1, 5];
        sampleN = 1000;
        data.lowdiv = generateGLVdataset(GLV_full_plusinv, divRange, sampleN, focal_ID, survivors_inds_full, invParams);

        fprintf('Mid-diversity samples...\n')
        divRange = [11, 15];
        sampleN = 1000;
        data.middiv = generateGLVdataset(GLV_full_plusinv, divRange, sampleN, focal_ID, survivors_inds_full, invParams);

        fprintf('High-diversity samples...\n')
        divRange = [21, 25];
        sampleN = 1000;
        data.highdiv = generateGLVdataset(GLV_full_plusinv, divRange, sampleN, focal_ID, survivors_inds_full, invParams);

        data.survivors_inds_full = survivors_inds_full;
        data.focal_strain_ind = focal_ID;

        fprintf('Saving model data...\n');
        save(['figureData/GLV_invader_' num2str(i) '.mat'],'data','params')
    end
    % make Pareto front
    fprintf('Collecting Pareto Front for GLV invader %d ...\n',i)
    generateGLVparetofront(i)
end
%% Run CR Model Experiments
clearvars params
% define parameters
params.invadersN = invadersN;
params.L = 40; % num traits each strain has (sets max num of possible coexisting strains)
params.M = 50; % num random strains to generate for finding coexisting community
params.p = 0.5; % density of traits in strains

params.K0 = 1e10;
params.dK = 0.1;
params.rho = 1;
params.T2h = @(T,env) env.rho ./ (1+max(1e-10,T)./env.K);

params.c = 0.1;
params.lambda = 0.5;
params.dchi = 0.01;

params.seed = 1;
params.nTrials = 1000;
params.simTime = 1e7;
% simulate dataset
rng(params.seed,'combRecursive');
generateCRMdataset(params);
% make Pareto front
for i = 1:params.invadersN
    fprintf('Collecting Pareto Front for CRM invader %d...\n',i)
    generateCRMparetofront(i)
end
end

function plotFigure(invadersN)
%% Panel A: average GGLV Pareto Fronts
clrs = [197/256, 90/256, 17/256;
    175/256, 170/256 170/256;
    143/256, 69/256, 199/256;];
cmap = myColorMap(clrs,150);
color_inds = floor(linspace(1,150,3));
clf;
font_sz = 11;
W = 16.51; % two column figure
H = 5.75;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[6 6 W H]);
x0 = 1.25; y0 = 1.1; w = 4; h = 3.5; dw = 1.1;
axA = axes('Units','centimeters','Position',[x0 y0 w h]);

% load data
for i = invadersN:-1:1
    load(['figureData\GLV_invader_' num2str(i) '_pareto_data.mat'],'GLV')
    desc_compl{i,1} = GLV.lowdiv.desc_complexity;
    desc_compl{i,2} = GLV.middiv.desc_complexity;
    desc_compl{i,3} = GLV.highdiv.desc_complexity;
    median_FVUs{i,1} = GLV.lowdiv.FVUmedian;
    median_FVUs{i,2} = GLV.middiv.FVUmedian;
    median_FVUs{i,3} = GLV.highdiv.FVUmedian;
end
fvuAvg = cell(1,3);
fvuStd = cell(1,3);
all_x = cell(1,3);
for i = 1:3 % loop over richness
    % find all values on complexity axis (x-axis)
    all_x{i} = unique([desc_compl{:,i}]);
    % interpolate each Pareto front at each of these x-values
    interp_fvu = nan(invadersN,length(all_x{i}));
    for j = 1:invadersN
        interp_fvu(j,:) = interp1(desc_compl{j,i},median_FVUs{j,i},all_x{i});
    end
    % average the actual and interpolated values
    fvuAvg{i} = median(interp_fvu,'omitmissing');
    fvuStd{i} = std(interp_fvu,'omitmissing');
end
% plot
x2plt = all_x{2}; y2plt = fvuAvg{2}; err2plt = fvuStd{2};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axA,x2plt,y2plt,cmap(color_inds(2),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{2}; y2plt = fvuAvg{2};
phmid = plot(axA,x2plt,y2plt,'-', 'Color',cmap(color_inds(2),:),'linewidth',1.5);

x2plt = all_x{1}; y2plt = fvuAvg{1}; err2plt = fvuStd{1};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axA,x2plt,y2plt,cmap(color_inds(1),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{1}; y2plt = fvuAvg{1};
plow = plot(axA,x2plt,y2plt,'-', 'Color',cmap(color_inds(1),:),'linewidth',1.5);

x2plt = all_x{3}; y2plt = fvuAvg{3}; err2plt = fvuStd{3};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axA,x2plt,y2plt,cmap(color_inds(3),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{3}; y2plt = fvuAvg{3};
phigh = plot(axA,x2plt,y2plt,'-', 'Color',cmap(color_inds(3),:),'linewidth',1.5);
box on
ylim([0 1.05])
% axis square
xlh = xlabel(axA,'Description Entropy');
xlh.Position(1) = xlh.Position(1) + 2.5;
ylabel('FVU')
set(gca,'FontSize',font_sz-1.25)
title({'Post-invasion Abundance', 'in GLV'},'FontSize',font_sz)
%% Panel B: average CRM Pareto Fronts
x0 = x0 + w + dw;
axB = axes('Units','centimeters','Position',[x0 y0 w h]);
% load data
for i = invadersN:-1:1
    load(['figureData/CRM_invader_' num2str(i) '_pareto_data.mat'],'CRM')
    desc_compl{i,1} = CRM.lowdiv.desc_complexity;
    desc_compl{i,2} = CRM.middiv.desc_complexity;
    desc_compl{i,3} = CRM.highdiv.desc_complexity;
    median_FVUs{i,1} = CRM.lowdiv.FVUmedian;
    median_FVUs{i,2} = CRM.middiv.FVUmedian;
    median_FVUs{i,3} = CRM.highdiv.FVUmedian;
end
fvuAvg = cell(1,3);
fvuStd = cell(1,3);
all_x = cell(1,3);
for i = 1:3 % loop over richness
    % find all values on complexity axis (x-axis)
    all_x{i} = unique([desc_compl{:,i}]);
    % interpolate each Pareto front at each of these x-values
    interp_fvu = nan(invadersN,length(all_x{i}));
    for j = 1:invadersN
        interp_fvu(j,:) = interp1(desc_compl{j,i},median_FVUs{j,i},all_x{i});
    end
    % average the actual and interpolated values
    fvuAvg{i} = median(interp_fvu,'omitmissing');
    fvuStd{i} = std(interp_fvu,'omitmissing');
end
% plot
x2plt = all_x{2}; y2plt = fvuAvg{2}; err2plt = fvuStd{2};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axB,x2plt,y2plt,cmap(color_inds(2),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{2}; y2plt = fvuAvg{2};
phmid = plot(axB,x2plt,y2plt,'-', 'Color',cmap(color_inds(2),:),'linewidth',1.5);

x2plt = all_x{1}; y2plt = fvuAvg{1}; err2plt = fvuStd{1};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axB,x2plt,y2plt,cmap(color_inds(1),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{1}; y2plt = fvuAvg{1};
plow = plot(axB,x2plt,y2plt,'-', 'Color',cmap(color_inds(1),:),'linewidth',1.5);

x2plt = all_x{3}; y2plt = fvuAvg{3}; err2plt = fvuStd{3};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axB,x2plt,y2plt,cmap(color_inds(3),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{3}; y2plt = fvuAvg{3};
phigh = plot(axB,x2plt,y2plt,'-', 'Color',cmap(color_inds(3),:),'linewidth',1.5);
box on
% axis square
% xlabel('Description Complexity')
% ylabel('FVU')
ylim([0 1.05])
xlim(axA.XLim)
set(gca,'FontSize',font_sz-1.25)
title({'Post-invasion Abundance', 'in CRM'},'FontSize',font_sz)
%% Panel C: average GGLV again, but with absolute variance
x0 = x0 + w + dw + 0.45;
axC = axes('Units','centimeters','Position',[x0 y0 w h]);
% load data
for i = invadersN:-1:1
    load(['figureData/CRM_invader_' num2str(i) '_pareto_data.mat'],'CRM')
    desc_compl{i,1} = CRM.lowdiv.desc_complexity;
    desc_compl{i,2} = CRM.middiv.desc_complexity;
    desc_compl{i,3} = CRM.highdiv.desc_complexity;
    median_AVUs{i,1} = CRM.lowdiv.AVUmedian;
    median_AVUs{i,2} = CRM.middiv.AVUmedian;
    median_AVUs{i,3} = CRM.highdiv.AVUmedian;
end
avuAvg = cell(1,3);
avuStd = cell(1,3);
all_x = cell(1,3);
for i = 1:3 % loop over richness
    % find all values on complexity axis (x-axis)
    all_x{i} = unique([desc_compl{:,i}]);
    % interpolate each Pareto front at each of these x-values
    interp_avu = nan(invadersN,length(all_x{i}));
    for j = 1:invadersN
        interp_avu(j,:) = interp1(desc_compl{j,i},median_AVUs{j,i},all_x{i});
    end
    % average the actual and interpolated values
    avuAvg{i} = median(interp_avu,'omitmissing');
    avuStd{i} = std(interp_avu,'omitmissing');
end
% plot
x2plt = all_x{2}; y2plt = avuAvg{2}; err2plt = avuStd{2};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axC,x2plt,y2plt,cmap(color_inds(2),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{2}; y2plt = avuAvg{2};
phmid = plot(axC,x2plt,y2plt,'-', 'Color',cmap(color_inds(2),:),'linewidth',1.5);

x2plt = all_x{1}; y2plt = avuAvg{1}; err2plt = avuStd{1};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axC,x2plt,y2plt,cmap(color_inds(1),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{1}; y2plt = avuAvg{1};
plow = plot(axC,x2plt,y2plt,'-', 'Color',cmap(color_inds(1),:),'linewidth',1.5);

x2plt = all_x{3}; y2plt = avuAvg{3}; err2plt = avuStd{3};
x2plt = [x2plt fliplr(x2plt)];
y2plt = [y2plt-err2plt fliplr(y2plt+err2plt)];
fill(axC,x2plt,y2plt,cmap(color_inds(3),:),'FaceAlpha',0.5,'EdgeColor','none')
hold on
x2plt = all_x{3}; y2plt = avuAvg{3};
phigh = plot(axC,x2plt,y2plt,'-', 'Color',cmap(color_inds(3),:),'linewidth',1.5);
box on
% axis squarer
xlabel('Description Complexity')
ylabel('AVU')
xlim(axB.XLim)
yticks([0 1e18 2e18])
yticklabels({'0','1 x 10^{18}', '2 x 10^{18}'})
% set(gca,'YScale','log')
% axC.YAxis.Exponent = 0;
ytickangle(90)
set(gca,'FontSize',font_sz-1.25)
title({'Post-invasion Abundance', 'in CRM (AVU)'},'FontSize',font_sz)
lhand = legend(gca,[plow phmid phigh],{'Low','Mid','High'},'Location','northeast','FontSize',font_sz-2);
lhand.Position(2) = lhand.Position(2) - 0.04;
text(lhand.Position(1)+1.8,1.9e18,lhand.Position(3),'Richness','FontSize',font_sz-1)
end