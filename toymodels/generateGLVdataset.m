function dataset = generateGLVdataset(GLV, divRange, sampleN, focal_strain_ID, survivors, params)
% choose subsets
sampleDiversity = randi(divRange, [1, sampleN]);
present = false(sampleN, GLV.S);
pop2sample = survivors;
pop2sample(pop2sample == focal_strain_ID) = [];
for ii=1:sampleN
    pick = randsample(pop2sample, sampleDiversity(ii)); % assemble communities without focal strain
    present(ii, pick) = true;
end
dataset.present = present;
abd_preinv = NaN(size(present));
% get abundance for each community
% if a divergence is ever encountered, stop.
for ii=1:sampleN
    thisAbd = GLV.getAbd(present(ii,:));
    if any(isnan(thisAbd))
        fprintf('Divergence detected\n');
        break;
    else
        abd_preinv(ii,:) = thisAbd;
    end
end
dataset.abd_preinv = abd_preinv;
% now invade with pathogen
present(:,focal_strain_ID) = true;
abd_postinv = NaN(size(present));
for ii=1:sampleN
    params.abd0(pop2sample) = abd_preinv(ii,pop2sample);
    GLV_invade = V2cflLotkaVolterra(params);
    thisAbd = GLV_invade.getAbd(present(ii,:));
    if any(isnan(thisAbd))
        fprintf('Divergence detected\n');
        break;
    else
        abd_postinv(ii,:) = thisAbd;
    end
end
dataset.abd_postinv = abd_postinv;
end