% Example script on applying the method to the E. coli iJO1366 model

%% Get and load the model
if exist('iJO1366.mat', 'file')
    fprintf(' iJO1366.mat found.\n');
else
    % download if not in the search path (need curl)
    [status_curl, result_curl] = system('curl --max-time 15 -s -k -L --head http://bigg.ucsd.edu/static/models/iJO1366.mat');
    % check if the URL exists
    if status_curl == 0 && ~isempty(strfind(result_curl, ' 200'))
        status_curlDownload = system('curl --max-time 60 -O -L http://bigg.ucsd.edu/static/models/iJO1366.mat');
        if status_curlDownload == 0
            fprintf(' + Downloaded:      http://bigg.ucsd.edu/static/models/iJO1366.mat\n');
        end
    else
        fprintf(' > The URL http://bigg.ucsd.edu/static/models/iJO1366.mat cannot be reached.\n');
    end
end

model = readCbModel('iJO1366.mat');
% add biomass metabolites
model = addMetabolite(model, 'biomass_core[e]');
model.S(findMetIDs(model, 'biomass_core[e]'), findRxnIDs(model, 'BIOMASS_Ec_iJO1366_core_53p95M')) = 1;
model = addMetabolite(model, 'biomass_wt[e]');
model.S(findMetIDs(model, 'biomass_wt[e]'), findRxnIDs(model, 'BIOMASS_Ec_iJO1366_WT_53p95M')) = 1;
% add export reactions
model = addReaction(model, {'EX_biomass_core(e)', 'Core biomass export'}, 'biomass_core[e] ->');
model = addReaction(model, {'EX_biomass_wt(e)', 'WT biomass export'}, 'biomass_wt[e] ->');

%% compute chemical formulae
% get all extracellular metabolites
metEx = any(model.S(:, sum(model.S ~= 0, 1) <= 1), 2);
% original chemical formulae
metFormulae = model.metFormulas;
% use all extracellular metabolites to infer all other metabolites
% (the best practice is to use as many known metabolites as one can to
%  reveal all hidden inconsistencies, but all extracellular metabolites +
%  sink/demand metabolites is sufficient to give defined results for all metabolites)
model.metFormulas(~metEx) = {''};
[model2, metCompute, ele, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model);
% compare with the original results
yn = compareMetFormulasLocal(metFormulae, model2.metFormulas);
fprintf('%d / %d metabolites have the same chemical formulae computed.\n', sum(yn), numel(model.mets));
fprintf('Mets with formulae different from the original model:\n');
format short
disp([{'mets', 'Original formulae', 'Computed formulae'}; model.mets(~yn), metFormulae(~yn), model2.metFormulas(~yn)])
% **** It can be seen that the differences all originate from the specific elements used to name
%      each conserved moieties in the network.

%% Look at the inconsistencies found in the results
fprintf('Minimum inconsistencies found when balancing each element:\n')
for j = 1:numel(solInfo.solEachEle)
    fprintf('%s: %.4f\n', strjoin(solInfo.elements(solInfo.eleConnect(:, j)), ', '), solInfo.solEachEle(j).sol.minIncon.obj)
end
% **** should see zero inconsistency for all elements

%% Check molecular weights of the biomass metabolites
bm = findMetIDs(model,{'biomass_core[e]'; 'biomass_wt[e]'});
fprintf('Molecular weight (g/mmol):\nbiomass_core[e]\t%.12f\nbiomass_wt[e]\t%.12f\n', ...
    getMolecularMass(model2.metFormulas{bm(1)}, 0, 1) / 1000, getMolecularMass(model2.metFormulas{bm(2)}, 0, 1) / 1000)
% **** You shall see that the molecular weight of the biomass produced in silico in the E. coli
%      model is pretty accurate.

%% Compute the range for the molecular weight of a metabolite

% Create an inconsistency in the model.

% Take for example the metabolite man6pglyc_c which has only two connecting reactions
% so that the algorithm would deem that it is equally likely for the inconsistency
% to lie in either of the reactions, but not only one of them.
r = findRxnIDs(model, 'MANPGH');
h = findMetIDs(model, 'h_c');
model.S(h, r) = 1;
[metMWrange, metForm, ele, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model, 'metMwRange', 'man6pglyc_c');
model.S(h, r) = 0;
% Look at the inconsistencies found in the results
fprintf('Minimum inconsistencies found when balancing each element:\n')
for j = 1:numel(solInfo.solEachEle)
    fprintf('%s: %.4f\n', strjoin(solInfo.elements(solInfo.eleConnect(:, j)), ', '), solInfo.solEachEle(j).sol.minIncon.obj)
end
fprintf('man6pglyc_c:\n')
fprintf('min MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrange(1), metForm{1}, strjoin(model.rxns(any(abs(rxnBal.minMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('max MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrange(2), metForm{2}, strjoin(model.rxns(any(abs(rxnBal.maxMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
% Note that the the difference between metMWrange(1) and metMWrange(2) is
% 1.0XXX. 1 due to the difference in the hydrogen atom, .0XXX due to the
% tolerance allowed when optimizing for each element.
% And the unbalanced reaction in the min case is different from that in the max case.

% But if we take a metabolite, e.g. glc__D_c (glucose) or fru_c (fructose) with many connecting
% reactions, the reaction with the one single inconsistency will be pinned down
r = findRxnIDs(model, 'XYLI2');
h = findMetIDs(model, 'h_c');
model.S(h, r) = 1;
[metMWrangeGlc,metFormGlc,ele,metEle,rxnBalGlc, ~, solGlc] = computeMetFormulae(model, 'metMwRange', 'glc__D_c');
[metMWrangeFru,metFormFru,~,~,rxnBalFru, ~, solFru] = computeMetFormulae(model, 'metMwRange', 'fru_c');
model.S(h, r) = 0;
fprintf('Minimum inconsistencies found:\nelement\tglc__D_c\tfru_c\n')
for j = 1:numel(solInfo.solEachEle)
    fprintf('%s:\t%.4f\t%.4f\n', strjoin(solInfo.elements(solInfo.eleConnect(:, j)), ', '), ...
        solGlc.solEachEle(j).sol.minIncon.obj, solFru.solEachEle(j).sol.minIncon.obj)
end

fprintf('\nglc__D_c:\n')
fprintf('min MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeGlc(1), metFormGlc{1}, strjoin(model.rxns(any(abs(rxnBalGlc.minMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('max MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeGlc(2), metFormGlc{2}, strjoin(model.rxns(any(abs(rxnBalGlc.maxMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('\nfru_c:\n')
fprintf('min MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeFru(1), metFormFru{1}, strjoin(model.rxns(any(abs(rxnBalFru.minMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
fprintf('max MW case:\n   MW: %.4f\n   Formula: %s\n   Unbalanced reaction: %s\n', ...
    metMWrangeFru(2), metFormFru{2}, strjoin(model.rxns(any(abs(rxnBalFru.maxMW)>1e-5 & sum(model.S ~= 0, 1) > 1, 1)), ', '))
% Note that here the differences between max and min MW are smaller and << 1,
% only due to the tolerance allowed. And the unbalanced reaction is unique.

% One can use the above method to check the MW of the biomass produced in
% silico regardless of the existing inconsistency

function [yn,element,metEle] = compareMetFormulasLocal(metFormulas_1, metFormulas_2)
% Given two cell arrays of the same length of chemical formulas 'form_1' and
% 'form_2', determine if they are the same.
%
% USAGE:
%    [yn,ele,metEle] = compareMetFormulas(metFormulas_1, metFormulas_2)
%
% INPUTS:
%    metFormulas_1:   a cell array of chemical formulae
%    metFormulas_2:   a second cell array of chemical formulae for comparison, same length as metFormulas_1
%
% OUTPUTS:
%    yn:        logical vector, true if equal, false if unequal
%    ele:       cell array of elements in the formulas
%    metEle:    numel(form1) x numel(ele) x 2 matrix of metFormulas_1 and metFormulas_2
%               metEle(i,j,k) is the stoichiometry of element{j} in form_i

if iscell(metFormulas_1) && iscell(metFormulas_2) && numel(metFormulas_1) ~= numel(metFormulas_2)
    error('The number of chemical formulas in each cell is different!')
end
[metEle1, ele1] = getElementalComposition(metFormulas_1);
[metEle2, ele2] = getElementalComposition(metFormulas_2);
if iscell(metFormulas_1) && ischar(metFormulas_2)
    metEle2 = repmat(metEle2,numel(metFormulas_1),1);
elseif ischar(metFormulas_1) && iscell(metFormulas_2)
    metEle1 = repmat(metEle1,numel(metFormulas_2),1);
end
nForm = size(metEle1,1);



element = union(ele1, ele2);
[~,id] = ismember({'C';'H';'N';'O';'P';'S'}, element);
id2 = setdiff(1:numel(element), id);
element = element([id(id~=0); id2(:)]);
id = strcmp(element,'Charge');
if any(id)
    element = [element(~id); element(id)];
end

metEle = zeros(nForm,numel(element), 2);
[~,idEl1] = ismember(ele1,element);
[~,idEl2] = ismember(ele2,element);
metEle(:,idEl1,1) = metEle1;
metEle(:,idEl2,2) = metEle2;

yn = all(abs(metEle(:,:,1) - metEle(:,:,2)) < 1e-6, 2);
end
