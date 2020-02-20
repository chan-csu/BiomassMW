%Need COBRA toolbox

%The example model is the intermediate Bacteroides thetaiotaomicron iAH991 model after
%correcting the reactions (except proton imbalance, which can be identified 
%by the procedure) but before renormalizing the biomass reaction.
model = readCbModel('example2_model.mat');
%Compute the chemical formulae for unknown metabolites
[model2, metCompute, ele, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model);
%IDs for unknown mtabolites
metU = findMetIDs(model, solInfo.metUnknown);
%IDs for known metabolites
metK = setdiff(1:numel(model.mets), metU);
%Generic conserved moieties
CMgen = find(~any(solInfo.N(metK,:), 1));
for j = 1:numel(CMgen)
    fprintf('Generic conserved moiety %d:\n', j);
    m = find(solInfo.N(:,CMgen(j)) ~= 0);
    for k = 1:numel(m)
        fprintf('\t%s\t%s\t%s\n', model2.mets{m(k)}, model2.metNames{m(k)}, model2.metFormulas{m(k)});
    end
end
bmId = findMetIDs(model,'biomass[e]');
eleCharge = strcmp(ele,'Charge');
bmForm = elementalMatrixToFormulae(metEle(bmId,~eleCharge), ele(~eleCharge));
fprintf('Biomass formula: %s\nCharge: %.6f\n', bmForm{:}, full(metEle(bmId,eleCharge)));
fprintf('\nFind the range for biomass MW (g/mmol):\n');
[metMWrange,metForm] = computeMetFormulae(model, 'metMwRange', 'biomass[e]');
fprintf('[%.6f, %.6f]\n', metMWrange(:) / 1000);

% unbalanced reactions
rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & sum(model.S ~= 0, 1) > 1);
fprintf('\nCurrently %d unbalanced reactions:\n%s\n', numel(rxnImbal), strjoin(rxnImbal,'\n'));

% fill the predicted proton imbalance from 'S_fill'
h = findMetIDs(model, {'h[c]'; 'h[e]'});
rxnFill = find(S_fill);
for j = 1:numel(rxnFill)
    metJ = model.mets(model.S(:, rxnFill(j)) ~= 0);
    comp = unique(regexp(metJ, '\[\w*\]$', 'match', 'once'));
    if numel(comp) == 1
        if strcmp(comp, '[c]')
            model.S(h(1), rxnFill(j)) = model.S(h(1), rxnFill(j)) + S_fill(rxnFill(j));
        elseif strcmp(comp, '[e]')
            model.S(h(2), rxnFill(j)) = model.S(h(2), rxnFill(j)) + S_fill(rxnFill(j));
        end
    else
        surfNet(model, model.rxns(rxnFill(j)))
        fprintf('Which proton should be used to fill the imbalance?\n1:  h[c]\n2:  h[e]\n')
        s = input('Enter [1] or [2]:  ', 's');
        while str2double(s) ~= 1 && str2double(s) ~= 2
            s = input('Enter [1] or [2]:  ', 's');
        end
        model.S(h(str2double(s)), rxnFill(j)) = model.S(h(str2double(s)), rxnFill(j)) + S_fill(rxnFill(j));
    end
end
fprintf('Reactions after filling proton imbalance:\n')
surfNet(model, model.rxns(rxnFill), 's', 0)

[model2, metCompute, ele, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model);
rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & sum(model.S ~= 0, 1) > 1);
fprintf('\n%d unbalanced reactions after proton filling:\n%s\n', numel(rxnImbal), strjoin(rxnImbal,'\n'));

% the only unbalanced reaction is AHMMPS, which is a common gapfilling
% reaction but with charge imbalance
