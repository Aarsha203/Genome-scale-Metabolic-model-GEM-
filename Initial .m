% Find the index of the 'EX_glc(e)' reaction in the model
rxnIndex = find(ismember(model.rxns, 'EX_glc(e)'));

% Set the new upper bound value
newUpperBoundValue = -0.1;

% Update the upper bound of the reaction with the new value
model.ub(rxnIndex) = newUpperBoundValue;
%%%%% % Find the index of the reaction ADK1 in model.rxns
reactionIndex = find(strcmp(model.rxns, 'CO2tm'));
% Check if the reaction ADK1 exists in the model
if ~isempty(reactionIndex)
% Get the gene associated with the reaction ADK1 from model.gene
associatedGene = model.genes(reactionIndex);
% Print the associated gene
fprintf('The gene associated with the reaction ADK1 is: %s\n', associatedGene{:});
else
fprintf('The reaction ADK1 does not exist in the model.\n');
end
