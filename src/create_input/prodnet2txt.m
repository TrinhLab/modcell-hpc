function prodnet2txt(prodnet_path, design_objective, output_dir)

% Converts key model information in prodnet datastructure to plain text files. Including a file with candidates for deletion and a file for each production network in MPS format.

% Notes:
% - The appropriate objective, based on the overall design objective (e.g. wGCP) is set for the models.
% - MPS assumes minimization so all objectives are negated since ModCell assumes maximization.

% Todo:
% - Make a map between mps safe ids and actual ids.
% - Write txt of values used to normalize objectives. It must follow the format of prod_id, value (since the production network parsing order is unreliable)

%% parse inputs
if ~exist('design_objective', 'var')
	design_objective = 'wGCP';
end

if ~exist('output_dir', 'var')
	output_dir = '.';
end

switch design_objective
case {'wGCP', 'sGCP', 'lsGCP'}
	state='growth';
case {'NGP'}
	state='nongrowth';
end

pn = getfield(load(prodnet_path), 'prodnet');

%% Set state
pn.set_mip_state(design_objective, state, 0, 1); % Why is growth specified?
pn.set_deletion_type('reactions', design_objective);

%% Make IDs safe for MPS format
% Currently the only manipulation is truncating the strings to a max name
all_ids = {};
for i =1:length(pn.prod_id)
    all_ids = union(all_ids, pn.model_array(i).rxns);
end
max_length = 8;
new_ids = cellfun(@(x)(truncatestr(x, max_length)),all_ids, 'UniformOutput',false);
if length(new_ids) ~= length(unique(new_ids))
    error('Reaction ID truncation led to duplicates, add some salt to he truncated ids to avoid this')
end

idmap = containers.Map(all_ids, new_ids);
for i =1:length(pn.prod_id)
    for j =1:length(pn.model_array(i).rxns)
        pn.model_array(i).rxns{j} = idmap(pn.model_array(i).rxns{j});
    end
end

writetable(table(all_ids, new_ids), fullfile(output_dir, 'rxnidmap.csv'))


%% Write candidates
fileID = fopen(fullfile(output_dir, 'cand'),'w');
for i = 1:length(pn.cand_ind)
	fprintf(fileID,'%s\n', pn.parent_model.rxns{pn.cand_ind(i)});
end
fclose(fileID);

%% Write models and model related info
fprintf('Writting model info.... \n')
for i =1:length(pn.prod_id)
	fprintf('\t%s\n', pn.prod_id{i})
	model = pn.model_array(i);

	[Contain OK] = mcBuildMPS([], [], model.S, model.b, -model.c, model.lb, model.ub, pn.prod_id{i}, 'VarNames', model.rxns, 'MPSfilename', [pn.prod_id{i}, '.mps']);
	if ~OK
		warning('Something went wrong when building mps file for model %s', pn.prod_id{i})
	end

	% Write non candidates
	fileID = fopen([pn.prod_id{i}, '.ncand'],'w');
	ncand = setdiff(pn.parent_model.rxns(pn.cand_ind), model.rxns(model.candk));
	for j = 1:length(ncand)
		fprintf(fileID,'%s\n', ncand{j});
	end
	fclose(fileID);

	% Write parameters
	fileID = fopen([pn.prod_id{i}, '.param'],'w');
	fprintf(fileID,'prod_rxn_name=%s\n', model.rxns{model.product_secretion_ind});
	fprintf(fileID,'max_prod_growth=%f\n', pn.max_product_rate_growth(i));
	fclose(fileID);
end
end

function newstr = truncatestr(instr, max_length)
if length(instr) < max_length
    newstr = instr;
else
    newstr = instr(1:max_length);
end
end
