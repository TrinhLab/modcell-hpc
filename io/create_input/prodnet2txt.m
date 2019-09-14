function prodnet2txt(prodnet_path, design_objective, output_dir)

% Converts key model information in prodnet datastructure to plain text files. Including a file with candidates for deletion and a file for each production network in MPS format.

% Notes:
% - The appropriate objective, based on the overall design objective (e.g. wGCP) is set for the models.
% - MPS assumes minimization so all objectives are negated since ModCell assumes maximization.


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

%% Comply to MPS format absurd restrictions

%% 1)Change problem ids to be a maximum length of 7 characters
og_model_ids = pn.prod_id;
new_model_ids = {};
for i=1:length(pn.prod_id)
	new_model_ids{i,1} = truncatestr_nosalt(og_model_ids{i}, 7);
end
if length(og_model_ids) ~= length(unique(new_model_ids))
    error('Truncated model ids non-unique, add some salt')
end
writetable(table(og_model_ids, new_model_ids), fullfile(output_dir, 'modelidmap.csv'))
pn.prod_id = new_model_ids;

%% 2)Make IDs safe for MPS format
% - Max length is truncated to 8
% - An integer between 10 and 99 is added at the end to avoid repetition.
% - No other checks are performed.(e.g. illegal characters)

% Generate new ids
all_ids = {};
for i =1:length(pn.prod_id)
    all_ids = union(all_ids, pn.model_array(i).rxns);
end
max_length = 8;
global salt_counter;
salt_counter = 10;
new_ids = cellfun(@(x)(truncatestr(x, max_length)),all_ids, 'UniformOutput',false);
if length(new_ids) ~= length(unique(new_ids))
    error('Reaction ID truncation led to duplicates, add more salt to he truncated ids to avoid this')
end

% Update model reactions and other places where ids appear
idmap = containers.Map(all_ids, new_ids);

	% Parent model
for j =1:length(pn.parent_model.rxns)
	pn.parent_model.rxns{j} = idmap(pn.parent_model.rxns{j});
end
	% Production networks
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
% Modifying the end of the string keeps reaction names readable.

global salt_counter;
if salt_counter > 99
	salt_counter = 10;
end

newstr = instr;
if length(instr) >= max_length
    newstr = newstr(1:max_length);
    newstr(end-1:end) = num2str(salt_counter);
elseif length(instr)  == (max_length - 1)
    newstr(end:end+1) = num2str(salt_counter);
else
    newstr(end+1:end+2) = num2str(salt_counter);
end

salt_counter = salt_counter+1;
end


function newstr = truncatestr_nosalt(instr, max_length)
%
newstr = instr;
if length(instr) >= max_length
    newstr = newstr(1:max_length);
end
end
