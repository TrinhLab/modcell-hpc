% Test that deletions and associated objective values calculated by
% modcell-hpc are consistent with modcell2
function test_objectives(solution_path, prodnet_path, design_objective)
%
% Inputs:
%   - solution_path: .csv from modcell-hpc
%   - prodnet_path: path to prodnet matlab structure associated with the problem.
%

tol = 0.0001;

if ~exist('design_objective', 'var')
	design_objective = 'wGCP';
end
pn = getfield(load(prodnet_path), 'prodnet');

T_in = readtable(solution_path);
T_calc = calc_objs(solution_path, pn, design_objective);

%fprintf('Sol-idx \t Obj-id \t objval-read \t objval-calc.\n')
inconsitency_found = 0;
for i = 1:length(T_in.SolutionIndex)
	for k = 1:length(pn.prod_id)
		T_in_idx = contains(T_in.Properties.VariableNames, [pn.prod_id{k},'_objective_']);
		T_calc_idx = contains(T_calc.Properties.VariableNames, [pn.prod_id{k},'_objective_calc']);
		if (abs(T_in{i, T_in_idx} - T_calc{i, T_calc_idx}) > tol)
			fprintf("Sol idx:%d \t Obj:%s \t In obj.: %2.2f \t Calc. obj.: %2.2f\t In-Out: %2.2f\n", ...
			T_in.SolutionIndex(i), pn.prod_id{k}, T_in{i, T_in_idx}, T_calc{i, T_calc_idx}, T_in{i, T_in_idx} - T_calc{i, T_calc_idx});
			inconsitency_found = 1;
		end
	end
end
if(~inconsitency_found)
	fprintf("No issues found!\n\n")
end
end


function [T_out, PF, design_vars] = calc_objs(output_file_path, prodnet, design_objective)
% Calculates the pareto front of the deletions and module reactions
% Based on format_output (for MIP in ModCell2)

prodnet.set_deletion_type('reactions', design_objective);

T = readtable(output_file_path);
PF = zeros(length(T.Deletion_id),length(prodnet.prod_id));
design_vars = [];
for i = 1:length(T.SolutionIndex)

	%% Parse deletions
	deletions = parse_list(T.Deletion_id{i});

	y = false(1,length(prodnet.cand_ind));
	[~,model_del_ind] = intersect(prodnet.parent_model.rxns, deletions,'stable');
	[~,cand_del_ind] = intersect(prodnet.cand_ind,model_del_ind,'stable');
	y(cand_del_ind) = true;

	%% Parse module reactions
	Z = false(length(prodnet.prod_id), length(prodnet.cand_ind));
	var_module_id = T.Properties.VariableNames(contains(T.Properties.VariableNames, '_module'));
	module_prod_id = cellfun(@(x)(strrep(x,'_module_','')),var_module_id,'UniformOutput',false);
	for k = 1:length(prodnet.prod_id)
		if ~isempty(intersect(module_prod_id, prodnet.prod_id{k}))
            modules = T{i,[prodnet.prod_id{k},'_module_']};
            if ~isempty(modules) || isnan(modules)
                module_rxn = {};
            else
                module_rxn = parse_list(T{i,[prodnet.prod_id{k},'_module_']}{1});
            end
			[~,model_del_ind] = intersect(prodnet.parent_model.rxns, module_rxn,'stable');
			[~,cand_del_ind] = intersect(prodnet.cand_ind,model_del_ind,'stable');
			Z(k, cand_del_ind) = true;
	end
end

%% Compute Pareto front
prodnet.set_module_and_deleted_variables(Z, y)
PF(i,:) = prodnet.calc_design_objectives(design_objective);
design_vars(i).Z = Z;
design_vars(i).y = y;
end
out_prod_id = safe_prod_id(prodnet.prod_id);
T_PF = array2table(PF,...
'VariableNames',cellfun(@(x)([x,'_objective_calc']),out_prod_id, 'UniformOutput',false));
T_out = [T, T_PF];

end

function cell_out = parse_list(list_str)

if isempty(list_str)
	cell_out = {};
else
	cell_out = textscan(list_str,'%s','Delimiter',',');
	cell_out = cell_out{:};
end
end

function out_prod_id = safe_prod_id(prod_id)
% Matlab tables can deal with numbers in ids..
out_prod_id = prod_id;
for k=1:length(prod_id)
	%out_prod_id{k} = strrep(out_prod_id{k}, '_', 'U');
    if regexp(prod_id{k}, '^\d')
		out_prod_id{k} = ['z', prod_id{k}];
    end
end
end
