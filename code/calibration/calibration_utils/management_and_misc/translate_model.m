% Translate between character and number based model definitions
%
% THe numerical representation is needed due to how FERUM expects limit
% state functions to be defined (only numerical arguments)

function tmodel = translate_model(model)

% -------------------------------------------------------------------------
% Dictionary
% -------------------------------------------------------------------------
% Resistance model
% Resistance model
resi_model_str = {'ec2_codified_2019', 'ec2_pre_2021', 'mc2010_level_ii_codified_2019'};
resi_model_num = [1, 2, 3];

% Load combination
load_comb_str = {'ec0_simple', 'ec0_advanced'};
load_comb_num = [101          , 102];

model_str = [resi_model_str, load_comb_str];
model_num = [resi_model_num, load_comb_num];

% -------------------------------------------------------------------------
% Translate
% -------------------------------------------------------------------------
if ischar(model)
    idx     = strcmpi(model, model_str);
    if any(idx)
        tmodel  = model_num(idx);
    else
        error(['Unknown string-defined model:', model])
    end
elseif isnumeric(model)
    idx     = model == model_num;
    if any(idx)
        tmodel  = model_str{idx};
    else
        error(['Unknown number-defined model:', model])
    end
else
    error('Unknown model variable type!')
end

end