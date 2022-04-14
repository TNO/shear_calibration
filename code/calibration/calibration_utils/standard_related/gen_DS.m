% Generate design situations
%
%SYNOPSYS
% [Prob, DS] = GEN_DS(free_par, Prob, DS, Options)
%
%INPUT
%
%
%OUTPUT
%
%

%NOTES:
% * the for loop could be parallelized
%

function [Prob, DS] = gen_DS(free_par, Prob, Prob_actions, DS, Options)

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
verbose         = Options.verbose;

load_combs      = Options.load_combs;
n_load_comb     = length(load_combs);

d_ds            = DS.Range.d;
f_cck_ds        = DS.Range.f_cck;
chi1_ds         = DS.Range.chi1;
chi2_ds         = DS.Range.chi2;
rho_ds          = DS.Range.rho;
d_lower_ds      = DS.Range.d_lower;
a_to_d_ratio_ds = DS.Range.a_to_d_ratio;

% create all possible combinations for the design scenarios (overkill,
% should be moved outside of the obj_fun)
p_ds_vectors = {chi1_ds, chi2_ds, d_ds, f_cck_ds, rho_ds, d_lower_ds, a_to_d_ratio_ds};
p_ds_2_actions = cartesian(p_ds_vectors{:});

% to avoid looping over chi2 if there is only one variable action
% the indices must be harmonized with the position of chi2
p_ds_1_action   = cartesian(p_ds_vectors{[1, 3:end]});
% insert zero chi2
p_ds_1_action = [
    p_ds_1_action(:,1),...
    zeros(size(p_ds_1_action,1),1),...
    p_ds_1_action(:,2:end)
];

p_ds_all        = [];
weights         = [];
load_combs_all  = cell(n_load_comb*size(p_ds_2_actions,1), 1);

% this should be solved in a more robust way + create list with the variable load set of design scenarios
for ii = 1:n_load_comb
    lc_ii       = load_combs{ii};
    W_b         = DS.Weight.(lc_ii).W;
    Chi1_b      = DS.Weight.(lc_ii).Chi1;
    Chi2_b      = DS.Weight.(lc_ii).Chi2;
    
    if strcmpi(lc_ii, 'traffic')
       p_ds_ii   = p_ds_1_action;
       chi12     = p_ds_ii(:,1:2);
       weights_ii= assign_weight(chi12, Chi1_b, Chi2_b, W_b);
    else
       p_ds_ii   = p_ds_2_actions;
       chi12     = p_ds_ii(:,1:2);
       weights_ii= assign_weight(chi12, Chi1_b, Chi2_b, W_b);
    end
    
    loc1 = size(p_ds_all,1) + 1;
    loc2 = loc1 + size(p_ds_ii,1) - 1;
    load_combs_all(loc1:loc2) = load_combs(ii);
    
    p_ds_all   = [p_ds_all; p_ds_ii];
    weights  = [weights; weights_ii];  
end

idx             = cellfun(@isempty, load_combs_all);
load_combs_all(idx) = []; 

p_ds_all_colnames = {'chi1', 'chi2', 'd', 'f_cck', 'rho', 'd_lower', 'a_to_d_ratio'};
n_ds            = size(p_ds_all, 1);

DS.weights_combis = weights;

Prob_keep       = Prob;

VR = zeros(n_ds, 1);

%--------------------------------------------------------------------------
% Semi-probabilistic design & save the design scenarios
%--------------------------------------------------------------------------
% loop over the design scenarios
% for consistency the mean values are also matched with the representative
% values
t_ds = tic;
for ii = 1:n_ds
    Prob_ii             = Prob_keep;
    
    chi1_ii             = p_ds_all(ii,1);
    chi2_ii             = p_ds_all(ii,2);
    d_repr_ii           = p_ds_all(ii,3);
    f_cck_ii            = p_ds_all(ii,4);
    rho_repr_ii         = p_ds_all(ii,5);
    d_lower_ii          = p_ds_all(ii,6);
    a_to_d_ratio_ii     = p_ds_all(ii,7);
   
    % concrete compressive strength
    Prob_ii.f_cc.repr   = f_cck_ii;
    [~, Prob_ii.f_cc]   = Prob_ii.f_cc.repr2mean(f_cck_ii, Prob_ii.f_cc);

    % effective depth
    Prob_ii.d.repr      = d_repr_ii;
    [~, Prob_ii.d]      = Prob_ii.d.repr2mean(d_repr_ii, Prob_ii.d);

    % rebar area
    b_repr_ii           = Prob_ii.b.repr;
    Asl_repr_ii         = rho_repr_ii * b_repr_ii * d_repr_ii;
    [~, Prob_ii.Asl]    = Prob_ii.Asl.repr2mean(Asl_repr_ii, Prob_ii.Asl);

    % 
    Prob_ii.d_lower.repr = d_lower_ii;
    Prob_ii.d_lower.mean = d_lower_ii;
    Prob_ii.a_to_d_ratio.repr = a_to_d_ratio_ii;
    Prob_ii.a_to_d_ratio.mean = a_to_d_ratio_ii;

    fix_par             = [chi1_ii, chi2_ii];
    
    [Prob_ii, VR_ii]    = inv_design(free_par, fix_par, Prob_ii, Prob_actions, Options, load_combs_all{ii});
    Prob(ii)            = Prob_ii;
    VR(ii)              = VR_ii;
end
t_elapsed_sec = toc(t_ds);

if verbose > 1
    disp('Generation of design scenarios has been completed (semi-probabilistic design).')
    disp(['Excecution time [sec]: ', num2str(t_elapsed_sec)])
    disp('-------------------------------------------------------------------------')
end
    
%--------------------------------------------------------------------------
% Collect results
%--------------------------------------------------------------------------
DS.p_ds_all           = p_ds_all;
DS.p_ds_all_colnames  = p_ds_all_colnames;
DS.t_elapsed_sec      = t_elapsed_sec;
DS.load_combs_all     = load_combs_all;
DS.VR                 = VR;

end