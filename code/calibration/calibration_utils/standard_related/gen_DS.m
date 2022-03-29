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
n_lc            = length(load_combs);

d_ds            = DS.Range.d;
f_cck_ds        = DS.Range.f_cck;
chi1_ds         = DS.Range.chi1;
chi2_ds         = DS.Range.chi2;
rho_ds          = DS.Range.rho;
d_lower_ds      = DS.Range.d_lower;
a_to_d_ratio_ds = DS.Range.a_to_d_ratio;

% create all possible combinations for the design scenarios (overkill,
% should be moved outside of the obj_fun)
combis1         = combvec(d_ds, f_cck_ds, chi1_ds, chi2_ds, rho_ds, d_lower_ds, a_to_d_ratio_ds)';
[~, ia]         = unique(combis1(:,[1,2,3,5,6,7]), 'rows');
combis2         = combis1(ia,:);
combis2(:,4)    = 0;
combis          = [];
weights         = [];
load_combs_all  = cell(n_lc*size(combis1,1), 1);

% this should be solved in a more robust way + create list with the variable load set of design scenarios
for ii = 1:n_lc
    lc_ii       = load_combs{ii};
    W_b         = DS.Weight.(lc_ii).W;
    Chi1_b      = DS.Weight.(lc_ii).Chi1;
    Chi2_b      = DS.Weight.(lc_ii).Chi2;
    
    if strcmpi(lc_ii, 'traffic')
       combis_ii = combis2;
       chi12     = combis_ii(:,3:4);
       weights_ii= assign_weight(chi12, Chi1_b, Chi2_b, W_b);
    else
       combis_ii = combis1;
       chi12     = combis_ii(:,3:4);
       weights_ii= assign_weight(chi12, Chi1_b, Chi2_b, W_b);
    end
    
    loc1 = size(combis,1) + 1;
    loc2 = loc1 + size(combis_ii,1) - 1;
    load_combs_all(loc1:loc2) = load_combs(ii);
    
    combis   = [combis; combis_ii];
    weights  = [weights; weights_ii];  
end

idx             = cellfun(@isempty, load_combs_all);
load_combs_all(idx) = []; 

combis_colnames = {'d', 'f_cck', 'chi1', 'chi2', 'rho', 'd_lower', 'a_to_d_ratio'};
n_ds            = size(combis, 1);

DS.weights_combis = weights;

Prob_keep       = Prob;

%--------------------------------------------------------------------------
% Semi-probabilistic design & save the design scenarios
%--------------------------------------------------------------------------
% loop over the design scenarios
% for consistency the mean values are also matched with the representative
% values
t_ds = tic;
for ii = 1:n_ds
    Prob_ii             = Prob_keep;
    
    d_repr_ii           = combis(ii,1);
    f_cck_ii            = combis(ii,2);
    chi1_ii             = combis(ii,3);
    chi2_ii             = combis(ii,4);
    rho_repr_ii         = combis(ii,5);
    d_lower_ii          = combis(ii,6);
    a_to_d_ratio_ii     = combis(ii,7);
   
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
    
    Prob_ii             = inv_design(free_par, fix_par, Prob_ii, Prob_actions, Options, load_combs_all{ii});
    Prob(ii)            = Prob_ii;
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
DS.combis           = combis;
DS.combis_colnames  = combis_colnames;
DS.t_elapsed_sec    = t_elapsed_sec;
DS.load_combs_all   = load_combs_all;

end