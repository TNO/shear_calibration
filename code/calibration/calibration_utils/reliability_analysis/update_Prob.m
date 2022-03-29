% Update the Prob structure
%
%SYNOPSYS
%
%
%
%TODO
% * check input
% * Probably we do not need these
%   - add standard deviation field
%   - add representative value field
%   - add design value field

function Prob = update_Prob(Prob, verbose)

% -------------------------------------------------------------------------
% INITIALIZE % PRE-PROCESS
% -------------------------------------------------------------------------

if nargin < 2
    verbose = 0;
end

% -------------------------------------------------------------------------
% UPDATE THE PAR-NONPAR DISTRIBUTIONS AND INITIALIZE ALL
% -------------------------------------------------------------------------
rv_names    = fieldnames(Prob);
n_rv        = length(rv_names);
marg        = nan(n_rv, 9);

for ii = 1:n_rv
    name    = rv_names{ii};
    X       = Prob.(name);
    if ~isfield(X, 'cov')
        X.cov = NaN;
    end
    if ~isfield(X, 'std')
        X.std = NaN;
    end
    if ~isfield(X, 'P_repr')
        X.P_repr = NaN;
    end
    if ~isfield(X, 'gamma')
        X.gamma = NaN;
    end

    if X.dist == 32
        load(['tmp\vector_distr_', num2str(X.dist_ID), '.mat'], 'x_grid', 'pdf')

        X.mean      = trapz(x_grid, pdf.*x_grid);
        X.std       = sqrt(trapz(x_grid, pdf.*(x_grid-X.mean).^2));
        X.cov       = X.std/X.mean;
        X.shift     = 0;
        X.scale     = 1;

        X.mean_base = X.mean;
        X.std_base  = X.std;
        X.cov_base  = X.std/X.mean;

        marg(ii,:) = [X.dist,  X.mean,  X.mean*X.cov,  X.mean,  X.dist_ID,  X.shift,  X.scale,  NaN, 0];
    elseif X.dist == 33
        [~, x_grid, pdf]    = hardcoded_pdf(1, X.dist_ID);

        mean_i      = trapz(x_grid, pdf.*x_grid);
        stdv_i      = sqrt(trapz(x_grid, pdf.*(x_grid-mean_i).^2));
        X.mean_base = mean_i;
        X.std_base  = stdv_i;
        X.cov_base  = stdv_i/mean_i;
        X.shift     = 0;
        X.scale     = 1;

        X.mean      = X.scale*mean_i + X.shift;
        X.std       = X.scale*stdv_i;
        X.cov       = X.std/X.mean;
    else
        % add dummy
        X.shift     = NaN;
        X.scale     = NaN;
        X.dist_ID   = NaN;

        % check if for non-parametric distributions only one of `cov` and
        % `std` is different than zero
        if ~isnan(X.cov) && ~isnan(X.std)
            error('Only one of the `X.cov` `X.std` pairs must differ from NaN.')
        end
    end
    Prob.(name)     = X;
end

% TODO: check if marg is needed here at all
% needed for the par-nonpar distributions (FERUM like)
marg = marg(~isnan(marg(:,1)),:);
distribution_parameter(marg);


% -------------------------------------------------------------------------
% REPRESENTATIVE VALUES
% -------------------------------------------------------------------------

for ii = 1:n_rv
    name    = rv_names{ii};
    X       = Prob.(name);
    if isfield(X, 'mean2repr')
        [~, X] = X.mean2repr(X.mean, X);
    else
        if ~isfield(X, 'repr')
            X.repr = NaN;
        end
    end

    Prob.(name)     = X;
end

% -------------------------------------------------------------------------
% COLLECT % POST-PROCESS
% -------------------------------------------------------------------------

if verbose > 0
    Tab = make_table(Prob);
    disp('------------------------------------------------------------')
    disp('Probabilistic models:')
    disp(Tab)
end

end