% Update the `DS` (design scenario) structure.
% 
% It does:
% - drop design scenarios that are not relevant for the selected resistance
%   model
%

function DS = update_DS(DS, Options)

verbose = Options.verbose;
resistance_model = Options.resistance_model;

if strcmp(resistance_model, 'ec2_codified_2019')
    % size of the maximum sieve gird used for aggregates
    d_lower = DS.Range.d_lower;
    if length(d_lower) > 1
        DS.Range.d_lower = d_lower(1);
        
        if verbose > 0
            warning(['`d_lower` does not affect `', resistance_model, '`; ' ...
                'therefore, its discretization is dropped from forming ' ...
                'design scenarios (DS) to reduce the computational ' ...
                'effort.'])
        end
    end
    
    % a-d ratio
    a_to_d_ratio = DS.Range.a_to_d_ratio;
    if length(a_to_d_ratio) > 1
        DS.Range.a_to_d_ratio = a_to_d_ratio(1);
        
        if verbose > 0
            warning(['`a_to_d_ratio` does not affect `', resistance_model, '`; ' ...
                'therefore, its discretization is dropped from forming ' ...
                'design scenarios (DS) to reduce the computational ' ...
                'effort.'])
        end
    end
end

end