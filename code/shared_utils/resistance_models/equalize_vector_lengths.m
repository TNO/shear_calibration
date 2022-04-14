function varargout = equalize_vector_lengths(varargin)
    args = varargin;
    n = nargin;
    [n_max, arg_max] = max(cellfun(@length, args));
    reference_arg = args{arg_max};

    outs = cell(size(args));
    for ii = 1:n
        arg = args{ii};
        if length(arg) == n_max
            outs{ii} = arg;
        elseif length(arg) == 1
            outs{ii} = arg * ones(size(reference_arg));
        else
            error(['Each input argument must have length 1 or the length' ...
                'of the longest input argument to be able to equalize ' ...
                'their lenghts.'])
        end
    end

    varargout = outs;
end