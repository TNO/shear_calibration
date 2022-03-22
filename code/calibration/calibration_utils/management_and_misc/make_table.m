% Display the updadted probabilistic models in the Command Window
%
%
%TODO: 
% * sub_fields should be the only the common field names on all the first
% level fields

function Tab = make_table(Struct)

fields      = fieldnames(Struct);
n_f         = length(fields);

% assuming that the same subfield names are present for all fields
field_ii = 2;
sub_fields  = fieldnames(Struct.(fields{field_ii}));

n_sf = length(sub_fields);
bm_keep = false(n_sf, 1);
% keep only the numerical ones
for jj = 1:n_sf
    sub_field_value = Struct.(fields{field_ii}).(sub_fields{jj});
    if isa(sub_field_value, 'double')
        bm_keep(jj) = true;
    end
end
sub_fields = sub_fields(bm_keep);

n_sf        = length(sub_fields);

M           = nan(n_f, n_sf);

% put the values into a matrix
for ii = 1:n_f
    for jj = 1:n_sf
        M(ii,jj) = Struct.(fields{ii}).(sub_fields{jj});
    end
end

Tab = array2table(M, 'RowNames', fields, 'Variablenames', sub_fields);

end