function  plot_beta_chi(Results, DS, Options)

    betas       = Results.beta;
    chi1        = DS.combis(:,3);
    chi2        = DS.Range.chi2;
    combis      = DS.combis;
    
    nprop       = [];
    ds_names    = fieldnames(DS.Range);
    for ii = 1:length(ds_names)
        nprop = [nprop; length(DS.Range.(ds_names{ii}))];
    end
    n_ds        = sum(nprop);
    n_dpar      = size(combis,2);
    
    load_combs  = Options.load_combs;
    
    figure
    hold on
    for jj = 1:length(load_combs)
        plot(chi1(jj * n_ds - nds+1 : jj * n_ds), betas(jj*n_ds - n_ds+1 : jj*n_ds))
    end
    
end

