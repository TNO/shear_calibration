% Prepare input for reliability analysis in Matlab (FERUM)
% - You have to have the pacal results in the `results` folder
%
% Downscaling:
% - done in Gumbel space
% - evenly spaced points on the -log(-log(cdf)) interval
% - interpolate to x
%
% This file is meant to be run from its location.

clearvars
close all
clc

action      = {'wind', 'snow', 'dummy1'};
ID          = [201, 301, 1001];

prec        = 8;

% relative to this file's location
pacal_res_dir = 'results\';
ferum_tmp_dir = '..\..\calibration\tmp\';

n_a         = length(action);

% number of discretization points
n_downscale = 1e3;

for ii = 1:n_a
    % .....................................................................
    % LOAD
    % .....................................................................
    action_ii = action{ii};
    ID_ii       = ID(ii);
    if ID_ii == 1001
        mean_d1 = 1;
        cov_d1  = 0.4;
        xx = linspace(0, 20, 1e3);
        M(:,1) = xx;
        M(:,2) = lognormpdf(xx, mean_d1, cov_d1);
        M(:,3) = lognormcdf(xx, mean_d1, cov_d1);
    else
        % load python PaCaL results
        M           = dlmread([pacal_res_dir, 'unit_', action_ii, '_pacal.txt']);
    end
    % .....................................................................
    % REMOVE DUPLICATE
    % .....................................................................
    x_grid      = M(:,1);
    cdf         = M(:,3);
    pdf         = M(:,2);
    
    [~, idx]    = unique(pdf);
    x_grid      = x_grid(idx);
    pdf         = pdf(idx);
    cdf         = cdf(idx);
    
    [~, idx]    = unique(cdf);
    x_grid      = x_grid(idx);
    pdf         = pdf(idx);
    cdf         = cdf(idx);

    mm          = trapz(x_grid, pdf.*x_grid);
    vv          = trapz(x_grid, pdf.*(x_grid-mm).^2);
    
    % .....................................................................
    % CHARACTERISTIC VALUE
    % .....................................................................
    % WARNING - SAME INPUT AS IN PYTHON-PACAL!!
    switch lower(action_ii)
        case 'wind'
            % Gumbel (qref1)
            X1_mean     = 1.00;
            X1_cov      = 0.27;
            X1_P_char   = 0.98;
            X1_char     = gumbelinvcdf(X1_P_char, X1_mean, X1_mean*X1_cov, 'mom');

            % Lognormal (ce)
            X2_mean     = 1.00;
            X2_cov      = 0.15;
            X2_P_char   = 0.94;
            X2_char     = lognorminv(X2_P_char, X2_mean, X2_cov);
            
            % Gumbel (cpe)
            X3_mean     = 1.00;
            X3_cov      = 0.20;
            X3_P_char   = 0.78;
            X3_char     = gumbelinvcdf(X3_P_char, X3_mean, X3_mean*X3_cov, 'mom');
            
            % Lognormal (cd)
            X4_mean     = 1.00;
            X4_cov      = 0.15;
            X4_P_char   = 0.50;
            X4_char     = lognorminv(X4_P_char, X2_mean, X2_cov);
            
            x_char      = X1_char*X2_char*X3_char*X4_char;
            P_char      = interp1(x_grid, cdf, x_char);
        case 'snow'
            % Gumbel (sground1) 
            X1_mean     = 1.00;
            X1_cov      = 0.60;
            X1_P_char   = 0.98;
            X1_char     = gumbelinvcdf(X1_P_char, X1_mean, X1_mean*X1_cov, 'mom');

            % normal (ce)
            X2_mean     = 1.0;
            X2_cov      = 0.15;
            X2_P_char   = 0.50;
            X2_char     = norminv(X2_P_char, X2_mean, X2_cov*X2_mean);
            
            x_char      = X1_char*X2_char;
            P_char      = interp1(x_grid, cdf, x_char);
        case 'dummy1'
            P_char      = 0.98;
            x_char      = lognorminv(P_char, mean_d1, cov_d1);
        otherwise
            error(['Unknown action: ', action_ii])
    end
    
    disp('----------------------')
    disp(['action:', action_ii])
    disp(['mean:', num2str(mm)])
    disp(['mean:', num2str(sqrt(vv))])
    disp(['x_char:', num2str(x_char)])
    disp(['P_char:', num2str(P_char)])
    
    % .....................................................................
    % DOWNSCALE
    % .....................................................................
    
    tcdf        = -log(-log(cdf));
    
    idx         = imag(tcdf) | isinf(tcdf) | isnan(tcdf);
    tcdf        = tcdf(~idx);
    
    tp1         = min(tcdf);
    tp2         = max(tcdf);
    down_tcdf   = linspace(tp1, tp2, n_downscale);
    down_cdf    = exp(-exp(-down_tcdf));
    down_x_grid = interp1(cdf, x_grid, down_cdf, 'linear', 'extrap');
    down_pdf    = interp1(x_grid, pdf, down_x_grid, 'linear', 'extrap');
    
    % enforce precision already here to avoid duplicates - lame
    ee          = floor(log10(abs(down_x_grid)));
    tx          = down_x_grid./(10.^ee);
    tx          = round(tx*10^prec)/10^prec;
    down_x_grid = tx.*(10.^ee);
    
    [~, idx]    = unique(down_x_grid);
    down_x_grid = down_x_grid(idx);
    down_pdf    = down_pdf(idx);
    down_cdf    = down_cdf(idx);
    
    
%     [~, idx]    = unique(down_pdf);
%     down_x_grid = down_x_grid(idx);
%     down_pdf    = down_pdf(idx);
%     down_cdf    = down_cdf(idx);
%     
    % .....................................................................
    % VISUALIZE ERROR DUE TO DOWNSCALING
    % .....................................................................
    
    figure
    plot(x_grid, -log(-log(cdf)))
    hold on
    plot(down_x_grid, -log(-log(down_cdf)), '--')
    
    e_pdf = interp1(down_x_grid, down_pdf, x_grid) - pdf;
    e_cdf = interp1(down_x_grid, down_cdf, x_grid) - cdf;
    
    figure
    subplot(2,1,1)
    plot(x_grid, e_pdf)
    xlim([1,max(x_grid)])
    xlabel('x')
    ylabel('abolute error (downscaled vs full)')
    title(['ID: ', num2str(ID_ii), ' pdf'])
    
    subplot(2,1,2)
    plot(x_grid, e_pdf./pdf)
    xlim([1,max(x_grid)])
    xlabel('x')
    ylabel('relative error (downscaled vs full)')
    
    figure
    subplot(2,1,1)
    plot(x_grid, e_cdf)
    xlim([1,max(x_grid)])
    xlabel('x')
    ylabel('abolute error (downscaled vs full)')
    title(['ID: ', num2str(ID_ii), ' cdf'])
    
    subplot(2,1,2)
    plot(x_grid, e_cdf./cdf)
    xlim([1,max(x_grid)])
    xlabel('x')
    ylabel('relative error (downscaled vs full)')

    % .....................................................................
    % SAVE
    % .....................................................................
    
    pdf         = down_pdf(:);
    cdf         = down_cdf(:);
    x_grid      = down_x_grid(:);
    
    x_med       = interp1(cdf, x_grid, 0.5);
    idx_left    = x_grid <= x_med;
    idx_right   = x_grid > x_med;
    
    save([ferum_tmp_dir, 'vector_distr_' num2str(ID_ii) '.mat'], 'x_grid', 'pdf', 'cdf', 'x_char', 'P_char')
    
    % .....................................................................
    % SAVE DOWNSCALED VERSION
    % .....................................................................
    filename = [ferum_tmp_dir, num2str(ID_ii), '_x_grid.txt'];
    dlmwrite(filename,x_grid, 'delimiter', ' ', 'precision', '%.8e') 
    
    filename = [ferum_tmp_dir, ferum_tmp_dir, num2str(ID_ii), '_pdf.txt'];
    dlmwrite(filename, pdf, 'delimiter', ' ', 'precision', '%.8e')

    filename = [ferum_tmp_dir, num2str(ID_ii), '_cdf_left.txt'];
    dlmwrite(filename, cdf(idx_left), 'delimiter', ' ', 'precision', '%.8e')
    
    filename = [ferum_tmp_dir, num2str(ID_ii), '_cdf_right.txt'];
    dlmwrite(filename, 1-cdf(idx_right), 'delimiter', ' ', 'precision', '%.8e')
%     fid = fopen(filename, 'wt');
%     if fid ~= -1
%         for row = 1 : size(A, 1)
%             fprintf(fid, '%.20e, %.20e;\n', A(row, 1), A(row, 2));
%         end
%         fclose(fid);
%     end

end

close all
% interp1(x_grid, cdf, 1)