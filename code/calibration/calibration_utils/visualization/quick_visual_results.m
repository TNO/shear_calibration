function quick_visual_results(Results, beta_t, title_str)

if nargin < 3
    title_str = '';
end

nn = length(Results.beta);

figure
plot(Results.beta, '.-', 'MarkerSize', 20)
hold on
plot(1:nn, beta_t*ones(nn,1), 'red')
ylabel('$\beta$')
xlabel('Design situation')
title(title_str)

prettify(gcf)

end