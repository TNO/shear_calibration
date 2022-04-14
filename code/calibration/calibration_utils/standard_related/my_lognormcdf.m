function P = my_lognormcdf(x, x_mean, x_cov)
    mu_lognorm      = log(x_mean.^2./sqrt((x_mean.*x_cov).^2 + x_mean.^2));
    sigma_lognorm   = sqrt(log((x_mean.*x_cov).^2./x_mean.^2 + 1));
    P               = logncdf(x, mu_lognorm, sigma_lognorm);
end