function [jequaloneDistz] = discretize_normal(z_grid_log,sigma_z1)

% z_grid: Grid for productivity shock z
% sigma_z1: Standard deviation
% Assumption: log(z_{t+1})=rho*log(z_t)+eps, eps~N(0,sigma_z1^2),t=1,..,T
% log(z_1) is itself a Normal random variable with mean 0 and standard
% deviation sigma_z1. This function discretizes the distribution of log(z_1)

n_z = length(z_grid_log);
jequaloneDistz = zeros(n_z,1);

z_first =  0.5*(z_grid_log(1)+z_grid_log(2));
jequaloneDistz(1) = normcdf(z_first/sigma_z1);
for ii=2:n_z-1
    z_lb = 0.5*(z_grid_log(ii)+z_grid_log(ii-1))/sigma_z1;
    z_ub = 0.5*(z_grid_log(ii+1)+z_grid_log(ii))/sigma_z1;
    jequaloneDistz(ii) = normcdf(z_ub)-normcdf(z_lb);
end
z_last =  0.5*(z_grid_log(n_z-1)+z_grid_log(n_z));
jequaloneDistz(n_z) = 1-normcdf(z_last/sigma_z1);

end %end function