function F = f_util(c,hs,theta,sigma)

kappa = 0;
uinner=(c^theta)*((kappa+hs)^(1-theta)); 
F=(uinner^(1-sigma))/(1-sigma);


end