function [yvec] = mean_k_width(xvec,k)

% y(1) is mean(x(1:5)
% y(2) is mean(x(6:10)
% y(3) is mean(x(11:15))
% etc

n_bin = fix(numel(xvec)/k);

% k=5
yvec = zeros(n_bin,1);
j = -(k-1);
for bin_c = 1:n_bin
    j = j+k;
    yvec(bin_c,1) = mean(xvec(j:j+(k-1)));
end



end %end function