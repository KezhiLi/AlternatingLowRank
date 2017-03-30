function [hh] = func_lsqhank(hin)
    [m,n]=size(hin);
    hin=fliplr(hin);
    tmp = zeros(m+n-1,1);
    for k=-(m-1):(n-1)
        tmp(k+m) = mean(diag(hin,k));
    end
    tmp=flipud(tmp);
    hh = hankel(tmp(1:m), tmp(m:end));
end