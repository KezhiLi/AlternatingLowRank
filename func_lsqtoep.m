
%--------------------------------------
function [X] = func_lsqtoep(X)

[n,m] = size(X);
k = min(n,m);
x = zeros(m,1);
y = zeros(n,1);

for i = 1:(k-1)
    sy = 0;
    sx = 0;
    for j = 1:i
        sy = sy + X(n-i+j,j);
        sx = sx + X(j,m-i+j);
    end
    y(n+1-i) = sy/i;
    x(m+1-i) = sx/i;
end


if n <= m
    for i = 1:(m-n+1)
        s = 0;
        for j = 1:n
            s = s + X(j,i+j-1);
        end
        x(i) = s/n;
    end
    y(1) = x(1);
end

if m < n
    for i = 1:(n-m+1)
        s = 0;
        for j = 1:m
            s = s + X(i+j-1,j);
        end
        y(i) = s/m;
    end
    x(1) = y(1);
end
X = toeplitz(y,x);

end
