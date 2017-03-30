%Generate pseduoinverse of hankel_struct
%for n times p hankel

function iS = hankel_struct_inv(n,p)

N = n+p-1;  %Signal length

    iS = sparse(N,n*p);
    for k = 1:N
        if k < n
            for j = 1:k
                i = k + 1 - j;
                iS(k,i+(j-1)*n) = 1/k;
            end
        end
        if (k >= n) && (k <= p)
            for i = 1:n
                j = k + 1 - i;
                iS(k,i+(j-1)*n) = 1/n;
            end
        end
        if k > p
         	for j = (k+1-n):p
                i = k + 1 - j;
                iS(k,i+(j-1)*n) = 1/(n+p-k);
            end
        end
    end
end