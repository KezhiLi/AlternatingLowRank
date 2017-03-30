%Least square (Frobeinius) Hankel approximation
%Martin Sundin, 2014-03-24

function hhat = hankel_proj(H)

[N,M] = size(H);
x = zeros(N,1);
y = zeros(M,1);
k = min(M,N);
for i = 1:(k-1)
    s = 0;
    for j = 1:i
        s = s + H(i-j+1,j);
    end
    x(i) = s/i;
end
for i = 1:(k-1)
    s = 0;
    for j = 1:i
        s = s + H(N+1-j,M-i+j);
    end
    y(M+1-i) = s/i;
end
if N >= M
    for i = M:N
        s = 0;
       	for j = 1:M
          	s = s + H(1+i-j,j);
        end
      	x(i) = s/M;
    end
   	y(1) = x(N);
end
if M > N
   	for i = 1:(M-N+1)
      	s = 0;
      	for j = 1:N
          	s = s + H(N+1-j,i+j-1);
        end
   	y(i) = s/N;
    end
    x(N) = y(1);
end

hhat = [x;y(2:end)];


end
