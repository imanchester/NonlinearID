function [Y,yss] = empirical_hankel(y,i)

[k,s] = size(y);

j = s - 2*i + 1;

% Form the Hanekl matrices
Y = zeros(2*k*i,j);

for p = 1:j   
    Y(:,p) = y((p-1)*k+1:(p-1)*k+2*i*k)';          
end

yss = y(:,i+1:s-i+1);