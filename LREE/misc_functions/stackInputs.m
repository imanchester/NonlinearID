function us = stackInputs(u,n)

[nu,T] = size(u);

us = zeros(nu*n,T);

for i = 1:n
    us((i-1)*nu+1:nu*i,:) = [zeros(nu,i-1), u(:,1:end-i+1)];
end


end

