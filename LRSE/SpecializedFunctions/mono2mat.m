function [pow,co,nmonos] = mono2mat(monos,x)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[nx,~] = size(x);

nmonos = 0; % Number of monomials

del = 0.1;

pow = [];
co = [];

    for i = 1:length(monos)

    % Extract monomials:
        tmp = monos{i};

    % Keep track of number of monomials:
%         nm = length(tmp);
%         nmonos = nmonos + nm;
        
%   Unfortunately, we need some trickery:
        tmp = [tmp, del*x'];

    % Decompose monomials:    
        [v,p,c] = decomp(tmp); 

%   Now we find the columns of the coefficients that are just del, these can be eliminated:
        fkrs = find(sum(c,1)==del);
        p(fkrs,:) = [];
        
    % Better method:
        nm = size(p,1);
        nmonos = nmonos + nm;        
        
    % Store useful powers:
        pow = [pow;p]; 

    % Create the coefficient mapping:
        tc = zeros(nm*nx,nm);
        tc(i:(nm+1)*nx:end-nx + i) = 1;

        co{i} = tc;

    end

    co = blkdiag(co{:});


end

