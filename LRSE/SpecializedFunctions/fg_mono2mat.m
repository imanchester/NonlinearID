function [pow,co,nmonos] = fg_mono2mat(monos,xu2v,x,u)

% 20/11/2016 - one function for both f and g.

[nx,~] = size(x);
[nu,~] = size(u);
nv = nx + nu;

v = xu2v*[x;u];

nmonos = 0; % Number of monomials

n_ = length(monos); % Wildcard: nx for f, ny for g

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
        tmp = [tmp, del*v'];

    % Decompose monomials:    
        [~,p,c] = decomp(tmp); 

%   Now we find the columns of the coefficients that are just del, these can be eliminated:
        fkrs = find(sum(c,1)==del);
        p(fkrs,:) = [];
        
    % Better method:
        nm = size(p,1);
        nmonos = nmonos + nm;            
        
    % Store useful powers:
        pow = [pow;p]; 

    % Create the coefficient mapping:
        tc = zeros(nm*n_,nm);
        tc(i:(nm+1)*n_:end-n_ + i) = 1;

        co{i} = tc;

    end

    co = blkdiag(co{:});


end

