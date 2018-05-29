function H = lr_sosHess_hybridSOS_2(p,pstar,sos_opt,varargin)
% 

%% Precomputed quantities

Almi = sos_opt.Almi;
blmi = sos_opt.blmi;
% Msos = sos_opt.s2v2s;

Mbar = sos_opt.Mbar;

%  = sos_opt.nQ;
N = sos_opt.N;
nphi = sos_opt.nphi;
lmitol = sos_opt.lmitol;

%% Compute Zinv if necessary

nv = length(varargin);

if nv == 1
    
    Zinv = varargin{1}; % User supplied Zinv 
    
else

    pf = pstar + N*p; % This is the full parameter vector
    
%   Compute Zinv for oursleves:
    Z = mss_v2s(Almi*pf + blmi) - lmitol*eye(nphi); % Note the tolerance
    Zinv = inv(Z);

end

%% Compute the Hessian

% Let's assume that we have inv(Z) at this point...
% [nz,~] = size(Zinv);
nz = nphi;
% 
H = zeros(nz^2);
% 
% % I think using two for loops is actually faster than alternatives...
% for i = 1:nz
%     for j = 1:nz      
% %       Install the i,jth block  
%         H((i-1)*nz+1:i*nz,(j-1)*nz+1:j*nz) = Zinv(:,j)*Zinv(i,:);
%     end
% end

for i = 1:nz
    for j = 1:i
%       Compute the i,j th block entry
        tmp = Zinv(:,j)*Zinv(i,:);
%       Install the i,jth block  
        H((i-1)*nz+1:i*nz,(j-1)*nz+1:j*nz) = tmp;
%       Install the j,i th block
        H((j-1)*nz+1:j*nz,(i-1)*nz+1:i*nz) = tmp' ;
    end
end

% Does C make it any faster? Maybe a bit.
% It depends: for the 4rd order system it was faster, for the 5th order it
% was slower; go figure.
% H = buildhess_mex(Zinv,nz); 

% A little trickery required for the Hessian:
% Abar = Msos*Asos;
% H = Abar'*H*Abar;
% Mbar = Msos*Almi*N;

% H = sparse(H);

H = Mbar'*H*Mbar;

end

