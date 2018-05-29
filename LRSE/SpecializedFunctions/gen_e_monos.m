function monos = gen_e_monos(x,pow,co,nx,nmono)

% This is the one that I always used to use:
% monos = reshape(co*prod(repmat(x',nmono,1).^pow,2),nx,nmono);

% Update: 6/12/2017
% For some reason this stopped working (although, I only noticed it for
% scalar systems, so there is the possibility that it never actually worked
% for scalar systems).

if size(x,1) > 1
    monos = reshape(co*prod(repmat(x',nmono,1).^pow,2),nx,nmono);
else
    monos = reshape(co*repmat(x',nmono,1).^pow,nx,nmono);
end

end

