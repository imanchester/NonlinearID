function res = buildModel_2(prg,x,u,e_monos,f_monos,g_monos)

% 29/11/2016 - the main purpose of this function is to allow direct
% feedthrough of control inputs (i.e. inputs in the output map).


%% Extract dimensions

nx = length(e_monos);
ny = length(g_monos);


%% Remove redundant monomials

res_prepro = lr_nonlin_preprocessModel_3(prg,x,u,e_monos,f_monos,g_monos,nx,0,ny); % JU2911

spot_monos_e = res_prepro.e_monos;
spot_monos_f = res_prepro.f_monos;
spot_monos_g = res_prepro.g_monos;

%% Generate necessary matrices
% -------------------------------------------------------------------------
xu2v = map_xu2v(x,u);

[epow,eco,nem] = mono2mat(spot_monos_e,x);
[fpow,fco,nfm] = fg_mono2mat(spot_monos_f,xu2v,x,u);
[gpow,gco,ngm] = fg_mono2mat(spot_monos_g,xu2v,x,u);

% A couple of quick checks:
if isempty(eco)
    error('Please review your choice of monomials for e.')
end
if isempty(fco)
    error('Please review your choice of monomials for f.')
end
if isempty(gco)
    error('Please review your choice of monomials for g.')
end

% Build anonymous functions
% -------------------------------------------------------------------------
e_monos = @(x) gen_e_monos(x,epow,eco,nx,nem);
f_monos = @(x,u) gen_e_monos(xu2v*[x;u],fpow,fco,nx,nfm);
g_monos = @(x,u) gen_e_monos(xu2v*[x;u],gpow,gco,ny,ngm); % JU2911

% e_monos(x)
% f_monos(x,u)
% g_monos(x)

% Build models as a function of the coefficients:
% -------------------------------------------------------------------------
e = @(x,c) e_monos(x)*c;
f = @(x,u,c) f_monos(x,u)*c;
g = @(x,u,c) g_monos(x,u)*c;

%% Results

res.e = e;
res.f = f;
res.g = g;

res.e_monos = e_monos;
res.f_monos = f_monos;
res.g_monos = g_monos;

end

