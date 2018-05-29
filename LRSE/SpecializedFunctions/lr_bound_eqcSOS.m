function J = lr_bound_eqcSOS(pn,pstar,blks,sos_opt)
% Function to simply compute the bound.
% Probably redundant, but we'll use it for now...

if isfield(blks,'lrnorm')
    lrnorm = blks.lrnorm;
else
    lrnorm = 1;
end

% Convert pn to p
p = pstar + sos_opt.N*pn;

%% Stuff you can precompute...

eps_e = blks.eps_e;
eps_f = blks.eps_f;
eta_g = blks.eta_g;

Fbar_e = blks.Fbar_e;
Fbar_f = blks.Fbar_f;
Gbar_g = blks.Gbar_g;

% Smuggle y in here:
y = blks.yd;

% Number of coefficients:
ncoef(1) = length(eps_e);
ncoef(2) = length(eps_f);
ncoef(3) = length(eta_g);

%% Assemble block matrices
% -------------------------------------------------------------------------

% Build Fbar, eps
% ---------------

% Slight mod:
    Fbar = Fbar_e{1}*p(1);
    eps = eps_e{1}*p(1);

% Contribution: e
    for i = 2:ncoef(1)
        Fbar = Fbar + Fbar_e{i}*p(i);
        eps = eps + eps_e{i}*p(i);
    end

% Contribution: f
    for i = 1:ncoef(2)
        Fbar = Fbar + Fbar_f{i}*p(ncoef(1)+i);
        eps = eps + eps_f{i}*p(ncoef(1)+i);
    end

% Build Gbar, eta
% ---------------

    Gbar = Gbar_g{1}*p(sum(ncoef(1:2))+1);
    eta = eta_g{1}*p(sum(ncoef(1:2))+1) - y(:);

% Contribution: g
    for i = 2:ncoef(3)
        Gbar = Gbar + Gbar_g{i}*p(sum(ncoef(1:2))+i);
        eta = eta + eta_g{i}*p(sum(ncoef(1:2))+i);
    end


%% Solve for Delta*
% -------------------------------------------------------------------------

% For now let's assume that the multiplier is simply 2*Delta, i.e. Q = 2I.
    Delta = ((-(eps' + eta'*Gbar))/(Gbar'*Gbar - (Fbar + Fbar')))';
    
%% Objective function
% ------------------------------------------------------------------------- 
    J = ((Gbar*Delta + eta)'*(Gbar*Delta + eta) - 2*Delta'*(Fbar*Delta - eps))/lrnorm;  


end

