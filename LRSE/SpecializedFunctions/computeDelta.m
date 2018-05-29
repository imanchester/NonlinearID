function delta = computeDelta(res_lr,y)
% This function assumes that you have already run the algorithm

    pk = res_lr.theta_ld; % The (potentially) lower dimensional theta
    pstar = res_lr.theta_star; 
    
    blks = res_lr.blks;
    sos_opt = res_lr.sos;
    
    eps_e = blks.eps_e;
    eps_f = blks.eps_f;
    eta_g = blks.eta_g;

    Fbar_e = blks.Fbar_e;
    Fbar_f = blks.Fbar_f;
    Gbar_g = blks.Gbar_g;

    % Number of coefficients:
    ncoef(1) = length(eps_e);
    ncoef(2) = length(eps_f);
    ncoef(3) = length(eta_g);

    Fbar_e_full = blks.Fbar_e_full;
    eps_e_full = blks.eps_e_full;
    Fbar_f_full = blks.Fbar_f_full;
    eps_f_full = blks.eps_f_full;
    Gbar_g_full = blks.Gbar_g_full;
    eta_g_full = blks.eta_g_full;
    d1 = blks.d1;
    d2 = blks.d2;    
    
    nth = sum(ncoef);

    % SOS
    Almi = sos_opt.Almi;
    blmi = sos_opt.blmi;
    N1 = sos_opt.N1;
    N2 = sos_opt.N2;
    N = sos_opt.N;
    dimNL = sos_opt.dimNL;
    nphi = sos_opt.nphi;    


    delta = compDelta(pk);

%% Compute Delta - straight from the normal functions

    function [Delta] = compDelta(pn)

% FIRST compute Delta
% -------------------------------------------------------------------------

% First obtain p (which represents the model parameters):      
    p = (pstar(1:end-dimNL) + N1*pn);
       
% Construct block matrices: it seems like the loops are actually faster.   
    Fbar = Fbar_e{1}*p(1);
    
% Contribution: e
    for i = 2:ncoef(1)
        Fbar = Fbar + Fbar_e{i}*p(i);
    end
    
    for i = 1:ncoef(2)
        Fbar = Fbar + Fbar_f{i}*p(ncoef(1)+i);
    end

    eps = [eps_e_full eps_f_full]*p(1:sum(ncoef(1:2))); 

    Gbar = Gbar_g{1}*p(sum(ncoef(1:2))+1);
% Contribution: g
    for i = 2:ncoef(3)
        Gbar = Gbar + Gbar_g{i}*p(sum(ncoef(1:2))+i);
    end
    
    eta = eta_g_full*p(sum(ncoef(1:2))+1:sum(ncoef)) - y(:);

% Solve for Delta*
    Delta = ((-(eps' + eta'*Gbar))/(Gbar'*Gbar - (Fbar + Fbar')))';     
    
    end



end

