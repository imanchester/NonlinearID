function res = lrse_minftau(x,y,u,p0,pstar,tau,B0,sos_opt,blks,varargin)


%% Precomputed quantities

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

nth = sum(ncoef);

% SOS
Almi = sos_opt.Almi;
blmi = sos_opt.blmi;
N1 = sos_opt.N1;
N2 = sos_opt.N2;
N = sos_opt.N;
dimNL = sos_opt.dimNL;
nphi = sos_opt.nphi;

dim_pmod = length(pstar) - dimNL;

lrnorm = blks.lrnorm;

% Tolerance for Q
lmitol = sos_opt.lmitol;

Fbar_e_full = blks.Fbar_e_full;
eps_e_full = blks.eps_e_full;
Fbar_f_full = blks.Fbar_f_full;
eps_f_full = blks.eps_f_full;
Gbar_g_full = blks.Gbar_g_full;
eta_g_full = blks.eta_g_full;
d1 = blks.d1;
d2 = blks.d2;

% Trace constraints
Atr = sos_opt.Atr;
trlim = sos_opt.trlim;

%% Search properties

nv = length(varargin);

if nv == 1
    options = varargin{1};
else
    options = [];
end

% Default parameters for the Wolfe line search   
c1 = 1e-4;
c2 = 0.9;
LS_interp = 2;
LS_multi = 0;
LS_max = 25;
LS_debug = 0;
LS_doPlot = 0;

if isfield(options,'maxIters')
    maxIters = options.maxIters;
else
    maxIters = 1e4;
end

tolmod = 1;
tolmod2 = 1e5;

if isfield(options,'progTol')
    progTol = options.progTol;
else
    progTol = tolmod*1e-15;
end

if isfield(options,'optTol')
    optTol = options.optTol;
else
    optTol = tolmod*1e-15;
end

if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 0;
end

% These should always be set:
trcnst = options.trcnst;
earlyterm = options.earlyterm;

%% Initialize 

p_k = p0;

% Function handle for line search
funObj = @(p) target(p);


%% Minimize f_tau

Js = zeros(1,100);
J_old = inf; 

% Store the times per iteration,
newtontimes = zeros(1,1e3);

% Store no. of line search iterations/line search function evaluations
lsiters = zeros(1,1e3); 

searchtimes = zeros(1,1e3);
hessiantimes = zeros(5,1e3);

for k = 1:maxIters

% Stat Newton timer:
tnewt = tic;
    
if k == 1   
    
% Gradient of the Lagrangian:
% Method 1:
    [J_k,dJ_k] = lagrangian(p_k);

% Gradient of the barrier:
    [bar_k,dbar_k] = barrier(p_k);
    
% Total gradient:
    g_k = dJ_k + tau*dbar_k;
    f_k = J_k + tau*bar_k;
    
end

% Compute the Hessian of the barrier function (analytically):
ddbar_k = lr_sosHess_hybridSOS_2(p_k,pstar,sos_opt);

% Update Hessian with trace barrier:
if trcnst
    % ddbar_k(1:nth,1:nth) = ddbar_k(1:nth,1:nth) + tracehessian(p_k);
    ddbar_k = ddbar_k + tracehessian(p_k);
    % Compute the trace of Gram matrix
    fprintf('\t\tTrace of Gram: %.5e.\n',trace(mss_v2s(Almi*(pstar + N*p_k) + blmi)));
end

[ddJ_k,btime] = lagrangianHessian(p_k);

hessiantimes(:,k) = btime;


if verbose >= 2
    fprintf('\tddJ rcond: %.5e, Max param: %.5e.\n',rcond(ddJ_k),max(abs(p_k)))
end

% Form the `total Hessian' (adding the Lagrangian and the barrier
% function):
hess_k = tau*ddbar_k;
% hess_k(1:nth,1:nth) = hess_k(1:nth,1:nth) + ddJ_k;
hess_k = hess_k + ddJ_k;

t1srch = tic;

d_k = - hess_k\g_k; % Approx. Newton step

searchtimes(k) = toc(t1srch);

gtd = g_k'*d_k; % Approx. Newton decrement

% Compute the step length using a Wolfe line search:
t = 1; % Initialize the step length
[t,f_new,g_new,lsfunevals] = WolfeLineSearch(p_k,t,d_k,f_k,g_k,gtd,c1,c2,LS_interp,LS_multi,LS_max,progTol,LS_debug,LS_doPlot,1, funObj);

% Store the number of line search iterations:
% NOTE: we're actually storing the no. of function evaluations in the line
% search.
lsiters(k) = lsfunevals;

% Update the parameter value
p_new = p_k + t*d_k;

% Store the value of the lag
J_new = lagrangian(p_new);
Js(k) = J_new;

% Store iteration time:
newtontimes(k) = toc(tnewt);

if t == 0 
    if verbose >= 2; fprintf('\tNewton Stopping: t = 0 !\n'); end;
    break
end

% Check termination conditions:
if max(abs(g_new)) <= optTol
    if verbose >= 2; fprintf('\tNewton stopping: max abs gradient less than opTol.\n'); end;
    break
end

if max(abs(t*d_k)) <= progTol;
    if verbose >= 2; fprintf('\tNewton stopping: step size below progTol.\n'); end;
    break
end

if (f_k - f_new) <= progTol;
% if abs(f_new-f_k) <= progTol
    if verbose >= 2; fprintf('\tNewton stopping: change in function value less than progTol.\n'); end;
    break
end


if (earlyterm)
    if J_old - J_new <= 1e-10; % Original
        if verbose >= 2; fprintf('\tNewton stopping: change in J less than tolerance.\n'); end;
        break;
    end
    J_old = J_new;
end


% Update gradient, function value and current parameter estimate:
g_k = g_new;
f_k = f_new;
p_k = p_new;
% dJ_k = dJ_new;
% Delta_k = Delta_new;

end

if verbose >= 2; fprintf('\tBFGS iterations: %d.\n',k); end;
if verbose >= 2; fprintf('\tBFGS objective: %.7e.\n',f_new); end;

%% Output results

res.p = p_new;
res.B = 0;
res.ftau = f_new;
res.times = newtontimes(1:k); % Return time/number of iterations
res.lsiters = lsiters(1:k);

res.hessiantimes = hessiantimes(:,1:k);
res.searchtimes = searchtimes(1:k);
% res.btimes = btimes(1:k);

%% Target function
    function [f,df] = target(p)

% NB: gradients are w.r.t. p (corrected in sub-functions).        
        
% Evaluate Lagrangian and gradient:
    [J,dJ] = lagrangian(p);
        
% Evaluate barrier and gradient:
    [bar,dbar] = barrier(p);
            
% Objective:
    f = J + tau*bar;
%     f = tau*J + bar;
        
% Gradient:
    df = dJ + tau*dbar;
%     df = tau*dJ + dbar;
                
    end


%% Lagrangian

    function [J,dJ] = lagrangian(pn)
       
% Full parameter vector
% ------------------------------------------------------------------------- 
% This function accepts pn (spanning null(A)) and returns derivative w.r.t
% pn.
    
% First obtain p (which represents the model parameters):      
    p = (pstar(1:end-dimNL) + N1*pn);
        
        
% Assemble block matrices
% -------------------------------------------------------------------------

% Slight mod:
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
% -------------------------------------------------------------------------
    Delta = ((-(eps' + eta'*Gbar))/(Gbar'*Gbar - (Fbar + Fbar')))';
    
% Objective function
% ------------------------------------------------------------------------- 
    J = ((Gbar*Delta + eta)'*(Gbar*Delta + eta) - 2*Delta'*(Fbar*Delta - eps))/lrnorm;    

% Compute gradient of the objective function (no barrier)
% -------------------------------------------------------------------------
    dJ = zeros(dim_pmod,1);

    dJ(1:ncoef(1)) = -2*(Delta'*(reshape(Fbar_e_full*Delta,d1,ncoef(1)) - eps_e_full))'; 
    dJ(ncoef(1)+1:sum(ncoef(1:2))) = -2*(Delta'*(reshape(Fbar_f_full*Delta,d1,ncoef(2)) - eps_f_full))'; 
    dJ(sum(ncoef(1:2))+1:sum(ncoef)) = 2*(Gbar*Delta + eta)'*(reshape(Gbar_g_full*Delta,d2,ncoef(3)) + eta_g_full);

% This is the gradient w.r.t to p.
% -------------------------------------------------------------------------
% We wish to return the gradient w.r.t pn:
    dJ = (N1'*dJ)/lrnorm;

    end

%% Barrier

    function [bar,dbar] = barrier(pn)

% Full parameters        
% -------------------------------------------------------------------------        
        p = pstar + N*pn; % This is the full parameter vector
        
% However, this might be useful...
        Z = mss_v2s(Almi*p + blmi); % Note the tolerance

% Compute the gradient of -log det Z:
        Zinv = inv(Z);
        Zinv(1:nphi+1:nphi^2) = 0.5*Zinv(1:nphi+1:nphi^2);

        db = -Almi'*2*mss_s2v(Zinv);  
        
% Correct the gradient:
        db = N'*db;        
        
        bar = logdetbar(Z - lmitol*eye(nphi));
        dbar = db;
    
% Optional: trace constraint.
    if trcnst
    	[trb,trdb] = tracebarrier(pn);
    	bar = bar + trb;
        dbar = dbar + trdb;
    end
    
    end


%% log
    function l = logbar(a)
        
        if a <= 0 
            l = inf;
        else
            l = -log(a);
        end
    end

%% logdetbar
    function l = logdetbar(a)
        
        minev = min(eig(a));
        
        if minev <= 0     
            l = inf;      
        else            
            l = -log(det(a));            
        end
        
    end

%% lagrangian_hessian

    function [ddJ,lrhtimes] = lagrangianHessian(pn,varargin)

% NOTE: for debugging purposes, we will return a bunch of computation times        
      
% Input:
%       - pn, full parameter vector
%       - B, some approximation of the Jacobian of Delta(theta).
% NOTE: if B is not present then we will calculate the Jacobian of Delta
        
% First obtain p (which represents the model parameters):      
    p = (pstar(1:end-dimNL) + N1*pn);

% NEW: test code for normalized Hessian
    if ~isempty(varargin)
       p_tmp = p;
       p(1:ncoef(1)) = p(1:ncoef(1))/p_tmp(1); % Not the most sophisticated normalization
       p(ncoef(1)+1:sum(ncoef(1:2))) = p(ncoef(1)+1:sum(ncoef(1:2)))/p_tmp(1);        
    end
    
% This is the Hessian that we will populate:
    ddJ = zeros(nth);
    
% Construct block matrices: it seems like the loops are actually faster.   
    Fbar = Fbar_e{1}*p(1);
    
    lrhts1 = tic;
    
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

    lrht1 = toc(lrhts1);
    
% Solve for Delta*
    M = Gbar'*Gbar - (Fbar + Fbar');
    Delta = ((-(eps' + eta'*Gbar))/(M))';
    
% Jacobian of Delta
% -------------------------------------------------------------------------     
    
% For convenience,
    io = sum(ncoef(1:2));

% We construct dM one column at a time...each column corresponds to the
% derivative w.r.t. an element of theta
    dM = zeros(d1,nth);

    lrhts2 = tic;
    
% Handle e,    
    for i = 1:ncoef(1)
        dM(:,i) = -(Fbar_e{i} + Fbar_e{i}')*Delta;
    end
    
% Handle f,
    for i = 1:ncoef(2)
        dM(:,ncoef(1)+i) = -(Fbar_f{i} + Fbar_f{i}')*Delta;
    end
    
% Handle g,
    for i = 1:ncoef(3)
        tmp = Gbar'*Gbar_g{i};
        dM(:,io+i) = (tmp + tmp')*Delta;
    end    
    
% Now we construct dm.
% We begin by building dm/dg
    dm = [eps_e_full, eps_f_full, blks.m_g*(kron(speye(ncoef(3)),p(io+1:io+ncoef(3)))) + blks.m_gc];
    
    
    lrht2 = toc(lrhts2);
    
% Now we solve for the Jacobian

% Out of curiosity, how expensive is this?
    lrhts3 = tic;
    
    B = M\(-dM-dm); % dDelta
    
    lrht3 = toc(lrhts3);
      
% Hessian: Part TWO
% -------------------------------------------------------------------------    

% It is most convenient to break this up into parameters associated with
% e,f,g.

    lrhts4 = tic;

% For e,
    for i = 1:ncoef(1)      
        tmp = 1*(-Delta'*(Fbar_e{i} + Fbar_e{i}') + eps_e{i}'); % Gradient of dJ_i w.r.t. Delta
%         for j = 1:ncoef(1)
        for j = 1:sum(ncoef)
            ddJ(i,j) = tmp*B(:,j);
        end
%         ddJ(i,1:ncoef(1)) = tmp*B(:,1:ncoef(1));
    end

% For f,    
    for i = 1:ncoef(2)      
        tmp = 1*(-Delta'*(Fbar_f{i} + Fbar_f{i}') + eps_f{i}'); % Gradient of dJ_i w.r.t. Delta:
%         for j = 1:ncoef(2)
%             ddJ(ncoef(1)+i,ncoef(1)+j) = tmp*B(:,ncoef(1)+j);
        for j = 1:sum(ncoef)
            ddJ(ncoef(1)+i,j) = tmp*B(:,j);
        end
    end
    
% For g,    
    for i = 1:ncoef(3)
        tmp = Gbar'*Gbar_g{i}; 
        tmp = tmp + tmp';
        tmp = 1*(Delta'*tmp + eta'*Gbar_g{i} + eta_g{i}'*Gbar);        
%         for j = 1:ncoef(3)
%             ddJ(io+i,io+j) = tmp*B(:,io+j);    
        for j = 1:sum(ncoef)
            ddJ(io+i,j) = tmp*B(:,j);
        end
    end

% Make sure the matrix is symmetric (notice, we have left out a 2 in all of 
% the above expressions, so we don't have to divide by 2 again).
    ddJ = ddJ + ddJ'; 
    
    lrht4 = toc(lrhts4);
    
% Hessian: Part ONE
% -------------------------------------------------------------------------

    lrhts5 = tic;

% This is a big matrix: [...G_i*Delta+eta_i...] (as a columns)
   GiDeta = reshape(Gbar_g_full*Delta,d2,ncoef(3)) + eta_g_full;
   
   for i = 1:ncoef(3)
       for j = 1:i-1
           tmp = 2*GiDeta(:,i)'*GiDeta(:,j);
%            tmp = 2*(Gbar_g{i}*Delta+eta_g{i})'*(Gbar_g{j}*Delta+eta_g{j});
           ddJ(io+i,io+j) = ddJ(io+i,io+j) + tmp;
           ddJ(io+j,io+i) = ddJ(io+j,io+i) + tmp;
       end
       tmp = 2*GiDeta(:,i)'*GiDeta(:,i);
%        tmp = 2*(Gbar_g{i}*Delta+eta_g{i})'*(Gbar_g{i}*Delta+eta_g{i});
       ddJ(io+i,io+i) = ddJ(io+i,io+i) + tmp;
   end    
   
   lrht5 = toc(lrhts5);
   
% Store all the times
    lrhtimes = [lrht1;lrht2;lrht3;lrht4;lrht5];
    
% Corrected Hessian output
    dimx = dim_pmod - nth;
    
    ddJ = N1'*[ddJ, zeros(nth,dimx); zeros(dimx,nth+dimx)]*N1; 
    
    
    end


%% Compute Delta

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

%% FDA: Delta

    function dDelta = fdaDelta(pn)
      
        fdd = 1e-4;        
        delta0 = compDelta(pn);
        dDelta = zeros(length(delta0),nth);
        
        for i = 1:nth
           
            tmp = pn;
            tmp(i) = tmp(i) + fdd; 
            delta1 = compDelta(tmp);
            dDelta(:,i) = (delta1-delta0)/fdd;
            
        end
           
    end

%% FDA: LR GRADIENT

% Finite difference approx. of LR
    function dJ = fdal(pn,fdd)

% Most naive FDA possible
        dJ = zeros(length(pn),1);
        J0 = lagrangian(pn);
        for i = 1:length(pn)
            tmp = pn;
            tmp(i) = tmp(i) + fdd;
            J1 = lagrangian(tmp);
            dJ(i) = (J1-J0)/fdd;
        end
        
    end


%% FDA: LR HESSIAN

% Finite difference approximation of Hessian

    function hess = fdalh(pn,fdd)

        hess = zeros(nth);

        [~,dJ0] = lagrangian(pn);
        dJ0 = dJ0(1:nth);

        for pindex = 1:nth
            tmp = pn;
            tmp(pindex) = tmp(pindex) + fdd;
            [~,dJ1] = lagrangian(tmp);

            hess(:,pindex) = (dJ1(1:nth) - dJ0)/fdd;

        end

        % Perhaps we should symmetrize:
        hess = (hess + hess')/2;

%         hess = hess - 1.1*min([0,min(eig(hess))])*eye(nth); % PD 
    end


%% FDA: barrier GRADIENT

% Finite difference approx. of barrier
    function db = fdab(pn,fdd)

% Most naive FDA possible
        db = zeros(length(pn),1);
        b0 = barrier(pn);
        for i = 1:length(pn)
            tmp = pn;
            tmp(i) = tmp(i) + fdd;
            b1 = barrier(tmp);
            db(i) = (b1-b0)/fdd;
        end
        
    end

%% FDA: barrier HESSIAN

% Finite difference approx. of barrier
    function db = fdabh(pn,fdd)

% Most naive FDA possible
        db = zeros(length(pn));
        [~,b0] = barrier(pn);
        for i = 1:length(pn)
            tmp = pn;
            tmp(i) = tmp(i) + fdd;
            [~,b1] = barrier(tmp);
            db(:,i) = (b1-b0)/fdd;
        end
        
    end

%% EXPERIMENTAL: trace barrier and gradient

    function [trb,trdb] = tracebarrier(pn)
       
        p = pstar + N*pn; % This is the full parameter vector
%         p = (pstar(1:end-dimNL) + N1*pn);
%         p = p(1:sum(ncoef));
        
        tmp = trlim - Atr*p;
        
        trb = logbar(tmp);
        trdb = -(1/tmp)*Atr';
        
    end

%% EXP: trace Hessian

    function trh = tracehessian(pn)
       
        p = pstar + N*pn; % This is the full parameter vector
%         p = (pstar(1:end-dimNL) + N1*pn);
%         p = p(1:sum(ncoef));
        
        tmp = trlim - Atr*p;
        
        trh = (1/tmp^2)*(Atr'*Atr);
        
    end

end







