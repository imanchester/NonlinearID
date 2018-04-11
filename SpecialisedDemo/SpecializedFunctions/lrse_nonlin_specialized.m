function res = lrse_nonlin_specialized(nlid_data,varargin)
% Identify a dynamical system via Lagrangian relaxation with a specialized
% interior point algorithm.


%% Initial Processing

% Extract data:
xd = nlid_data.x;
yd = nlid_data.y;
ud = nlid_data.u;

% Dimensions:
[nx,~] = size(xd);
[ny,~] = size(yd);
[nu,~] = size(ud);

% Monomials:
e_monos = nlid_data.e_monos;
f_monos = nlid_data.f_monos;
g_monos = nlid_data.g_monos;


%% Process varargin

if ~isempty(varargin)
    options = varargin{1};
else
    options = [];
end

% Different levels
% 0 - absolutely nothing
% 1 - essentials (probably most useful)
% 2 - full (everything you might want)

if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 1;
end

warning('on','MATLAB:nearlySingularMatrix')
if verbose <= 1
    warning('off','MATLAB:nearlySingularMatrix')
end

if isfield(options,'recordse')
    recse = options.recordse; % Record simulation error
else
    recse = 0;
end

if isfield(options,'trcnst')
    trcnst = options.trcnst; % Impose trace constraint
else
    trcnst = 0;
end

if isfield(options,'earlyterm')
    earlyterm = options.earlyterm; % Terminate inner iterations based on actual objective
else
    earlyterm = 0;
end


%% Problem set-up makes use of SPOT functions.
% -------------------------------------------

clear pr
pr = spotsosprog;

% NEW set-up:
% -----------

% Indeterminate variables:
[pr,x] = pr.newIndeterminate('x',nx);
[pr,u] = pr.newIndeterminate('u',nu);

% Now we define polynomials for each of these coefficients
e_pol = e_monos(x);
f_pol = f_monos(x,u);
g_pol = g_monos(x,u); % CHANGED: JU2911

% Coefficients
[pr,ce] = pr.newFree(size(e_pol,2)); 
[pr,cf] = pr.newFree(size(f_pol,2)); 
[pr,cg] = pr.newFree(size(g_pol,2)); 

th = [ce;cf;cg];
nth = length(th);

% Store the number of coefficients:
ncoef(1) = length(ce);
ncoef(2) = length(cf);
ncoef(3) = length(cg);

% Build e,f,g
e = e_pol*ce;
f = f_pol*cf;
g = g_pol*cg;

% If we want to compute the simulation error...
if recse
    e_fnc = @(x,c) e_monos(x)*c;
    f_fnc = @(x,u,c) f_monos(x,u)*c;
    g_fnc = @(x,u,c) g_monos(x,u)*c; % CHANGED: JU2911
end

if verbose >= 2; fprintf('Number of model parameters: %d\n',sum(ncoef)); end;


% Encode SOS constraint
% ---------------------
if verbose; fprintf('\tEncoding SOS constraint:\n'); end;
t1 = tic;
% res_sos = lr_sos2lmiPlusLin(pr,e,f,g,x,th);
ops_sos = [];
ops_sos.verbose = verbose;
res_sos = lr_sos2lmiPlusLin_2(pr,e,f,g,x,th,ops_sos);
if verbose; fprintf(' done (%.2f sec).\n',toc(t1)); end;

% Generate block matrices
% -----------------------
if verbose; fprintf('\tBuilding block matrices:'); end;
t1 = tic;
% res_blk = lr_genblks(e_pol,f_pol,g_pol,x,u,xd,ud,yd);
% res_blk = lr_genblks_2(e_pol,f_pol,g_pol,x,u,xd,ud,yd);

if isfield(options,'blks')
    res_blk = options.blks;
else
    res_blk = lr_genblks_4(e_pol,f_pol,g_pol,x,u,xd,ud,yd);
end
res_blk.lrnorm = 1; % We don't want to normalize...for now...

if verbose; fprintf(' done (%.2f sec).\n',toc(t1)); end;

% Initialization
% --------------

if verbose; fprintf('\tFinding initial model:'); end;
t1 = tic;
res_init = lr_init_lmiPLusLin(res_sos); % used this for the paper
% res_init = lr_init_lmiPLusLin_2(res_sos);
% res_init = lr_init_simple_lmiPLusLin(pr,e,x,th,res_sos);
if verbose; fprintf(' %s (%.2f sec).\n',res_init.sol.info,toc(t1)); end;

if (res_init.sol.problem ~= 0) || (res_init.tol < 1e-3) 
    fprintf('\n\nAborting: sorry, we were unable to find a feasible (or well conditioned) initial model.\n')
    fprintf('Please consider revising your monomial selection.\n')
    fprintf('e.g., try making sure the maximum degree of the polynomial in each component of e is the same.\n')
    error('Unable to find a feasible (or well conditioned) initial model.')
end

% Original
% -------------------------------------------------------------------------
% NB: This theta_star (from the init function) is:
% theta = [p;Q], i.e. model parameters and Gram matrix.
% What we will actually search for, shall be a vector that spans the
% nullspace of A (in the constraint Ax+b). 

theta_star = res_init.p;
theta = zeros(size(res_sos.N,2),1);
theta_prev = theta;

% max(theta_star)
% res_init.tol
% 
% dum = 1;



% Tau loop
% ----------

% Note: theta contains everything that we must search for (model
% parameters, stability certificate and additional parameters for sos
% constraint).

% How about printing the actual number of decision variables:
ndecvar = length(theta);
if verbose >= 2; fprintf('\tNumber of decision variables: %d.\n',ndecvar); end;

tau = 1e4;
iters = 50;
Js = zeros(iters,1); % value of bound
ftaus = zeros(iters,1); % f_tau
J_prev = inf;
bailout = 0;

np = length(theta);
B_i = 0; % Code: compute B_1 inside function


% NEW (23/11/2016)
bfgs_options.verbose = verbose;
bfgs_options.earlyterm = earlyterm;
bfgs_options.trcnst = trcnst;

if verbose; fprintf('\nBeginning tau loop...\n\n'); end;

if recse
    ses = zeros(1,iters);
    se_best = inf;
    theta_best = 0;
end

%% Start the tau loop

% Before we begin, we will store some stuff
res.ndecvar = ndecvar;
res.blks = res_blk;
res.sos = res_sos; 
res.theta_star = theta_star;
res.blks = res_blk; % Pass some preprocessing out for reuse...

% Record the delta* after each iteration
deltas = zeros(length(xd)*nx,iters);

ts = tic;
for i = 1:iters
     
    if verbose >= 2; fprintf('\nIter: %d.\n',i); end;
    
% Optional: don't terminate early on first iteration
% -------------------------------------------------------------------------
    if i == 1
        bfgs_options.earlyterm = 0;
    else
        bfgs_options.earlyterm = earlyterm;
    end  
    
% Run inner iterations
% -------------------------------------------------------------------------

    res_bfgs = lrse_minftau(xd,yd,ud,theta,theta_star,tau,B_i,res_sos,res_blk,bfgs_options);

% Update the parameters:
    theta = res_bfgs.p;
    
% Update the Hessian:    
    B_i = res_bfgs.B; % Standard method
%     B_i = fdh(1e-4,theta,theta_star,res_sos,res_blk);

% Store the Newton times:
    newtontimes{i} = res_bfgs.times;
    
% Store lne search iterations
    lsiters{i} = res_bfgs.lsiters;

    hessiantimes{i} = res_bfgs.hessiantimes;
    searchtimes{i} = res_bfgs.searchtimes;
        
% For interest: also store the Delta, to see how it evolves
% I don't really think this is necessary anymore, right? It wasn't that
% interesting.
    res.theta_ld = theta; %
    deltas(:,i) = computeDelta(res,yd);
     
            
    J = lr_bound_eqcSOS(theta,theta_star,res_blk,res_sos);
    
    if verbose >= 2
        fprintf('Objective: %.5e\n',J*res_blk.lrnorm); 
        fprintf('Tau: %.5e\n',tau);
    end
    Js(i) = J;
    ftaus(i) = res_bfgs.ftau; % Now store the f_tau for itnerest
        
    tau = tau/5; % Update tau
    

% Optional: check simulation error at the conclusion of each tauiteration. 
    if recse
        theta_tmp = theta_star + res_sos.N*theta;
        ec_bfgs_new_u = theta_tmp(1:ncoef(1));
        ec_bfgs_tmp = ec_bfgs_new_u/ec_bfgs_new_u(1);
        fc_bfgs_tmp = theta_tmp(ncoef(1)+1:sum(ncoef(1:2)))/ec_bfgs_new_u(1);
        gc_bfgs_tmp = theta_tmp(sum(ncoef(1:2))+1:sum(ncoef(1:3)));

        [se_tmp,~] = se_nonlinear_2(ud,yd,xd(:,1),e_fnc,f_fnc,g_fnc,ec_bfgs_tmp,fc_bfgs_tmp,gc_bfgs_tmp); % JU2911  
        
        ses(i) = se_tmp; % Store sim error
        
        if verbose; fprintf('Sim error: %.5e\n',se_tmp); end
        
        if se_tmp < se_best
            se_best = se_tmp;
            res.ec_best = ec_bfgs_tmp;
            res.fc_best = fc_bfgs_tmp;
            res.gc_best = gc_bfgs_tmp;
        end
        
    end

    if bailout
        i
        break
    end
    
    if tau < 1e-15
        if verbose >= 2; fprintf('Bailing out. tau = %.2e, i = %d\n',tau,i); end;
        break;
    end

    if abs(J - J_prev) < 1e-10
        if verbose >= 2; fprintf('Breaking out (no improv). tau = %.2e, i = %d\n',tau,i); end;
        break;
    end
    
    if J < 0 
        theta = theta_prev;
        J = J_prev;
        break
    end
    
    theta_prev = theta;
    J_prev = J;    
end

t_lr_custom = toc(ts);

if verbose >= 2
figure
semilogy(1:iters,Js)
title('Convergence: Quasi-Newton')
xlabel('Iteration')
ylabel('J_{\lambda}')

figure
semilogy(1:iters,abs(ftaus))
title('Convergence: Quasi-Newton')
xlabel('Iteration')
ylabel('f_{\tau}')
end

theta_ld = theta;

% Recover the coefficients:
theta = theta_star + res_sos.N*theta;
ec_bfgs_new_u = theta(1:ncoef(1));
fc_bfgs_new_u = theta(ncoef(1)+1:sum(ncoef(1:2)));
gc_bfgs_new = theta(sum(ncoef(1:2))+1:sum(ncoef(1:3)));

ec_bfgs_new = ec_bfgs_new_u/ec_bfgs_new_u(1);
fc_bfgs_new = fc_bfgs_new_u/ec_bfgs_new_u(1);

if verbose
% fprintf('\n------------------------------------------------\n')
fprintf('\n')
fprintf('LRSE Specialized results \n')
fprintf('------------------------------------------------\n')
fprintf('\tFinal bound: %.8e \n',J*res_blk.lrnorm)
fprintf('\tLoop time: %.3f sec\n',t_lr_custom)
fprintf('\tTau iterations: %d\n\n',i)
end


%% Output results

res.ec = ec_bfgs_new;
res.fc = fc_bfgs_new;
res.gc = gc_bfgs_new;

res.theta = theta; % This is the full theta

res.Js = Js(1:i);
res.Jf = J;

res.solvertime = t_lr_custom;
res.ndecvar = ndecvar;

% Some extra stuff that can be useful for continuing with a cutting plane
% method:
res.theta_ld = theta_ld; % This is the lower dimensional theta

% Return the times for each Newton iteration, for each tau iteration
res.newtontimes = newtontimes;

% Return no. of line search iterations (for each inner iteration)
res.lsiters = lsiters;

res.hessiantimes = hessiantimes;
res.searchtimes = searchtimes;

% Store the delta*'s after each iteration,
res.deltas = deltas(:,1:i);

if recse
    res.ses = ses(1:i);
end



end

