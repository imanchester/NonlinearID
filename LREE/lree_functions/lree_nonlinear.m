function res = lree_nonlinear(nlid_data,varargin)
% LR of equation error (LREE)

% Input arguments:
% - nlid_data: a structure with the following fields,
%   .x - [nx,T] array of states
%   .y - [ny,T] array of outputs
%   .u - [nu,T] array of inputs
%   .e_monos - a function handle, mapping x to an [nx,ne] array of ne
%   monomials, such that e(x) = e_monos(x)*ec, where ec is a vector of
%   coefficients
%   .f_monos - a function handle, mapping (x,u) to an [nx,nf] array of nf
%   monomials, such that f(x,u) = f_monos(x,u)*fc, where fc is a vector of
%   coefficients
%   .g_monos - a function handle, mapping (x,u) to an [ny,ng] array of ng
%   monomials, such that g(x,u) = g_monos(x,u)*gc, where gc is a vector of
%   coefficients

% Outputs:
% - res: a structure containing the following fields,
%   .ec - vector of coefficients such that e(x) = e_monos(x)*ec
%   .fc - vector of coefficients such that f(x,u) = f_monos(x,u)*fc
%   .gc - vector of coefficients such that g(x,u) = g_monos(x)*gc
%   .solvertime - computation time (in SDP solver)
%   .sol - information from SDP parser/solver
%   .objective - value of the bound on equation error 



%% Initial Processing

% Extract data:
xd = nlid_data.x;
yd = nlid_data.y;
ud = nlid_data.u;

% Dimensions:
[nx,~] = size(xd);
[ny,T] = size(yd);
[nu,~] = size(ud);

% Monomials:
e_monos = nlid_data.e_monos;
f_monos = nlid_data.f_monos;
g_monos = nlid_data.g_monos;

%% Options

if length(varargin) >= 1
    options = varargin{1};
else
    options = [];
end

if isfield(options,'verbose')    
    verbose = options.verbose;
else
    verbose = 0;
end

if isfield(options,'stability')   
    
    if strcmp(options.stability,'none')
        stability = 0;
    elseif strcmp(options.stability,'state')
        stability = 1;
    elseif strcmp(options.stability,'output')
        stability = 2;
    end
    
else
    stability = 1; % Default is state-stability
end


if isfield(options,'multi')    
    multi = options.multi;
else
    multi = 1e-4;
end


%% Problem formulation

del = 1e-5;
yalmip('clear')

%% Decision variables

if stability == 2
    
% Custom SOS: search over a parametrization of all models that satisfy the 
% linear equality constraints.

% For convenience, we use SPOT in exactly the same way as in the
% specialized algorithm.
clear pr
pr = spotsosprog;

% Indeterminate variables:
[pr,x] = pr.newIndeterminate('x',nx);
[pr,u] = pr.newIndeterminate('u',nu);

% Now we define polynomials for each of these coefficients
e_pol = e_monos(x);
f_pol = f_monos(x,u);
g_pol = g_monos(x,u);

% Coefficients
[pr,ce] = pr.newFree(size(e_pol,2)); 
[pr,cf] = pr.newFree(size(f_pol,2)); 
[pr,cg] = pr.newFree(size(g_pol,2)); 

th = [ce;cf;cg];

% Build e,f,g
e = e_pol*ce;
f = f_pol*cf;
g = g_pol*cg;

% Encode SOS constraint
res_sos = lr_sos2lmiPlusLin_2(pr,e,f,g,x,th);

% Now that's done we can define some decision variables:
res_init = lr_init_lmiPLusLin(res_sos);

thetabar = sdpvar(size(res_sos.N,2),1); % These are the actual SDP var decision variables

theta = res_init.p + res_sos.N*thetabar;

% Extract the coefficients:
nec = size(e_pol,2);
nfc = size(f_pol,2);
ngc = size(g_pol,2);
npc = (nx+1)*nx/2;

ec = theta(1:nec);
fc = theta(nec+1:nec+nfc);
gc = theta(nec+nfc+1:nec+nfc+ngc);
pc = theta(nec+nfc+ngc+1:nec+nfc+ngc+npc);

clear pr x u

x = sdpvar(nx,1);
u = sdpvar(nu,1);

ev = e_monos(x);
fv = f_monos(x,u);
gv = g_monos(x,u);

% Storage function,
P = mss_v2s(pc);
constraints = [];

else

x = sdpvar(nx,1);
u = sdpvar(nu,1);

ev = e_monos(x);
fv = f_monos(x,u);
gv = g_monos(x,u);    
    
gc = sdpvar(size(gv,2),1);
fc = sdpvar(size(fv,2),1);
ec = sdpvar(size(ev,2),1);

constraints = [];

end

%% Define the system (i.e. polynomials)

% Functions
e = ev*ec;
f = fv*fc;
g = gv*gc;

% Jacobians
E = jacobian(e,x);
F = jacobian(f,x); 
G = jacobian(g,x);


%% Very basic implementation of LREE

% 1. Least squares fit of output map (g)
if verbose; fprintf('Beginning LREE problem set-up...'); end
ts = tic;

% Equation error: output
by = 0;
ay = 0;
for t = 1:T
    my = g_monos(xd(:,t),ud(:,t));
    
    by = by - 2*yd(:,t)'*my;
    ay = ay + my'*my;
    
end

eta = (yd(:)'*yd(:)) + by*gc + gc'*ay*gc;


% 2. LR of equation error:

% Choose a multiplier,
Q = multi*eye(nx);

% Slack variables
s = sdpvar(T-1,1);

for t = 1:T-1                     

% Equation errors:    
    e_x = replace(f,[x;u],[xd(:,t);ud(:,t)]) - replace(e,x,xd(:,t+1));

% Jacobians
    E_x = replace(E,x,xd(:,t+1)); 

% LREE
    tmp = [s(t), e_x'*Q; Q'*e_x, Q'*E_x + E_x'*Q - eye(nx)] >= 0;
    
    constraints = [constraints, tmp];

end

objective = sum(s) + eta;


%% Add stability constraint

if verbose; fprintf('done (%.2f sec)\n',toc(ts)); end

if stability == 0 % No stability, just well-posedness
    
    Hsos = E + E' - 1e-8*eye(nx);
    dum = sdpvar(nx,1);   
    
    constraints = [constraints, sos(dum'*Hsos*dum)]; 
    
    ops = sdpsettings('solver','mosek');
    ops.verbose = verbose;
    sol = solvesos(constraints,objective,ops,[ec;fc;gc]);  
    
    
elseif stability == 1
    
    P = sdpvar(nx);
    c_p = P >= del*eye(nx); 
    constraints = [constraints, c_p]; 
    
    Hsos = [E + E' - P - 1e-9*eye(nx), F';
            F,P];   
    dum = sdpvar(2*nx,1);
    
    constraints = [constraints, sos(dum'*Hsos*dum)]; 
    
    ops = sdpsettings('solver','mosek');
    ops.verbose = verbose;
    sol = solvesos(constraints,objective,ops,[ec;fc;gc;P(:)]);     
    
    
elseif stability == 2
    
    Z = mss_v2s(res_sos.Almi*theta + res_sos.blmi);
    constraints = [constraints, Z >= 0];
    
    options = sdpsettings('solver','mosek');
    options.verbose = verbose;
    sol = solvesdp(constraints,objective,options); 
    
end
             


%% output results

ec_id = double(ec);
fc_id = double(fc);
gc_ee = double(gc);

ec_ee = ec_id/ec_id(1);
fc_ee = fc_id/ec_id(1);

% Identified coefficients
res.ec = ec_ee;
res.fc = fc_ee;
res.gc = gc_ee;

res.solvertime = sol.solvertime;
res.sol = sol;
res.objective = double(objective);

yalmip('clear')

end

