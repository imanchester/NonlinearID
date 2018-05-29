function res = lree_sos_nonlinear(nlid_data,varargin)
% LR of equation error, but with sos instead of linearization

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

% if isfield(options,'stability')    
%     stability = options.stability;
% else
%     stability = 0;
% end

if isfield(options,'stability')   
    
    if strcmp(options.stability,'none')
        stability = 0;
    elseif strcmp(options.stability,'state')
        stability = 1;
    elseif strcmp(options.stability,'output')
        stability = 2;
    end
    
%     stability = options.stability;
else
    stability = 1; % Default is state-stability
end


if isfield(options,'multi')    
    multi = options.multi;
else
    multi = 1e-4;
end

% if isfield(options,'stability')    
%     stability = options.stability;
% else
%     stability = 0;
% end

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
% res_sos = lr_sos2Axb(pr,e,f,g,x,th);
res_sos = lr_sos2lmiPlusLin_2(pr,e,f,g,x,th);

% Now that's done we can define some decision variables:
% res_init = lr_init_eqcSOS(res_sos);
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
% qc = theta(end-res_sos.nQ+1:end);

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

constraints = []; % Old notation

end

%% Define the system (i.e. polynomials)

dummy = 1;

fprintf('Passed break point\n')

e = ev*ec;
f = fv*fc;
g = gv*gc;

E = jacobian(e,x);
F = jacobian(f,x); % Might need these for stability... 
G = jacobian(g,x);


%% Very basic implementation of LREE

fprintf('Beginning LREE problem set-up...')
ts = tic;

% In the past, I think minimization of equation error was a little
% unreliable in Yalmip? Seemed to get different results from different
% formulations...


% 1. So, first we will try the equation error formulation from the fast ee
% function:

% Equation error: output
by = 0;
ay = 0;
for t = 1:T
%     tmp = yd(:,t) - replace(g,[x;u],[xd(:,t);ud(:,t)]);
    my = g_monos(xd(:,t),ud(:,t));
    
    by = by - 2*yd(:,t)'*my;
    ay = ay + my'*my;
    
end

% by = by/T;
% ay = ay/T;

eta = (yd(:)'*yd(:)) + by*gc + gc'*ay*gc;


% 2. Then we add LR of state ee:

% Choose a `multiplier',
% Q = 1e-4*eye(nx);
Q = multi*eye(nx);

% Slack variables
s = sdpvar(T-1,1);

for t = 1:T-1                     

% % Equation errors:    
% %     e_x = replace(f,[x;u],[xd(:,t);ud(:,t)]) - replace(e,x,xd(:,t+1));
% %     e_y = replace(g,[x;u],[xd(:,t);ud(:,t)]) - yd(:,t);
% % 
% % Jacobians
% %     E_x = replace(E,x,xd(:,t));
% % 
% % LREE
% %     tmp = [s(t), e_x'*Q; Q'*e_x, Q'*E_x + E_x'*Q - eye(nx)] >= 0;
    
% SOS constraint
    tmp = s(t) - ( (x-xd(:,t+1))'*(x-xd(:,t+1)) - 2*(x-xd(:,t+1))'*(e - replace(f,[x;u],[xd(:,t);ud(:,t)])) );
            
    constraints = [constraints, sos(tmp)];

end

objective = sum(s) + eta;
fprintf('done (%.2f sec)\n',toc(ts))



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
             
% Storage function,
%     P = sdpvar(nx);
%     c_p = P >= del*eye(nx); 
%     constraints = [constraints, c_p];
    
% What if we just put in the sos constraint?
%     Hsos = [E + E' - P - 1e-9*eye(nx), F', G';
%             F,P,zeros(nx,ny);
%             G,zeros(ny,nx),eye(ny)];
%     dum = sdpvar(2*nx+ny,1);



if verbose
%     fprintf('Finished with setup, passing to MOSEK\n')
end

% Perform optimization
% options = sdpsettings('solver','mosek');
% options.verbose = verbose;
% sol = solvesdp(constraints,objective,options);

%%

ec_id = double(ec);
fc_id = double(fc);
gc_ee = double(gc);

ec_ee = ec_id/ec_id(1);
fc_ee = fc_id/ec_id(1);

% fprintf('\nTrue e : Identified e:\n')
% % sdisplay(ev*ec_nrm)
% [ec_tr, ec_ee]
% 
% fprintf('\nTrue f : Identified f:\n')
% % sdisplay(fv*fc_nrm)
% [fc_tr, fc_ee]
% 
% fprintf('\nTrue g : Identified g:\n')
% % sdisplay(gv*gc_id)
% [gc_tr, gc_ee]

res.ec = ec_ee;
res.fc = fc_ee;
res.gc = gc_ee;

res.solvertime = sol.solvertime;
res.sol = sol;
res.objective = double(objective);


yalmip('clear')

end

