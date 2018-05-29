function res = ee_nonlin(nlid_data,varargin)
% Equation error with well-posedness condition.

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
    verbose = 1;
end

if isfield(options,'stability')    
    stability = options.stability;
else
    stability = 0;
end

% Quick hack
if isfield(options,'linear')    
    linear = options.linear;
else
    linear = 0;
end

%% Problem formulation

yalmip('clear')

% Define indeterminate for state: the term in quotes ('x') gives the name 
% of the variables as they will appear when 'printed'.
x = sdpvar(nx,1);
u = sdpvar(nu,1);

ev = e_monos(x);
fv = f_monos(x,u);
% gv = g_monos(x);  
gv = g_monos(x,u);
    
gc = sdpvar(size(gv,2),1);
fc = sdpvar(size(fv,2),1);
ec = sdpvar(size(ev,2),1);

% Define the system (i.e. polynomials)

e = ev*ec;
f = fv*fc;
g = gv*gc;

E = jacobian(e,x);

fprintf('Starting problem setup...')
t1 = tic;

% Equation error: output
eta = 0;
for t = 1:T
    tmp = yd(:,t) - replace(g,[x;u],[xd(:,t);ud(:,t)]);
    eta = eta + tmp'*tmp;
end

% Equation error: state
eps = 0;
for t = 1:T-1
    tmp = replace(e,x,xd(:,t+1)) - replace(f,[x;u],[xd(:,t);ud(:,t)]);
    eps = eps + tmp'*tmp;
end

t2 = toc(t1);
fprintf('Done (%.5e sec)\n',t2)

%% Stability constraints

if stability
    E = jacobian(e,x);
    F = jacobian(f,x);
    G = jacobian(g,x);
    
    P = sdpvar(nx);

    H = [E + E' - P, F', G';
         F, P, zeros(nx,ny);
         G, zeros(ny,nx), eye(ny)];

    [sH,~] = size(H); 
    dum = sdpvar(sH,1); % Dummy variable for matrix SOS

     stbl = sos(dum'*H*dum);
     constraints = stbl;
     sosparam = [ec;fc;gc;mss_s2v(P)];
else    
    
    if linear        
        constraints = E + E' >= 1e-6*eye(nx);  
        
        options = sdpsettings('solver','mosek');
        options.verbose = verbose;   
        sol = solvesdp(constraints, eta + eps,options);        
        
    else
        % Well-posedness constraint
        constraints = sos(E+E');
        % sosparam = [ec;fc;K;gc];
        sosparam = [ec;fc;gc]; 
        
        options = sdpsettings('solver','mosek');
        options.verbose = verbose;   
        sol = solvesos(constraints, eta + eps,options,sosparam);        
        
    end
end
    
%% Solve the problem
% options = sdpsettings('solver','mosek');
% options.verbose = verbose;   
% sol = solvesos(constraints, eta + eps,options,sosparam);

ec_id = double(ec);
fc_id = double(fc);
gc_ee = double(gc);

ec_ee = ec_id/ec_id(1);
fc_ee = fc_id/ec_id(1);

% fprintf('\n------------------------------------------------\n')
fprintf('\n')
fprintf('EE results\n')
fprintf('------------------------------------------------\n')

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

yalmip('clear')

end

