function res = lr_init_lmiPLusLin(res_sos)
% Use SDP to find a valid initial model.
% Should probably just find the analytic centre of the constraints, but for
% now...

% Basically, we just want to choose a parameter vector `theta' such that
% Z = A_lmi*theta + b_lmi >= 0 
% A_eq*theta = b_eq;

%% Basic preprocessing

Almi = res_sos.Almi;
blmi = res_sos.blmi;
Aeq = res_sos.Aeq;
beq = res_sos.beq;
tol = res_sos.lmitol;
nphi = res_sos.nphi;


%% Yalmip implementation

yalmip('clear')

% Parameter vector:
p = sdpvar(size(Almi,2),1);
d = sdpvar(1);
% d = 1e-3;

% LMI constraint:
const = mss_v2s(Almi*p + blmi) >= d*eye(nphi);
const = [const, d >= tol];

% Linear constraint:
const = [const, Aeq*p == beq];

% Experimental trace constraint
const = [const, res_sos.trlim - res_sos.Atr*p >= 0];

% Solve:
options = sdpsettings('solver','mosek');
options.verbose = 0;
sol = solvesdp(const,-d,options);
% sol = solvesdp(const,p'*p,options);


% Some parameters will be undefined:
p = double(p);
p(isnan(p)) = 0;

%% Output results

res.p = p;
res.sol = sol;
res.tol = double(d);

yalmip('clear')

end

