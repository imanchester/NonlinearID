function res = lr_init_lmiPLusLin_spot(res_sos)
% Use SDP to find a valid initial model.

% We want to choose a parameter vector `theta' such that
% Z = A_lmi*theta + b_lmi >= 0 
% A_eq*theta = b_eq;

%% Basic preprocessing

Almi = res_sos.Almi;
blmi = res_sos.blmi;
Aeq = res_sos.Aeq;
beq = res_sos.beq;
tol = res_sos.lmitol;
nphi = res_sos.nphi;

%%

% fprintf('Using spot\n')

pr = spotprog;

[pr,d] = pr.newPos(1,1);

[pr,p] = pr.newFree(size(Almi,2),1);

pr = pr.withPSD(mss_v2s(Almi*p + blmi) - d*eye(nphi));

pr = pr.withEqs(Aeq*p - beq);

pr = pr.withPos(res_sos.trlim - res_sos.Atr*p);

opt = spot_sdp_default_options();
opt.verbose = 0;
sol_spot = pr.minimize(-d, @spot_mosek, opt); 
% sol_spot = pr.minimize(-d, @spot_sedumi, opt); 


d_spot = double(sol_spot.eval(d));

p_spot = double(sol_spot.eval(p));

clear sol

if d_spot >= 1e-3
    sol.info = 'Feasible.';
    sol.problem = 0;
else
    sol.info = 'Problematic.';
    sol.problem = 1;
end

%% Output results

res.p = p_spot;
res.tol = d_spot;
res.sol = sol;

yalmip('clear')

end

