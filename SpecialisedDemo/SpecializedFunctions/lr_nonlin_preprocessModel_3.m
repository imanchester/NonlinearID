function res = lr_nonlin_preprocessModel_3(pr,x,u,spot_monos_e,spot_monos_f,spot_monos_g,nx,nu,ny)

% 29/11/2016 - the main purpose of this function is to allow direct
% feedthrough of control inputs (i.e. inputs in the output map).

% The idea with this script is: monomials in, monomials out - keep all the 
% anonymous functions on the outside

%% Problem set-up makes use of SPOT functions.
% -------------------------------------------

xu2v = map_xu2v(x,u);

[epow,eco,nem] = mono2mat(spot_monos_e,x);
[fpow,fco,nfm] = fg_mono2mat(spot_monos_f,xu2v,x,u); %JU2911
[gpow,gco,ngm] = fg_mono2mat(spot_monos_g,xu2v,x,u); %JU2911

% Build anonymous functions
% -------------------------------------------------------------------------

e_monos = @(x) gen_e_monos(x,epow,eco,nx,nem);
f_monos = @(x,u) gen_e_monos(xu2v*[x;u],fpow,fco,nx,nfm);
g_monos = @(x,u) gen_e_monos(xu2v*[x;u],gpow,gco,ny,ngm); % JU2911

% Now we define polynomials for each of these coefficients
e_pol = e_monos(x);
f_pol = f_monos(x,u);
g_pol = g_monos(x,u);

% Coefficients
[pr,ce] = pr.newFree(size(e_pol,2)); 
[pr,cf] = pr.newFree(size(f_pol,2)); 
[pr,cg] = pr.newFree(size(g_pol,2)); 

th = [ce;cf;cg];
% nth = length(th);

% Store the number of coefficients:
ncoef(1) = length(ce);
ncoef(2) = length(cf);
ncoef(3) = length(cg);

% Build e,f,g
e = e_pol*ce;
f = f_pol*cf;
g = g_pol*cg;

% This is straight from the usual LMI encoding function
% -------------------------------------------------------------------------

%% Stability constraint

% First step is to generate the polynomial that encodes the stability
% constraint

% We'll introduce P inside this function:
[pr,P] = pr.newSym(nx);

% Augment the parameter vector with P.
theta = [th;mss_s2v(P)];
nth = length(theta);

% Jacobians:
E = diff(e,x);
F = diff(f,x);
G = diff(g,x);

% Matrix that must be nonnegative for all x,u
H = [E + E' - P, F', G'; F, P, zeros(nx,ny) ;G, zeros(ny,nx), eye(ny)];

% Dummy varaibles:
[pr,z] = pr.newIndeterminate('z',2*nx+ny);

% This polynomial must be SOS:
p = z'*H*z;

%% Monomials

% Use SPOT to generate candidate monomials.
opt = [];
decvar = pr.variables;
phi = spotsosprog.buildGramBasis(p,decvar,opt);

%% Remove redundant monomials from p

% So: some of the coefficients of monomials in p must be zero.
% Therefore these monomials do not appear in the expansion phi'*Q*phi.

% Operator to map from Q matrices (sym) to polynomials

nphi = length(phi); % Number of monomials
nQ = nphi*(nphi+1)/2; % Number of unique elements in Q
[pr,q] = pr.newFree(nQ);
Q = mss_v2s(q); 

% Update decision variables after 
decvar = pr.variables;

tmp = phi'*Q*phi; % Generic polynomial.

% sosCnst = p - tmp;

% -------------------------------------------------------------------------
% Representation of the original polynomial

Ap = diff(p,theta);
bp = subs(p,theta,0*theta);

% Now: p = Ap*theta + bp;

[vp,pp,cp] = decomp([bp Ap].');

% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Representation of the Gram polynomial

% Old version: could not handle non-control-affine systems
% Aq = diff(tmp,q);
% bq = subs(tmp,q,0*q);
% 
% % Now: p = Ap*theta + bp;
% [vq,pq,cq] = decomp([bq Aq].');

% NEW version: 12/10/2016
% Can handle...
tmp = tmp + sum(vp);

Aq = diff(tmp,q);
bq = subs(tmp,q,0*q);

% Now: p = Ap*theta + bp;
[vq,pq,cq] = decomp([bq Aq].');

% Search for artificial monomials
% If a monomial only appears in the first row of c, then we know that it
% has been artificially introduced.
artm = find(cq(1,:).*(sum(abs(cq),1) == 1));

% Now remove those monomials from pq:
pq(artm,:) = [];

% -------------------------------------------------------------------------

% For every monomial in p, check to see if it's covered by the gram poly.

rmv = [];

for i = 1:size(pp,1)
    
    md = repmat(pp(i,:),size(pq,1),1);
    chk = sum(sum(md == pq,2) == ones(size(pq,1),1)*size(pq,2));
    
    if ~chk
       
% Here's what we actually care about:
%       We want to check to see if there is a unique element in the constraint,
        nzrs = find(cp(2:end,i)' ~= 0);
        
%       Check to see if it's a straight-up simple set-to-zero constraint,        
        if (length(nzrs) == 1) && (cp(1,i) == 0)
            rmv = [rmv,nzrs];
        end
                 
    end

end

%% Remove monomials from model

spot_monos_e_mod = spot_monos_e; 
spot_monos_f_mod = spot_monos_f; 
spot_monos_g_mod = spot_monos_g; 

% We only wish to proceed if we have to...
if ~isempty(rmv)

e_updated = 0;
f_updated = 0;
g_updated = 0;
    
% We must make the changes in descending order, otherwise the indexes will no longer 
% be valid after we remove entries from the polynomials.
% Alternatively, we could probably just insert zeros?
% Update: looks like we might have to insert zeros...
rmv = sort(rmv,'descend');

e_pol_mod = e_pol;
f_pol_mod = f_pol;
g_pol_mod = g_pol;

%     for i = 1:length(rmv)
    for i = rmv

%         i_rmv = rmv(i);

        if i <= ncoef(1)
% This monomial effects e:
            e_pol_mod(:,i) = zeros(nx,1);  
            e_updated = 1; 
        elseif i <= sum(ncoef(1:2))
% This .. f:    
            f_pol_mod(:,i-ncoef(1)) = zeros(nx,1);
            f_updated = 1; %
        else
% This...g: this will probablt never get used...
            g_pol_mod(:,i-sum(ncoef(1:2))) = zeros(ny,1);
            g_updated = 1;
        end
        
    end
    
% Once we've removed all the monomials, we can generate new anonymous functions:

%% Modifying e

if e_updated
    
    for i = 1:nx    
        spot_monos_e{i} = e_pol_mod(i,:);
    end    
    
%     [epow,eco,nem] = mono2mat(spot_monos_e,x);
%     e_monos_mod = @(x) gen_e_monos(x,epow,eco,nx,nem);
    
end

%% Modifying f

if f_updated
    
    for i = 1:nx    
        spot_monos_f_mod{i} = f_pol_mod(i,:);
    end    
    
%     xu2v = map_xu2v(x,u);
%     [fpow,fco,nfm] = f_mono2mat(spot_monos_f,xu2v,x,u);
%     f_monos_mod = @(x,u) gen_e_monos(xu2v*[x;u],fpow,fco,nx,nfm);
    
end

%% Modifying g

if g_updated
    
    for i = 1:ny   
        spot_monos_g_mod{i} = g_pol_mod(i,:);
    end  
    
%     [gpow,gco,ngm] = g_mono2mat(spot_monos_g,x,ny);
%     g_monos_mod = @(x) gen_e_monos(x,gpow,gco,ny,ngm);
    
end

    
end

%% Return results

res.e_monos = spot_monos_e_mod;
res.f_monos = spot_monos_f_mod;
res.g_monos = spot_monos_g_mod;


end

