function res = lr_sos2lmiPlusLin_2(pr,e,f,g,x,theta,varargin)
% Function to encode stability condition (sos constraint) as a LMI.

% Inputs
% - pr: spotsosprog
% - e : mss_poly 
% - f : mss_poly
% - g : mss_poly

% - theta: (spot variables) vector of model parameters


% Outputs


%% Preprocessing

nthp = length(theta);

[nx,~] = size(x);
[ny,~] = size(g);

% if ~isempty(varargin)
%     options = varagin{1};
% else

if ~isempty(varargin)
    options = varargin{1};
else
    options = [];
end

if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 1;   
end

%% Stability constraint

% First step is to generate the polynomial that encodes the stability
% constraint

% We'll introduce P inside this function:
[pr,P] = pr.newSym(nx);

% Augment the parameter vector with P.
theta = [theta;mss_s2v(P)];
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

Aq = diff(tmp,q);
bq = subs(tmp,q,0*q);

% Now: p = Ap*theta + bp;

[vq,pq,cq] = decomp([bq Aq].');

% -------------------------------------------------------------------------

% For every monomial in p, check to see if it's covered by the gram poly.

uncvd = [];

Ae = [];
be = [];

mono_deduct = 0;

for i = 1:size(pp,1)
    
    md = repmat(pp(i,:),size(pq,1),1);
    chk = sum(sum(md == pq,2) == ones(size(pq,1),1)*size(pq,2));
    
    if ~chk
% This monomial is not covered by the Gram polynomial, and requires an 
% equality constraint:
        Ae = [Ae; cp(2:end,i)'];
        be = [be; -cp(1,i)];
        
% Build the monomials that we must deduct from p
        mono_deduct = mono_deduct + (cp(:,i)'*[1;theta])*prod(vp.^(pp(i,:)'));
        
        uncvd(end+1) = i;
    end

end

% Remove the redundant monomials:
p_mod = p - mono_deduct;

% -------------------------------------------------------------------------
% Representation of p_mod

Apm = diff(p_mod,theta);
bpm = subs(p_mod,theta,0*theta);

% Now: p = Ap*theta + bp;

[vpm,ppm,cpm] = decomp([bpm Apm].');

%% New code!

% So, there were some problems when Ae was empty...

if isempty(Ae)
    if verbose >= 1; fprintf('\n\t No. of equality constraints: %d\n',0); end;
    Ae = zeros(1,nth);
    be = 0;
else
    if verbose >= 1; fprintf('\t No. of equality constraints: %d\n',size(Ae,1)); end;
end
% -------------------------------------------------------------------------

%% Nullspace of the equality constraint

% NA = null(full(Ae));


%% L: Q -> p

% I'm pretty sure that we get the L oeprator for free...

L = cq(2:end,:)'; % We remove the first row because there are no constant terms in the Gram poly


%% Determine the nullspace of L:
% NL = null(full(L)); 
NL = spnull(sparse(L));
[~,dimNL] = size(NL);

if verbose >= 1; fprintf('\t Dimension of N(L): %d\n',dimNL); end;

% figure
% spy(spnull(sparse(L)))
% 
% fk=1;

%% Build A operator: theta -> Q

% REMEMBER: Q = A*theta + b (i.e. don't forget the constant).

A = zeros(nQ,nth);

% Build the matrix by considering each element of theta separately:

for i = 1:nth
       
% We can extract the monomial INDEXES from cpm
    monindxvec = cpm(i+1,:); % Note: i+1 to skip first row (constant).

    monindx = find(monindxvec ~= 0);

% Now, for each monomial we want to find the appropriate row of the L operator:
    for k = monindx
       
        m = ppm(k,:); % Monomial
        
        r = sum(repmat(m,size(pq,1),1) == pq,2) == size(pq,2); % Corresponding row of L
        
        cols = find(L(r,:) ~= 0); % Find an element of L in this row (corresponds to the element of Q)
        
        
%       How often do we have a choice?
%         if length(cols) > 1
%             look = 1;
%         end
        
        c = cols(1); % For simplicity, take the first one

        A(c,i) = monindxvec(k)/L(r,c);

    end
    
end

%% Now we have to handle the constant term: b.

b = zeros(nQ,1);

% We can extract the monomial INDEXES from cpm
monindxvec = cpm(1,:); % Note: first row (constant).

monindx = find(monindxvec ~= 0);

% Now, for each monomial we want to find the appropriate row of the L operator:
for k = monindx
    
    m = ppm(k,:); % Monomial

    r = sum(repmat(m,size(pq,1),1) == pq,2) == size(pq,2); % Corresponding row of L

    cols = find(L(r,:) ~= 0); % Find an element of L in this row (corresponds to the element of Q)

%       How often do we have a choice?
%         if length(cols) > 1
%             look = 1;
%         end    
    
    c = cols(1); % For simplicity, take the first one
    
    b(c) = monindxvec(k)/L(r,c);

end

%% This is unfortunate...

% Here is a matrix that converts a vec representing s2v(Z) to a vec
% representing Z(:)

s2v2s = zeros(nphi^2,nQ);

for i = 1:nphi^2
    
    tmp = zeros(nphi);
    tmp(i) = 1;
    tmp = tmp + tmp';
    
    tmp = double(mss_s2v(tmp)~=0);
    s2v2s(i,:) = tmp';
end

% Was this sparse originally? Probably not...
% res.s2v2s = sparse(s2v2s);
res.s2v2s = (s2v2s);


%% Output the results

% Perhaps we'll make the output a little excessive at first (can always
% cull unnecessary stuff later).

res.Almi = ([A, NL]);
res.blmi = (b);

res.Aeq = ([Ae zeros(size(Ae,1),dimNL)]);
res.beq = (be);

res.dimNL = dimNL;

% Parametrization of space satisfying full equality constraints:
Neq = (null(full(res.Aeq)));

% Note: Neq is probably a better name, but for backward compatibility we
% will use N (basically, it means we don't have to rewrite the `bound' func).
res.N = Neq;
res.N1 = Neq(1:end-dimNL,:); % Mapping to model parameters
% res.N1 = sparse(Neq(1:end-dimNL,:)); % Note: it's not clear (empirically)
% that sparse arithmetic is faster than dense (10/4/18)
res.N2 = Neq(end-dimNL+1:end); % Mapping to auxiliary params

res.lmitol = 1e-8; % Tolerance for LMI constraint
% res.lmitol = 1e-1;

res.nphi = nphi; % Number of monomials in basis

% Important: there is no reason to compute this at runtime
res.Mbar = sparse(res.s2v2s*res.Almi*res.N);

% res.AQ = A;
% res.N = N;
% res.A = [A, -N];
% res.b = b;
% res.dimN = dimN;
% res.trs = s;
% res.trlim = trlim;

%% Alternative output (sparse)

% res.Almi = sparse([A, NL]);
% res.blmi = sparse(b);
% 
% res.Aeq = sparse([Ae zeros(size(Ae,1),dimNL)]);
% res.beq = sparse(be);
% 
% res.dimNL = dimNL;
% 
% % Parametrization of space satisfying full equality constraints:
% Neq = sparse(null(full(res.Aeq)));
% 
% % Note: Neq is probably a better name, but for backward compatibility we
% % will use N (basically, it means we don't have to rewrite the `bound' func).
% res.N = Neq;
% res.N1 = Neq(1:end-dimNL,:); % Mapping to model parameters
% res.N2 = Neq(end-dimNL+1:end); % Mapping to auxiliary params
% 
% res.lmitol = 1e-8; % Tolerance for LMI constraint
% 
% res.nphi = nphi; % Number of monomials in basis
% 
% % Very important: there is no reason to compute this at runtime
% res.Mbar = res.s2v2s*res.Almi*res.N;

%% Trace barrier

% THERE IS A LOT ABOUT THIS THAT WON'T GENERALIZE, but for the time being, 
% we just want to see if this will help.

% First we need to somehow extract the right model parameters...
% hc = [1,38,78,122]; % Hard code for now
% 
% Atr = zeros(1,nthp);
% Atr(hc) = 1;
% trlim = nx*1e5;
% 
% res.Atr = sparse(Atr);
% res.trlim = trlim;

res.Atr = 0;
res.trlim = 0;


% What about...if we put a constraint on the trace of Z?
% Find all of the elements that contribute to the trace of Z.
QI = logical(mss_s2v(eye(nphi)));

% Now:
A_trtmp = res.Almi(QI,:);
b_trtmp = res.blmi(QI,:);

res.Atr = sparse(sum(A_trtmp));
res.trlim = nphi*1e3-sum(b_trtmp);

if verbose >= 1; fprintf('\t Maximum trace: %.5e.\n',nphi*1e3); end;

%% An alternative init strat

% QI = (mss_s2v(eye(nphi)));
% 
% Ar = A(QI,:);
% br = b(QI,:);
% 
% % Ai = [res.Almi;res.Aeq];
% % bi = [QI-b;be];
% 
% % Ai = [A;Ae];
% % bi = [QI-b;be];
% 
% initparam = Ar\br;
% 
% res.initparam = initparam;



% end

