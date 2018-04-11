function res = lr_genblks_4(e_pol,f_pol,g_pol,x,u,xd,ud,yd)

% Some speed up by eliminating for loops down the track...

% 03/11/2016 - includes preprocessing for exact computation of the Hessian
% of LR. More precisely, we compute the gradient of Delta, which is then
% used in the Hessian calculation.


% The gradient is given by dM*Delta + M*dDelta = dm.

%% Basic preprocessing

[nx,T] = size(xd);
% nu = size(ud,1);
ny = size(g_pol,1);

ncoef(1) = size(e_pol,2);
ncoef(2) = size(f_pol,2);
ncoef(3) = size(g_pol,2);

%% Precompute stuff

% Handle: e
% ---------

% Build matrices for equation error
% eps = f(x,u) - e(x)
Fbar_e_full = [];
eps_e_full = [];

% Fbar_stk = [];
Fbar_e_stk = [];

for i = 1:ncoef(1)
   
% [e_i(x_2) , e_i(x_3)...
    e_is = -dmsubs(e_pol(:,i),x,xd(:,2:end)); % Note the minus sign
    
% Now we can simply stack:
    eps_e{i} = [zeros(nx,1);e_is(:)];
    
% A little messy, but whatever...
% First, compute the jacobian
    E_i = diff(e_pol(:,i),x);
    E_is = dmsubs(E_i(:),x,xd);
    
% Now we have to loop:
%     tmp = zeros(T*nx,T*nx);
    tmp = sparse(T*nx,T*nx);
    
    for t = 1:T
        tmp((t-1)*nx+1:t*nx,(t-1)*nx+1:t*nx) = reshape(E_is(:,t),nx,nx);
    end
    
    Fbar_e{i} = sparse(tmp);
%     Fbar_e{i} = tmp; % For experimental purposes
    
% Experimental: 20/10/2016    
    Fbar_e_full = [Fbar_e_full; Fbar_e{i}];
    eps_e_full = [eps_e_full [zeros(nx,1);e_is(:)]]; 
    
%     Fbar_stk = [Fbar_stk,Fbar_e{i}(:)];
    Fbar_e_stk = [Fbar_e_stk, Fbar_e{i}(:)];
    
% Hessian stuff:
% -------------------------------------------------------------------------
% We first precompute quantities for dm = [dm/dth1, dm/dth2, ...]   
% Note that eps_e_full is exactly what we need...for the first part of dm.
    
end

res.Fbar_e_full = sparse(Fbar_e_full);
res.eps_e_full = eps_e_full;
res.d1 = T*nx;
% res.spi = sparse(eye(ncoef(1)));

% Handle: f
% ---------

% Build matrices for equation error
% eps = f(x,u) - e(x)

Fbar_f_full = [];
eps_f_full = [];

Fbar_f_stk = [];

for i = 1:ncoef(2)
   
% [f_i(x_1,u_1) , f_i(x_2,u_2)...
    f_is = dmsubs(f_pol(:,i),[x;u],[xd(:,1:end-1);ud(:,1:end-1)]);
    
% Now we can simply stack:
    eps_f{i} = [zeros(nx,1);f_is(:)];
    
% A little messy, but whatever...
% First, compute the jacobian
    F_i = diff(f_pol(:,i),x);
    F_is = dmsubs(F_i(:),[x;u],[xd(:,1:end-1);ud(:,1:end-1)]);
    
% Now we have to loop:
    tmp = zeros(T*nx,T*nx);
    
    for t = 1:T-1
        tmp(nx+(t-1)*nx+1:nx+t*nx,(t-1)*nx+1:t*nx) = -reshape(F_is(:,t),nx,nx);
    end
    
    Fbar_f{i} = sparse(tmp);  
%     Fbar_f{i} = tmp; 
    
    Fbar_f_full = [Fbar_f_full; Fbar_f{i}];
    eps_f_full = [eps_f_full [zeros(nx,1);f_is(:)]]; 
    
%     Fbar_stk = [Fbar_stk, Fbar_f{i}(:)];
    Fbar_f_stk = [Fbar_f_stk, Fbar_f{i}(:)];
    

% Hessian stuff:
% -------------------------------------------------------------------------
% We first precompute quantities for dm = [dm/dth1, dm/dth2, ...]   
% Note that eps_f_full is exactly what we need...for the 2nd part of dm.    
    
end

res.Fbar_f_full = Fbar_f_full;
res.eps_f_full = eps_f_full;
% res.Fbar_stk = Fbar_stk;
res.Fbar_e_stk = Fbar_e_stk;
res.Fbar_f_stk = Fbar_f_stk;

% Handle: g
% ---------

% Build matrices for equation error
% eta = g(x,u) - y

Gbar_g_full = [];
eta_g_full = [];

for i = 1:ncoef(3)
   
% [g_i(x_1,u_1) , g_i(x_2,u_2)...
    g_is = dmsubs(g_pol(:,i),[x;u],[xd;ud]);
    
% Now we can simply stack:
    eta_g{i} = g_is(:);
    
% A little messy, but whatever...
% First, compute the jacobian
    G_i = diff(g_pol(:,i),x);
    G_is = dmsubs(G_i(:),[x;u],[xd;ud]);
    
% Now we have to loop:
    tmp = zeros(T*ny,T*nx);
    
    for t = 1:T
        tmp((t-1)*ny+1:t*ny,(t-1)*nx+1:t*nx) = reshape(G_is(:,t),ny,nx);
    end
    
    Gbar_g{i} = sparse(tmp); 
%     Gbar_g{i} = tmp; 
    
    Gbar_g_full = [Gbar_g_full; Gbar_g{i}];
    eta_g_full = [eta_g_full, g_is(:)]; 
    
    
end

res.Gbar_g_full = Gbar_g_full;
res.eta_g_full = eta_g_full;
res.d2 = T*ny;

%% Completely new section for Hessian stuff

% profile on

% % 24/11/2016 - ALL THIS WAS WORKING JUST FINE (DO NOT MODIFY)
% 
% % The gradient is given by dM*Delta + M*dDelta = dm.
% % We first precompute quantities for dm = [dm/dth1, dm/dth2, ...]
% 
% % Terms related to e and f are handled above; below we compute the matrices
% % related to g.
% 
% m_g = zeros(T*nx,ncoef(3)^2); % 
% m_gc = zeros(T*nx,ncoef(3)); % Constant term
% 
% for t = 1:T % Block rows
%     
%     for k = 1:ncoef(3) % Block columns
%        
% % There are probably more efficient ways to compute this, but because we 
% % are doing this offline, we'll just use loops to make sure that we get it
% % right.
%         tmp_1 = [];
%         tmp_2 = [];
%         
%         g_k = subs(g_pol(:,k),[x;u],[xd(:,t);ud(:,t)]); % kth monomial evaluated at t
%         G_k = subs(diff(g_pol(:,k),x),[x;u],[xd(:,t);ud(:,t)]); % Jacobian of kth monomial 
%         
%         for j = 1:ncoef(3)           
%             tmp_1 = [tmp_1, subs(diff(g_pol(:,j),x),[x;u],[xd(:,t);ud(:,t)])'*g_k];
%             tmp_2 = [tmp_2, G_k'*subs(g_pol(:,j),[x;u],[xd(:,t);ud(:,t)])];            
%         end
%         
%         m_g((t-1)*nx+1:t*nx,(k-1)*ncoef(3)+1:k*ncoef(3)) = tmp_1 + tmp_2;
%         
% % Constant term:
%         m_gc((t-1)*nx+1:t*nx,k) = -G_k'*yd(:,t);
%         
%     end
%         
% end
% 
% % Now we can move onto the more challenging dM term. 
% % We need T*nx separate matrices...one for each row of M. Fortunately, they
% % should be pretty sparse.
% 
% res.m_g = m_g; % Should this be sparse?
% res.m_gc = m_gc;
% 
% profile off
% profile viewer
% dung = 1;


%% Try to compute the above faster

% % profile on
% 
% Many of these quantities have been computed above, but for now, let's
% just recompute...

% g_ks
% G_ks

for k = 1:ncoef(3)
   
% [g_i(x_1,u_1) , g_i(x_2,u_2)...
    g_ks{k} = dmsubs(g_pol(:,k),[x;u],[xd;ud]);

    G_k = diff(g_pol(:,k),x);
    G_ks{k} = dmsubs(G_k(:),[x;u],[xd;ud]);    

end


% The gradient is given by dM*Delta + M*dDelta = dm.
% We first precompute quantities for dm = [dm/dth1, dm/dth2, ...]

% Terms related to e and f are handled above; below we compute the matrices
% related to g.

m_g = zeros(T*nx,ncoef(3)^2); % 
m_gc = zeros(T*nx,ncoef(3)); % Constant term

for t = 1:T % Block rows
    
    for k = 1:ncoef(3) % Block columns
       
% There are probably more efficient ways to compute this, but because we 
% are doing this offline, we'll just use loops to make sure that we get it
% right.
        tmp_1 = zeros(nx,ncoef(3));
        
        g_k = g_ks{k}(:,t); % kth monomial evaluated at t
        G_k = reshape(G_ks{k}(:,t),ny,nx); % Jacobian of kth monomial 
        
        for j = 1:ncoef(3)           
            tmp_1(:,j) = reshape(G_ks{j}(:,t),ny,nx)'*g_k + G_k'*g_ks{j}(:,t);            
        end
        
        m_g((t-1)*nx+1:t*nx,(k-1)*ncoef(3)+1:k*ncoef(3)) = tmp_1;
        
% Constant term:
        m_gc((t-1)*nx+1:t*nx,k) = -G_k'*yd(:,t);
        
    end
        
end

% Now we can move onto the more challenging dM term. 
% We need T*nx separate matrices...one for each row of M. Fortunately, they
% should be pretty sparse.

res.m_g = m_g; % Should this be sparse?
res.m_gc = m_gc;
% 
% % profile off
% % profile viewer
% % dung = 1;


%% Output the results
res.eps_e = eps_e;
res.eps_f = eps_f;
res.eta_g = eta_g;
res.Fbar_e = Fbar_e;
res.Fbar_f = Fbar_f;
res.Gbar_g = Gbar_g;

res.lrnorm = T*nx;
res.yd = yd; % Would be nice to get rid of this...

end































