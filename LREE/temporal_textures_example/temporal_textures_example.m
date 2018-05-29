%% temporal textures example

clear all


%% Extract data

% Data available at: http://vismod.media.mit.edu/pub/szummer/temporal-texture/raw/

% [im, dims] = datread_jack('steam.y.sub.c', 96, 176, 120,'uchar');
[im, dims] = datread_jack('steam.y.sub2', 115, 170, 120,'uchar');

%% View a particular matrix as a grayscale image

imsel = 100;

I = mat2gray(im(:,:,imsel));

figure
imshow(I)

%% Extract all images

Is = zeros(dims);

for i = 1:dims(3)
    Is(:,:,i) = mat2gray(im(:,:,i));
end


%% Process data

[nr,nc,ni] = size(Is); 

ny = nr*nc;

y = zeros(ny,ni);

for i = 1:ni
    tmp = Is(:,:,i);
    y(:,i) = tmp(:);
end

% Remove offset?
y0 = mean(y,2);

%% Subspace identification

% Because we have no inputs, the subspace algorithm is basically just an
% SVD (i.e. we don't need any projection to compensate for the effect of
% inputs)

T = 80;
iss = 6;
nx = 40;

yd = y(:,1:T) - y0;

[Y,yss] = empirical_hankel(yd,iss);

[U,S,V] = svd(Y,'econ');

tmp = S*V';

xss = tmp(1:nx,:);

%% Fit a model with least squares

Als = xss(:,2:end)/xss(:,1:end-1);
Cls = yss/xss;

%% Fit stable model with lacy and bernstein

fprintf('Fitting model with Lacy & Bernstein method...\n')

yalmip('clear')

Q = sdpvar(nx,nx,'full');
P = sdpvar(nx,nx,'symmetric');

tmp = [P - 1e-10*eye(nx), Q; Q', P];

const = [tmp >= 0];

tmp = Als*P - Q;

obj = tmp(:)'*tmp(:);

t1 = tic;
options = sdpsettings('solver','mosek');
options.verbose = 0;
sol = optimize(const,obj,options);
time_lb = toc(t1);

fprintf('%s\n',sol.info)

Alb = double(Q)/double(P);


%% Fit a stable model with lree

fprintf('Fitting model with LREE...\n')


yalmip('clear')

% Form empirical covariance matrix
z = [xss(:,2:end);xss(:,1:end-1)];
Z = z*z';

nz = 2*nx;

E = sdpvar(nx,nx,'full');
F = sdpvar(nx,nx,'full');
P = sdpvar(nx,nx,'symmetric');
R = sdpvar(nz,nz,'symmetric');

H = 1e-4*eye(nx);
% H = 1;

obj = trace(R*Z);

EF = [E,-F];


tmp = [R,       EF'*H ;
       H'*EF,  H'*E + E'*H - eye(nx)];

const = [tmp >= 0];   

const = [const, [E + E' - P - 1e-12*eye(nx), F'; F, P] >= 0];

t1 = tic;
options = sdpsettings('solver','mosek');
options.verbose = 0;
sol = optimize(const,obj,options);
time_lr = toc(t1);

Alr = double(E)\double(F);

ealr = eig(Alr);
sr_lr = max(diag(sqrt(ealr*ealr')));
fprintf('sr (lree): %.5f\n',sr_lr)


%% Constraint generation

fprintf('Fitting model with constraint generation method...\n')


% First set-up the objective:
M = (zeros((size(xss,2)-1)*nx,nx*(nx)));
n = xss(nx+1:end)';
for i = 1:size(xss,2)-1    
    M((i-1)*nx+1:i*nx,:) = [kron(sparse(eye(nx)),xss(:,i)')];
end
P = (M'*M);
p = M'*n;
r = trace(xss(:,2:end)*xss(:,2:end)');

P = kron(eye(nx),xss(:,1:end-1)*xss(:,1:end-1)');
tmp = xss(:,1:end-1)*xss(:,2:end)';
p = tmp(:);
r = trace(xss(:,2:end)*xss(:,2:end)');

% Check the least squares soln (no need to recompute)
Acg = Als;

eals = eig(Als);
sr_cg = max(diag(sqrt(eals*eals')));

cg_iter = 1;
maxiters = 50;

G = [];
h = [];

t1 = tic;

while (sr_cg > 1) && (cg_iter < maxiters)
   
% Generate constraint
    [Ucg,Scg,Vcg] = svd(Acg);
    
%     tmp = Ucg(:,1)*Vcg(:,1)';
    tmp = Vcg(:,1)*Ucg(:,1)';
    G = [G;tmp(:)'];
    h = [h;1];
    
    yalmip('clear')
    
    a = sdpvar(nx*nx,1);
    
    obj = a'*P*a - 2*a'*p;
    
    constraint = [G*a <= h];
%     constraint = [];
    
    options = sdpsettings('solver','mosek');
    options.verbose = 0;
    sol = optimize(constraint,obj,options);  
    
    a = double(a);
    
    Acg = reshape(a,nx,nx)';
    
    eacg = eig(Acg);
    sr_cg = max(diag(sqrt(eacg*eacg')));
    
    fprintf('It: %d, sr: %.5f\n',cg_iter,sr_cg)
    
    cg_iter = cg_iter + 1;
    
end

time_cg = toc(t1);

%% Reconstruction error 
% Basically check the equation error of the identified models
% Seems like a poor metric?

tmp = xss(:,2:end) - Als*xss(:,1:end-1);
Q = tmp*tmp';
Q_ls = Q/(size(xss,2));
ee_ls = trace(Q)
eals = eig(Als);
sr_ls = max(diag(sqrt(eals*eals')))

tmp = xss(:,2:end) - Alr*xss(:,1:end-1);
Q = tmp*tmp';
Q_lr = Q/(size(xss,2));
ee_lr = trace(Q)
ealr = eig(Alr);
sr_lr = max(diag(sqrt(ealr*ealr')))

tmp = xss(:,2:end) - Alb*xss(:,1:end-1);
Q = tmp*tmp';
Q_lb = Q/(size(xss,2));
ee_lb = trace(Q)
ealb = eig(Alb);
sr_lb = max(diag(sqrt(ealb*ealb')))

tmp = xss(:,2:end) - Acg*xss(:,1:end-1);
Q = tmp*tmp';
Q_cg = Q/(size(xss,2));
ee_cg = trace(Q)
eacg = eig(Acg);
sr_cg = max(diag(sqrt(eacg*eacg')))

fprintf('Normalized reconstruction errors:\n')
e_lr = 100*(ee_lr - ee_ls)/ee_ls;
fprintf('\tLREE: %.5e\n',e_lr)

e_lb = 100*(ee_lb - ee_ls)/ee_ls;
fprintf('\tL&B: %.5e\n',e_lb)


e_cg = 100*(ee_cg - ee_ls)/ee_ls;
fprintf('\tCG: %.5e\n',e_cg)


%% Simulate to synthesize images

Ts = 1000;

ys_ls = zeros(ny,Ts);
ys_lr = zeros(ny,Ts);
ys_lb = zeros(ny,Ts);
ys_cg = zeros(ny,Ts);

% What kind of initial conditions are appropriate?
xt_ls = xss(:,1); %?
xt_lr = xss(:,1); %?
xt_lb = xss(:,1); %?
xt_cg = xss(:,1); %?

% Driven by white noise
w = 2*mvnrnd(zeros(nx,1),eye(nx),Ts)';

g_ls = Q_ls^0.5;
g_lr = Q_lr^0.5;
g_lb = Q_lb^0.5;
g_cg = Q_cg^0.5;


for t = 1:Ts
    ys_ls(:,t) = Cls*xt_ls;
    xt_ls = Als*xt_ls + g_ls*w(:,t);
    
    ys_lr(:,t) = Cls*xt_lr;
    xt_lr = Alr*xt_lr + g_lr*w(:,t); 
    
    ys_lb(:,t) = Cls*xt_lb;
    xt_lb = Alb*xt_lb + g_lb*w(:,t); 
    
    ys_cg(:,t) = Cls*xt_cg;
    xt_cg = Acg*xt_cg + g_cg*w(:,t);
end
    

%% Inspect images

t = 400;

It_ls = reshape(ys_ls(:,t) + y0,nr,nc);
It_lr = reshape(ys_lr(:,t) + y0,nr,nc);


figure
subplot(2,1,1)
imshow(It_ls)
title(['t = ' num2str(t) ', least squares'])

subplot(2,1,2)
imshow(It_lr)
title(['LR'])


%% Figure for paper

% % Compare LS and LR because it's cool?
% ts = [10,500]
% n = 2;
% 
% figure
% 
% t = ts(1);
% subplot(2,n,1)
% It_ls = reshape(ys_ls(:,t) + y0,nr,nc);
% imshow(It_ls)
% 
% p = get(gca,'position');
% p(1) = 1*p(1);
% p(2) = 1*p(2);
% p(3) = 1*p(3);
% p(4) = 1*p(4);
% set(gca,'position',p);
% 
% subplot(2,n,n+1)
% It_lr = reshape(ys_lr(:,t) + y0,nr,nc);
% imshow(It_lr)
% 
% p = get(gca,'position');
% p(1) = 1*p(1);
% p(2) = 1*p(2);
% p(3) = 1*p(3);
% p(4) = 2*p(4);
% set(gca,'position',p);
% xlabel(['$t=' num2str(ts(1)) '$'],'interpreter','latex')
% 
% t = ts(2);
% subplot(2,n,2)
% It_ls = reshape(ys_ls(:,t) + y0,nr,nc);
% imshow(It_ls)
% 
% p = get(gca,'position');
% p(1) = 0.85*p(1);
% p(2) = 1*p(2);
% p(3) = 1*p(3);
% p(4) = 1*p(4);
% set(gca,'position',p);
% 
% subplot(2,n,n+2)
% It_lr = reshape(ys_lr(:,t) + y0,nr,nc);
% imshow(It_lr)
% 
% p = get(gca,'position');
% p(1) = 0.85*p(1);
% p(2) = 1*p(2);
% p(3) = 1*p(3);
% p(4) = 2*p(4);
% set(gca,'position',p);
% xlabel(['$t=' num2str(ts(2)) '$'],'interpreter','latex')


%% Figure for paper

ts = [10,100,1000];
n = length(ts);

figure

for i = 1:n

t = ts(i);
subplot(2,n,i)
It_ls = reshape(ys_ls(:,t) + y0,nr,nc);
imshow(It_ls)

p = get(gca,'position');
p(1) = 0.8*p(1);
p(2) = 1*p(2);
p(3) = 1*p(3);
p(4) = 1*p(4);
set(gca,'position',p);

subplot(2,n,n+i)
It_lr = reshape(ys_lr(:,t) + y0,nr,nc);
imshow(It_lr)

p = get(gca,'position');
p(1) = 0.8*p(1);
p(2) = 1*p(2);
p(3) = 1*p(3);
p(4) = 2.6*p(4);
set(gca,'position',p);
xlabel(['$t=' num2str(ts(i)) '$'],'interpreter','latex')

end

%% Animation

if 1

figure(100)

pt = 0;

nsp = 4;

for t = 1:Ts
    
    It_ls = reshape(ys_ls(:,t) + y0,nr,nc);
    It_lr = reshape(ys_lr(:,t) + y0,nr,nc);
    It_lb = reshape(ys_lb(:,t) + y0,nr,nc);
    It_cg = reshape(ys_cg(:,t) + y0,nr,nc);
    
%     figure(100)
%     subplot(1,nsp,1)
%     imshow(It_ls)
% %     title(['t = ' num2str(t) ', least squares'])
%     title(['t = ' num2str(t) ', (1)'])

    figure(100)
    subplot(1,nsp,1)
    imshow(It_cg)
%     title(['t = ' num2str(t) ', least squares'])
    title(['t = ' num2str(t) ', constraint generation'])


    subplot(1,nsp,2)
    imshow(It_lr)
    title(['Lagrangian relaxation'])
%     title(['(2)'])
    
    subplot(1,nsp,3)
    tt = mod(t,ni) + 1;
    imshow(Is(:,:,tt))
    title(['True (time = ' num2str(tt) ')']) 

    
    subplot(1,nsp,4)
    imshow(It_lb)
    title(['Lacy&Bernstein'])   
%     title(['(3)']) 
    

    
    drawnow
    pause(pt)
    
end

end

%% Simulate forward on training data

Ts = 30;

ys_ls = zeros(ny,Ts);
ys_lr = zeros(ny,Ts);
ys_lb = zeros(ny,Ts);
ys_cg = zeros(ny,Ts);


% What kind of initial conditions are appropriate?
xt_ls = xss(:,end); %?
xt_lr = xss(:,end); %?
xt_lb = xss(:,end); %?
xt_cg = xss(:,end); %?

% Driven by white noise? 
w = 0*mvnrnd(zeros(nx,1),eye(nx),Ts)'; % Check for the zero
w_ls = mvnrnd(zeros(nx,1),Q_ls,Ts)';
w_lr = mvnrnd(zeros(nx,1),Q_lr,Ts)';
w_lb = mvnrnd(zeros(nx,1),Q_lb,Ts)';
w_cg = mvnrnd(zeros(nx,1),Q_cg,Ts)';

g_ls = Q_ls^0.5;
g_lr = Q_lr^0.5;
g_lb = Q_lb^0.5;
g_cg = Q_cg^0.5;



for t = 1:Ts
    ys_ls(:,t) = Cls*xt_ls;
    xt_ls = Als*xt_ls + g_ls*w(:,t);
    
    ys_lr(:,t) = Cls*xt_lr;
    xt_lr = Alr*xt_lr + g_lr*w(:,t); 
    
    ys_lb(:,t) = Cls*xt_lb;
    xt_lb = Alb*xt_lb + g_lb*w(:,t); 
    
    ys_cg(:,t) = Cls*xt_cg;
    xt_cg = Acg*xt_cg + g_cg*w(:,t);
end

y_comp = y(:,T:T+Ts-1) - y0;


%% Compute sim errors

fprintf('Simulation errors:\n')

y2 = y_comp(:)'*y_comp(:);

tmp = y_comp - ys_ls;
se_ls = tmp(:)'*tmp(:)/y2;
fprintf('\tLS: %.5e\n',se_ls)

tmp = y_comp - ys_lr;
se_lr = tmp(:)'*tmp(:)/y2;
fprintf('\tLREE: %.5e\n',se_lr)

tmp = y_comp - ys_lb;
se_lb = tmp(:)'*tmp(:)/y2;
fprintf('\tL&B: %.5e\n',se_lb)


tmp = y_comp - ys_cg;
se_cg = tmp(:)'*tmp(:)/y2;
fprintf('\tCG: %.5e\n',se_cg)



%% Compare simulated sequences

figure(101)

pt = 2;

nsp = 4;

for t = 1:Ts
% for t = 1

    
    It_ls = reshape(ys_ls(:,t) + y0,nr,nc);
    It_lr = reshape(ys_lr(:,t) + y0,nr,nc);
    It_lb = reshape(ys_lb(:,t) + y0,nr,nc);
    It_cg = reshape(ys_cg(:,t) + y0,nr,nc);
    It_tr = reshape(y_comp(:,t) + y0,nr,nc);
    
%     figure(100)
%     subplot(1,nsp,1)
%     imshow(It_ls)
% %     title(['t = ' num2str(t) ', least squares'])
%     title(['t = ' num2str(t) ', (1)'])

    figure(101)
    subplot(1,nsp,1)
    imshow(It_cg)
%     title(['t = ' num2str(t) ', least squares'])
    title(['t = ' num2str(t) ', constraint generation'])


    subplot(1,nsp,2)
    imshow(It_lr)
    title(['Lagrangian relaxation'])
%     title(['(2)'])
    
    subplot(1,nsp,3)
%     tt = mod(t,ni) + 1;
%     imshow(Is(:,:,tt))
%     title(['True (time = ' num2str(tt) ')']) 
    imshow(It_tr)    
    title(['True']) 

    
    subplot(1,nsp,4)
    imshow(It_lb)
    title(['Lacy&Bernstein'])   
%     title(['(3)']) 
    

    
    drawnow
    pause(pt)
    
end



















