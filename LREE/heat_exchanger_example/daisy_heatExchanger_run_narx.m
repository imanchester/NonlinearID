%% polynomial regressors

optionalInputRegressors = 1;

r1 = {};

if optionalInputRegressors

    r2 = {...
    'y1(t-1)^2',...
    'y1(t-2)^2',...
    'y1(t-3)^2',...
    'y1(t-1)*y1(t-2)',...
    'y1(t-1)*y1(t-3)',...
    'y1(t-2)*y1(t-3)',...
    'y1(t-1)*u1(t-1)',...
    'y1(t-2)*u1(t-1)',...
    'y1(t-3)*u1(t-1)',...
    'u1(t-1)^2',...
    'y1(t-1)*u1(t)',...
    'y1(t-2)*u1(t)',...
    'y1(t-3)*u1(t)',...
    'u1(t)^2',...
    };




    r3 = {...
    'y1(t-1)^2',...
    'y1(t-2)^2',...
    'y1(t-3)^2',...
    'y1(t-1)*y1(t-2)',...
    'y1(t-1)*y1(t-3)',...
    'y1(t-2)*y1(t-3)',...
    'y1(t-1)^3',...
    'y1(t-2)^3',...
    'y1(t-3)^3',...
    'y1(t-1)*y1(t-2)*y1(t-3)',...
    'y1(t-1)^2*y1(t-2)',...
    'y1(t-1)^2*y1(t-3)',...
    'y1(t-2)^2*y1(t-1)',...
    'y1(t-2)^2*y1(t-3)',...
    'y1(t-3)^2*y1(t-1)',...
    'y1(t-3)^2*y1(t-2)',...
    'y1(t-1)*u1(t-1)',...
    'y1(t-2)*u1(t-1)',...
    'y1(t-3)*u1(t-1)',...
    'y1(t-1)^2*u1(t-1)',...
    'y1(t-2)^2*u1(t-1)',...
    'y1(t-3)^2*u1(t-1)',...
    'y1(t-1)*u1(t-1)^2',...
    'y1(t-2)*u1(t-1)^2',...
    'y1(t-3)*u1(t-1)^2',...
    'y1(t-1)*y1(t-2)*u1(t-1)',...
    'y1(t-1)*y1(t-3)*u1(t-1)',...
    'y1(t-2)*y1(t-3)*u1(t-1)',...
    'u1(t-1)^2',...
    'u1(t-1)^3',...
    'y1(t-1)*u1(t)',...
    'y1(t-2)*u1(t)',...
    'y1(t-3)*u1(t)',...
    'y1(t-1)^2*u1(t)',...
    'y1(t-2)^2*u1(t)',...
    'y1(t-3)^2*u1(t)',...
    'y1(t-1)*u1(t)^2',...
    'y1(t-2)*u1(t)^2',...
    'y1(t-3)*u1(t)^2',...
    'y1(t-1)*y1(t-2)*u1(t)',...
    'y1(t-1)*y1(t-3)*u1(t)',...
    'y1(t-2)*y1(t-3)*u1(t)',...
    'u1(t)^2',...
    'u1(t)^3',...
    };

else


    r2 = {...
    'y1(t-1)^2',...
    'y1(t-2)^2',...
    'y1(t-3)^2',...
    'y1(t-1)*y1(t-2)',...
    'y1(t-1)*y1(t-3)',...
    'y1(t-2)*y1(t-3)',...
    'y1(t-1)*u1(t-1)',...
    'y1(t-2)*u1(t-1)',...
    'y1(t-3)*u1(t-1)',...
    'u1(t-1)^2',...
    'y1(t-1)*u1(t-2)',...
    'y1(t-2)*u1(t-2)',...
    'y1(t-3)*u1(t-2)',...
    'u1(t-2)^2',...
    };




    r3 = {...
    'y1(t-1)^2',...
    'y1(t-2)^2',...
    'y1(t-3)^2',...
    'y1(t-1)*y1(t-2)',...
    'y1(t-1)*y1(t-3)',...
    'y1(t-2)*y1(t-3)',...
    'y1(t-1)^3',...
    'y1(t-2)^3',...
    'y1(t-3)^3',...
    'y1(t-1)*y1(t-2)*y1(t-3)',...
    'y1(t-1)^2*y1(t-2)',...
    'y1(t-1)^2*y1(t-3)',...
    'y1(t-2)^2*y1(t-1)',...
    'y1(t-2)^2*y1(t-3)',...
    'y1(t-3)^2*y1(t-1)',...
    'y1(t-3)^2*y1(t-2)',...
    'y1(t-1)*u1(t-1)',...
    'y1(t-2)*u1(t-1)',...
    'y1(t-3)*u1(t-1)',...
    'y1(t-1)^2*u1(t-1)',...
    'y1(t-2)^2*u1(t-1)',...
    'y1(t-3)^2*u1(t-1)',...
    'y1(t-1)*u1(t-1)^2',...
    'y1(t-2)*u1(t-1)^2',...
    'y1(t-3)*u1(t-1)^2',...
    'y1(t-1)*y1(t-2)*u1(t-1)',...
    'y1(t-1)*y1(t-3)*u1(t-1)',...
    'y1(t-2)*y1(t-3)*u1(t-1)',...
    'u1(t-1)^2',...
    'u1(t-1)^3',...
    'y1(t-1)*u1(t-2)',...
    'y1(t-2)*u1(t-2)',...
    'y1(t-3)*u1(t-2)',...
    'y1(t-1)^2*u1(t-2)',...
    'y1(t-2)^2*u1(t-2)',...
    'y1(t-3)^2*u1(t-2)',...
    'y1(t-1)*u1(t-2)^2',...
    'y1(t-2)*u1(t-2)^2',...
    'y1(t-3)*u1(t-2)^2',...
    'y1(t-1)*y1(t-2)*u1(t-2)',...
    'y1(t-1)*y1(t-3)*u1(t-2)',...
    'y1(t-2)*y1(t-3)*u1(t-2)',...
    'u1(t-2)^2',...
    'u1(t-2)^3',...
    };

end


%% Polynomial model

if run_poly_p1

opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

% Regressors
    r = r1;

t1 = tic;
nlarx_poly = nlarx(narx_dataset_fit, narx_orders,...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly);
close(101)

t_narx_poly = t2;

ys_nlarx_poly = get(nlarx_poly_sim,'outputData')';

se_narx_poly = nse(yd(:,:),ys_nlarx_poly);  
fprintf('Simulation error: %.5e\n',se_narx_poly) 

clear res_tmp
res_tmp.nlarx = nlarx_poly;
res_tmp.se = se_narx_poly;
res_tmp.ys = ys_nlarx_poly;
res_tmp.time = t_narx_poly;

res_narx_poly_p{1} = res_tmp;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly),ys_nlarx_poly,'r--')    
title('NARX. reg = polynomial (linear), focus = prediction')

end


%% Polynomial model

if run_poly_p2

opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

% Regressors
    r = r2;

t1 = tic;
nlarx_poly = nlarx(narx_dataset_fit, narx_orders,...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly);
close(101)

t_narx_poly = t2;

ys_nlarx_poly = get(nlarx_poly_sim,'outputData')';

se_narx_poly = nse(yd(:,:),ys_nlarx_poly);  
fprintf('Simulation error: %.5e\n',se_narx_poly) 

clear res_tmp
res_tmp.nlarx = nlarx_poly;
res_tmp.se = se_narx_poly;
res_tmp.ys = ys_nlarx_poly;
res_tmp.time = t_narx_poly;

res_narx_poly_p{2} = res_tmp;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly),ys_nlarx_poly,'r--')    
title('NARX. reg = polynomial (quadratic), focus = prediction')

end



%% Polynomial model

if run_poly_p3

opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

% Regressors
    r = r3;

t1 = tic;
nlarx_poly = nlarx(narx_dataset_fit, narx_orders,...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly);
close(101)

t_narx_poly = t2;

ys_nlarx_poly = get(nlarx_poly_sim,'outputData')';

se_narx_poly = nse(yd(:,:),ys_nlarx_poly);  
fprintf('Simulation error: %.5e\n',se_narx_poly) 

clear res_tmp
res_tmp.nlarx = nlarx_poly;
res_tmp.se = se_narx_poly;
res_tmp.ys = ys_nlarx_poly;
res_tmp.time = t_narx_poly;

res_narx_poly_p{3} = res_tmp;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly),ys_nlarx_poly,'r--')    
title('NARX. reg = polynomial (cubic), focus = prediction')

end



%% Polynomial model

if run_poly_s1

% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

% Regressors
    r = r1;

t1 = tic;
nlarx_poly = nlarx(narx_dataset_fit, narx_orders,...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly);
close(101)

t_narx_poly = t2;

ys_nlarx_poly = get(nlarx_poly_sim,'outputData')';

se_narx_poly = nse(yd(:,:),ys_nlarx_poly);  
fprintf('Simulation error: %.5e\n',se_narx_poly) 

clear res_tmp
res_tmp.nlarx = nlarx_poly;
res_tmp.se = se_narx_poly;
res_tmp.ys = ys_nlarx_poly;
res_tmp.time = t_narx_poly;

res_narx_poly_s{1} = res_tmp;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly),ys_nlarx_poly,'r--')    
title('NARX. reg = polynomial (linear), focus = simulation')

end



%% Polynomial model

if run_poly_s2

% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

% Regressors
    r = r2;

t1 = tic;
nlarx_poly = nlarx(narx_dataset_fit, narx_orders,...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly);
close(101)

t_narx_poly = t2;

ys_nlarx_poly = get(nlarx_poly_sim,'outputData')';

se_narx_poly = nse(yd(:,:),ys_nlarx_poly);  
fprintf('Simulation error: %.5e\n',se_narx_poly) 

clear res_tmp
res_tmp.nlarx = nlarx_poly;
res_tmp.se = se_narx_poly;
res_tmp.ys = ys_nlarx_poly;
res_tmp.time = t_narx_poly;

res_narx_poly_s{2} = res_tmp;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly),ys_nlarx_poly,'r--')    
title('NARX. reg = polynomial (quadratic), focus = simulation')

end





%% Polynomial model

if run_poly_s3

% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

% Regressors
    r = r3;

t1 = tic;
nlarx_poly = nlarx(narx_dataset_fit, narx_orders,...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly);
close(101)

t_narx_poly = t2;

ys_nlarx_poly = get(nlarx_poly_sim,'outputData')';

se_narx_poly = nse(yd(:,:),ys_nlarx_poly);  
fprintf('Simulation error: %.5e\n',se_narx_poly) 

clear res_tmp
res_tmp.nlarx = nlarx_poly;
res_tmp.se = se_narx_poly;
res_tmp.ys = ys_nlarx_poly;
res_tmp.time = t_narx_poly;

res_narx_poly_s{3} = res_tmp;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly),ys_nlarx_poly,'r--')    
title('NARX. reg = polynomial (cubic), focus = simulation')

end





%%

if run_sig_p
    
opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

t1 = tic;

% nlarx_sig_ss = nlarx(narx_dataset_fit,[3 2 0], [sigmoidnet], 'nlreg', 'search', opt);
nlarx_sig_p = nlarx(narx_dataset_fit,[3 1 0], [sigmoidnet], opt);
t2 = toc(t1);
figure(101)
nlarx_optml_sim = compare(narx_dataset_fit,nlarx_sig_p);
close(101)

t_narx_sig_p = t2;

ys_narx = get(nlarx_optml_sim,'outputData')';

%     res_narxss = results{expindex}.narxss;

se_narx_sig_p = nse(yd(:,:),ys_narx(1,:));
fprintf('Simulation error (sig p): %.5e\n',se_narx_sig_p)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx),ys_narx,'r--')
title(['sig p Training error: ' num2str(se_narx_sig_p)])

end

%%

if run_wav_p

opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

t1 = tic;

% nlarx_sig_ss = nlarx(narx_dataset_fit,[3 2 0], [sigmoidnet], 'nlreg', 'search', opt);
nlarx_wav_p = nlarx(narx_dataset_fit,[3 1 0], [wavenet], opt);
t2 = toc(t1);
figure(101)
nlarx_wav_sim = compare(narx_dataset_fit,nlarx_wav_p);
close(101)

t_narx_wav_p = t2;

ys_narx_wav = get(nlarx_wav_sim,'outputData')';

se_narx_wav_p = nse(yd(:,:),ys_narx_wav(1,:));
fprintf('Simulation error (wav p): %.5e\n',se_narx_wav_p)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx_wav),ys_narx_wav,'r--')
title(['wav p Training error: ' num2str(se_narx_wav_p)])

end


%%

if run_sig_s
    
% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

t1 = tic;

% nlarx_sig_ss = nlarx(narx_dataset_fit,[3 2 0], [sigmoidnet], 'nlreg', 'search', opt);
nlarx_sig_s = nlarx(narx_dataset_fit,[3 1 0], [sigmoidnet], opt);
t2 = toc(t1);
figure(101)
nlarx_optml_sim = compare(narx_dataset_fit,nlarx_sig_s);
close(101)

t_narx_sig_s = t2;

ys_narx = get(nlarx_optml_sim,'outputData')';

%     res_narxss = results{expindex}.narxss;

se_narx_sig_s = nse(yd(:,:),ys_narx(1,:));
fprintf('Simulation error (sig s): %.5e\n',se_narx_sig_s)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx),ys_narx,'r--')
title(['sig s Training error: ' num2str(se_narx_sig_s)])

end


if run_wav_s

% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

t1 = tic;

% nlarx_sig_ss = nlarx(narx_dataset_fit,[3 2 0], [sigmoidnet], 'nlreg', 'search', opt);
nlarx_wav_s = nlarx(narx_dataset_fit,[3 1 0], [wavenet], opt);
t2 = toc(t1);
figure(101)
nlarx_wav_sim = compare(narx_dataset_fit,nlarx_wav_s);
close(101)

t_narx_wav_s = t2;

ys_narx_wav = get(nlarx_wav_sim,'outputData')';

se_narx_wav_s = nse(yd(:,:),ys_narx_wav(1,:));
fprintf('Simulation error (wav s): %.5e\n',se_narx_wav_s)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx_wav),ys_narx_wav,'r--')
title(['wav s Training error: ' num2str(se_narx_wav_s)])

end

%%

if run_sig_ss
    
% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

t1 = tic;

nlarx_sig_ss = nlarx(narx_dataset_fit,[3 1 0], [sigmoidnet], 'nlreg', 'search', opt);
% nlarx_sig_ss = nlarx(narx_dataset_fit,[3 2 0], [sigmoidnet], opt);
t2 = toc(t1);
figure(101)
nlarx_optml_sim = compare(narx_dataset_fit,nlarx_sig_ss);
close(101)

t_narx_sig_ss = t2;

ys_narx = get(nlarx_optml_sim,'outputData')';

%     res_narxss = results{expindex}.narxss;

se_narx_sig_ss = nse(yd(:,:),ys_narx(1,:));
fprintf('Simulation error (sig ss): %.5e\n',se_narx_sig_ss)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx),ys_narx,'r--')
title(['sig ss Training error: ' num2str(se_narx_sig_ss)])

end


if run_wav_ss

% opt = nlarxOptions('Focus','Prediction');
opt = nlarxOptions('Focus','Simulation');

t1 = tic;

nlarx_wav_ss = nlarx(narx_dataset_fit,[3 1 0], [wavenet], 'nlreg', 'search', opt);
% nlarx_wav_s = nlarx(narx_dataset_fit,[3 2 0], [wavenet], opt);
t2 = toc(t1);
figure(101)
nlarx_wav_sim = compare(narx_dataset_fit,nlarx_wav_ss);
close(101)

t_narx_wav_ss = t2;

ys_narx_wav = get(nlarx_wav_sim,'outputData')';

se_narx_wav_ss = nse(yd(:,:),ys_narx_wav(1,:));
fprintf('Simulation error (wav ss): %.5e\n',se_narx_wav_ss)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx_wav),ys_narx_wav,'r--')
title(['wav ss Training error: ' num2str(se_narx_wav_ss)])

end

%

if run_sig_ps
    
opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

t1 = tic;

nlarx_sig_ps = nlarx(narx_dataset_fit,[3 1 0], [sigmoidnet], 'nlreg', 'search', opt);
% nlarx_sig_ss = nlarx(narx_dataset_fit,[3 2 0], [sigmoidnet], opt);
t2 = toc(t1);
figure(101)
nlarx_optml_sim = compare(narx_dataset_fit,nlarx_sig_ps);
close(101)

t_narx_sig_ps = t2;

ys_narx = get(nlarx_optml_sim,'outputData')';

%     res_narxss = results{expindex}.narxss;

se_narx_sig_ps = nse(yd(:,:),ys_narx(1,:));
fprintf('Simulation error (sig ps): %.5e\n',se_narx_sig_ps)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx),ys_narx,'r--')
title(['sig ps Training error: ' num2str(se_narx_sig_ps)])

end


if run_wav_ps

opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

t1 = tic;

nlarx_wav_ps = nlarx(narx_dataset_fit,[3 1 0], [wavenet], 'nlreg', 'search', opt);
% nlarx_wav_s = nlarx(narx_dataset_fit,[3 2 0], [wavenet], opt);
t2 = toc(t1);
figure(101)
nlarx_wav_sim = compare(narx_dataset_fit,nlarx_wav_ps);
close(101)

t_narx_wav_ps = t2;

ys_narx_wav = get(nlarx_wav_sim,'outputData')';

se_narx_wav_ps = nse(yd(:,:),ys_narx_wav(1,:));
fprintf('Simulation error (wav ps): %.5e\n',se_narx_wav_ps)  

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_narx_wav),ys_narx_wav,'r--')
title(['wav ps Training error: ' num2str(se_narx_wav_ps)])

end



%% poly prediction

if run_poly_p

opt = nlarxOptions('Focus','Prediction');
% opt = nlarxOptions('Focus','Simulation');

% Regressors

r = {...
% 'y1(t-1)^7',...
'y1(t-1)^2',...
'y1(t-2)^2',...
'y1(t-3)^2',...
'y1(t-1)*y1(t-2)',...
'y1(t-1)*y1(t-3)',...
'y1(t-2)*y1(t-3)',...
'y1(t-1)^3',...
'y1(t-2)^3',...
'y1(t-3)^3',...
'y1(t-1)*y1(t-2)*y1(t-3)',...
'y1(t-1)^2*y1(t-2)',...
'y1(t-1)^2*y1(t-3)',...
'y1(t-2)^2*y1(t-1)',...
'y1(t-2)^2*y1(t-3)',...
'y1(t-3)^2*y1(t-1)',...
'y1(t-3)^2*y1(t-2)',...
'y1(t-1)*u1(t-1)',...
'y1(t-2)*u1(t-1)',...
'y1(t-3)*u1(t-1)',...
'y1(t-1)^2*u1(t-1)',...
'y1(t-2)^2*u1(t-1)',...
'y1(t-3)^2*u1(t-1)',...
'y1(t-1)*u1(t-1)^2',...
'y1(t-2)*u1(t-1)^2',...
'y1(t-3)*u1(t-1)^2',...
'y1(t-1)*y1(t-2)*u1(t-1)',...
'y1(t-1)*y1(t-3)*u1(t-1)',...
'y1(t-2)*y1(t-3)*u1(t-1)',...
'u1(t-1)^2',...
'u1(t-1)^3',...
'y1(t-1)*u1(t)',...
'y1(t-2)*u1(t)',...
'y1(t-3)*u1(t)',...
'y1(t-1)^2*u1(t)',...
'y1(t-2)^2*u1(t)',...
'y1(t-3)^2*u1(t)',...
'y1(t-1)*u1(t)^2',...
'y1(t-2)*u1(t)^2',...
'y1(t-3)*u1(t)^2',...
'y1(t-1)*y1(t-2)*u1(t)',...
'y1(t-1)*y1(t-3)*u1(t)',...
'y1(t-2)*y1(t-3)*u1(t)',...
'u1(t)^2',...
'u1(t)^3',...
};



t1 = tic;
nlarx_poly_p = nlarx(narx_dataset_fit,[3 1 0],...
              wavenet,... % This get's ignored...
              'CustomRegressors',...
              r,...
              'nlreg',....
              {[]},...  % because this says that no regs should go through the nonlinearity    
               opt); 
t2 = toc(t1);
figure(101)
nlarx_poly_sim = compare(narx_dataset_fit,nlarx_poly_p);
close(101)

t_narx_poly_p = t2;

ys_nlarx_poly_p = get(nlarx_poly_sim,'outputData')';

se_narx_poly_p = nse(yd(:,:),ys_nlarx_poly_p);  
fprintf('Simulation error: %.5e\n',se_narx_poly_p) 
    
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_nlarx_poly_p),ys_nlarx_poly_p,'r--')    

end
