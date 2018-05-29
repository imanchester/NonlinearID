%% Liquid heat exchanger

close all
clear variables
clc 


%% Methods to run

% State-space models:
run_lrse_u2     = 1; % Lagrangian relaxation of linearized simulation error (LRSE)
run_lree_u2     = 1; % Lagrangian relaxation of equation error (LREE)
run_stblee_u2   = 1; % Implicit equation error (with stability) (iEE-s)
run_ee_u2       = 1; % Implicit equation error (without stability) (iEE)

% Alternative methods (can be ignored):
run_lrse        = 0;
run_lree        = 0;
run_lree_ae     = 0;
run_stblee      = 0;
run_lree_ct1    = 0;
run_rie         = 0;
run_lrse_ae     = 0;

% Select state-space model structures, i.e., degree of polynomials in
% (e,f,g)
% 1. (1,1,1)
% 2. (3,1,1)
% 3. (3,3,1)
narx_structure = 0;
lr_model_sel = [1,2,3];


% NARX models:
narx_orders = [3 2 0]; % model orders

% Select model structures:
% 1. lienar
% 2. quadratic
% 3. cubic
poly_model_sel = [1,2,3]; 

% Prediction
run_poly_p1 = 1;
run_poly_p2 = 1;
run_poly_p3 = 1;

% Simulation
run_poly_s1 = 1;
run_poly_s2 = 1;
run_poly_s3 = 1;


% Alternative models (can be ignored):
run_wav_p = 0;
run_sig_p = 0;

% Prediction *
run_sig_ps = 0;
run_wav_ps = 0;

% Simulation
run_sig_s = 0;
run_wav_s = 0;

% Simulation *
run_sig_ss = 0;
run_wav_ss = 0;

run_poly_1 = 0;
run_poly_p = 0;


%% load problem data

% See, ftp://ftp.esat.kuleuven.be/pub/SISTA/data/process_industry/exchanger.txt
% for description and data.

load('heatExchanger_97_002.mat')

%% plot data

figure
subplot(2,1,1)
plot(1:length(data.u),data.u(1,:))
ylabel('Input_1')
title('Heat exchanger: DaISy 97-002')
subplot(2,1,2)
plot(1:length(data.y),data.y(1,:))
ylabel('Output')
xlabel('Time (sec)')

%% basic preprocessing

% Normalize the data
ops.u = 1; % normalize inputs
ops.y = 1; % normalize outputs

res_data = normalizeData(data,ops);
pdata = res_data.data;
nv = res_data.normvals;

%% Data for fitting

Tfit = 1000; % number of samples for training
Tstart = 550; % ignore initial datapoints

ud = pdata.u(:,Tstart:Tstart+Tfit-1);
yd = pdata.y(:,Tstart:Tstart+Tfit-1);

[nu,~] = size(ud);
[ny,~] = size(yd);
nx = 3;

dt = 1; % Sampling time
narx_dataset_fit = iddata(yd',ud',dt);

% Generate states - lags/narx
u_raw = ud;
y_raw = yd;

isf = diag(max(abs(u_raw),[],2));      % Input
% ssf = diag(max(abs(xn),[],2));      % State
osf = diag(max(abs(y_raw),[],2));      % Output

% Stack outputs for states
y_sc = y_raw;
xd = stackInputs(y_sc,nx);
yd = y_sc;

% figure
% subplot(2,1,1)
% plot(1:Tfit,ud(1,:))
% ylabel('Input_1')
% title('Training data')
% 
% subplot(2,1,2)
% plot(1:Tfit,yd(1,:))
% ylabel('Output')
% xlabel('Time (sec)')


%% Load state-space models

daisy_heatExchanger_models_nx3

%% Create data-objects for state-space identification

nlid_data.x = xd;
nlid_data.u = ud;
nlid_data.y = yd; 

ud2 = [ud;ud(2:end),0];

nlid_data_u2.x = xd;
nlid_data_u2.u = ud2;
nlid_data_u2.y = yd; 

%% Run state-space: TRAINING

daisy_heatExchanger_run_our_methods

%% Run narx methods: TRAINING

daisy_heatExchanger_run_narx


%% Validation data

Tval = 2000; % number of data points for validation
uv = pdata.u(:,Tstart+Tfit+1:Tstart+Tfit+Tval);
yv = pdata.y(:,Tstart+Tfit+1:Tstart+Tfit+Tval);

uv2 = [uv;uv(2:end),0];

narx_dataset_val = iddata(yv',uv',dt);

figure
subplot(2,1,1)
plot(1:Tval,uv(1,:))
ylabel('Input_1')
title('Validation data')

subplot(2,1,2)
plot(1:Tval,yv(1,:))
ylabel('Output')
xlabel('Time (sec)')

%% run state-space methods: VALIDATION

daisy_heatExchanger_val_ours

%% run narx methods: VALIDATION

daisy_heatExchanger_val_narx

%% print results

fprintf('\nTraining error:\n')
fprintf('------------------\n')

if run_lree

fprintf('\nLREE:\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end

if run_lree_u2

fprintf('\nLREE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end

if run_lree_ae

fprintf('\nLREE (AE):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree_ae{lrm_i};
    if ~isempty(res_tmp)
        fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)
    end
end
end

if run_lrse

fprintf('\nLRSE:\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lrse{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end

if run_lrse_u2

fprintf('\nLRSE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lrse_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end

if run_stblee_u2

fprintf('\nSTBLEE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_stblee_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end

if run_ee_u2

fprintf('\nEE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_ee_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end

if run_poly_p2

fprintf('\nPoly (pred):\n')
    
for lrm_i = poly_model_sel

    res_tmp = res_narx_poly_p{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end


if run_poly_s2

fprintf('\nPoly (sim):\n')
    
for lrm_i = poly_model_sel

    res_tmp = res_narx_poly_s{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.se)

end
end



%% print results

fprintf('\nValidation error:\n')
fprintf('------------------\n')

if run_lree

fprintf('\nLREE:\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end

if run_lree_u2

fprintf('\nLREE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end

if run_lree_ae

fprintf('\nLREE (AE):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree_ae{lrm_i};
    if ~isempty(res_tmp)
        fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)
    end

end
end

if run_lrse

fprintf('\nLRSE:\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lrse{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end

if run_lrse_u2

fprintf('\nLRSE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lrse_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end

if run_stblee_u2

fprintf('\nSTBLEE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_stblee_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end


if run_ee_u2

fprintf('\nEE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_ee_u2{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end



if run_poly_p2

fprintf('\nPoly (pred):\n')
    
for lrm_i = poly_model_sel

    res_tmp = res_narx_poly_p{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end


if run_poly_s2

fprintf('\nPoly (sim):\n')
    
for lrm_i = poly_model_sel 

    res_tmp = res_narx_poly_s{lrm_i};
    fprintf('\tM%d, error: %.5e\n',lrm_i,res_tmp.ve)

end
end



%% print results

fprintf('\nComp times:\n')
fprintf('------------------\n')

if run_lree

fprintf('\nLREE:\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)

end
end

if run_lree_u2

fprintf('\nLREE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree_u2{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)

end
end

if run_lree_ae

fprintf('\nLREE (AE):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lree_ae{lrm_i};
    if ~isempty(res_tmp)
        fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)
    end
end
end

if run_lrse

fprintf('\nLRSE:\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lrse{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)

end
end

if run_lrse_u2

fprintf('\nLRSE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_lrse_u2{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)

end
end

if run_stblee_u2

fprintf('\nSTBLEE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_stblee_u2{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)

end
end



if run_ee_u2

fprintf('\nEE (U2):\n')
    
for lrm_i = lr_model_sel 

    res_tmp = res_ee_u2{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.solvertime)

end
end


if run_poly_p2

fprintf('\nPoly (pred):\n')
    
for lrm_i = poly_model_sel 

    res_tmp = res_narx_poly_p{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.time)

end
end


if run_poly_s2

fprintf('\nPoly (sim):\n')
    
for lrm_i = poly_model_sel 

    res_tmp = res_narx_poly_s{lrm_i};
    fprintf('\tM%d, time: %.5e\n',lrm_i,res_tmp.time)

end
end


%% figure

fs = 12;

tdom = [1040 1140];

figure
subplot(2,1,1)
plot(1:Tval,uv,'color','k','linewidth',1)
ylabel('Input','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'XTickLabel',[])
set(gca,'fontsize',fs)
axis([tdom -1 1.1])

p = get(gca,'position');
p(1) = 1*p(1);
p(2) = 1*p(2);
p(3) = 1*p(3);
p(4) = 0.2*p(4);
set(gca,'position',p);


subplot(2,1,2)
plot(1:Tval,yv,'color',[0.5,0.5,0.5,0.3],'linewidth',5)
hold on
plot(1:Tval,res_lree_u2{3}.yv,'color','b','linewidth',1.2)
% plot(1:Tval,res_lrse_u2{3}.yv,'--r','linewidth',1.2)
plot(1:Tval,res_narx_poly_s{2}.yv,'--r','linewidth',1.2)
plot(1:Tval,res_narx_poly_s{3}.yv,'-.','color',[0,0.7,0.3],'linewidth',1.2)
axis([tdom -0.9 0.25])
% legend({'Data','LREE (3,3,1)','LRSE (3,3,1)','narx-s (cubic)'},'interpreter','latex','location','southeast')
legend({'Data','LREE (3,3,1)','narx-s (quadratic)','narx-s (cubic)'},'interpreter','latex','location','southeast')
set(gca,'fontsize',fs)
xlabel('Time','interpreter','latex')
ylabel('Output','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
legend boxoff

p = get(gca,'position');
p(1) = 1*p(1);
p(2) = 2.7*p(2);
p(3) = 1*p(3);
p(4) = 0.8*p(4);
set(gca,'position',p);






%% generate synthetic inputs

% valiadtion data
% uvs = uv(1,500:end);
% 
% figure
% plot(1:length(uvs),uvs)
% ylabel('Input_1')
% title('Validation data')

% check autocorrelation
% figure
% [acf, lags] = xcorr(uvs - mean(uvs), 100, 'coeff');
% stem(lags(101:200), acf(101:200), 'Color', [117 112 179] / 256, 'LineWidth', 2);
% xlabel('lag'); 
% ylabel('ACF');

numtrials = 10;

Tv = 2e3;

ydum = zeros(1,Tv);
ydum(end) = 1;

umags = [2.2];

vesa_p = zeros(length(umags),numtrials);
vesa_s = zeros(length(umags),numtrials);


for j = 1:length(umags)

    umag = umags(j);

    fprintf('\n|U| = %.2f\n',umag)

for i = 1:numtrials

    uva = -1 + 2*rand(1,Tv);
    uva = uva*umag/(max(abs(uva)));

    narx_dataset_tmp = iddata(ydum',uva',dt);

 
% -------------------------------------------------------------------------   
    res_tmp = res_narx_poly_p{2};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_tmp,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(ydum(:,:),yv_narx_poly(1,:));
    fprintf('\tValidation error (p): %.5e\n',ve_narx_poly)  
    
    vesa_p(j,i) = ve_narx_poly;
    
    if i == 1
        figure
        plot(1:Tv,ydum)
        hold on
        plot(1:Tv,yv_narx_poly,'r--')
        ax=gca;
        set(ax,'YLim',[-2,2])
        title(['narx (p) Validation error: ' num2str(ve_narx_poly) ', |u| = ' num2str(umag)])
        pause(0.1)
    end


% -------------------------------------------------------------------------   
    res_tmp = res_narx_poly_s{2};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_tmp,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(ydum(:,:),yv_narx_poly(1,:));
    fprintf('\tValidation error (s): %.5e\n',ve_narx_poly)  
    
    vesa_s(j,i) = ve_narx_poly;
    
    if i == 1
        figure
        plot(1:Tv,ydum)
        hold on
        plot(1:Tv,yv_narx_poly,'r--')
        ax=gca;
        set(ax,'YLim',[-2,2])
        title(['narx (s) Validation error: ' num2str(ve_narx_poly) ', |u| = ' num2str(umag)])
        pause(0.1)    
    end

    
end

end

%%

fprintf('Number of divergent models:\n')

fprintf('\t narx (quadratic) - prediction focus:\n')
sum(isnan(vesa_p),2)

fprintf('\t narx (quadratic) - simulation focus:\n')
sum(isnan(vesa_s),2)





































