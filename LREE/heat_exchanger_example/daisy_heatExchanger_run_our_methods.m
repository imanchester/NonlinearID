%% Run LREE 

if run_lree

    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
    lrm = lrms_lree{lrm_i};
    
    nlid_data.e_monos = lrm.e_monos;
    nlid_data.f_monos = lrm.f_monos;
    nlid_data.g_monos = lrm.g_monos;       
    
    
    options_lree.stability = 'state';
    options_lree.multi = 1;

    res_tmp = lree_nonlinear(nlid_data,options_lree);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_2(ud,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (lree): %.5e\n',se_tmp)
    
    res_lree{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['LREE (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])

end

end



%% Run LRSE 

if run_lrse

    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
    lrm = lrms_lrse{lrm_i};
    
    nlid_data.e_monos = lrm.e_monos;
    nlid_data.f_monos = lrm.f_monos;
    nlid_data.g_monos = lrm.g_monos;       
    
    
    options_specialized = [];
    options_specialized.verbose = 0;
    options_specialized.recordse = 0; % Record sim error

    % Some new options:
    options_specialized.trcnst = 0; % Impose limit on the trace
    options_specialized.earlyterm = 1; % Terminate inner loop early

    res_tmp = lrse_nonlin_specialized(nlid_data,options_specialized);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_2(ud,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (lrse): %.5e\n',se_tmp)
    
    res_lrse{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['LRSE (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])

end

end


%% Run LREE (ae)

if run_lree_ae

    
for lrm_i = lr_model_sel 

    if lrm_i ~= 1
% Update the id_data structure with the appropriate model,
    lrm = lrms_lree_ae{lrm_i};
    
    nlid_data.e_monos = lrm.e_monos;
    nlid_data.f_monos = lrm.f_monos;
    nlid_data.g_monos = lrm.g_monos;       
    
    
% Fit model with LR:
    options_lree_ae = [];
    options_lree_ae.stability = 'state'; % Basic
    options_lree_ae.multi = 1;

    res_tmp = lree_nonlinear_ae(nlid_data,options_lree_ae);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_ae(ud,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (lree ae): %.5e\n',se_tmp)
    
    res_lree_ae{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['LREE AE (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])
    end
end

end



%% Run LREE U2

if run_lree_u2

    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
    lrm = lrms_lree_u2{lrm_i};
    
    nlid_data_u2.e_monos = lrm.e_monos;
    nlid_data_u2.f_monos = lrm.f_monos;
    nlid_data_u2.g_monos = lrm.g_monos;       
    
    
    options_lree.stability = 'state';
    options_lree.multi = 1;

    res_tmp = lree_nonlinear(nlid_data_u2,options_lree);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_2(ud2,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (lree u2): %.5e\n',se_tmp)
    
    res_lree_u2{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['LREE U2 (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])

end

end



%% Run stblee u2

if run_stblee_u2

    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
    lrm = lrms_lree_u2{lrm_i};
    
    nlid_data_u2.e_monos = lrm.e_monos;
    nlid_data_u2.f_monos = lrm.f_monos;
    nlid_data_u2.g_monos = lrm.g_monos;       
    
    
    options_ee = []; 
    options_ee.stability = 1; 
    options_ee.verbose = 0; 
    res_tmp = ee_nonlin_fast(nlid_data_u2,options_ee);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_2(ud2,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (stblee u2): %.5e\n',se_tmp)
    
    res_stblee_u2{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['STBLEE U2 (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])

end

end




%% Run ee u2

if run_ee_u2

    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
    lrm = lrms_lree_u2{lrm_i};
    
    nlid_data_u2.e_monos = lrm.e_monos;
    nlid_data_u2.f_monos = lrm.f_monos;
    nlid_data_u2.g_monos = lrm.g_monos;       
    
    
    options_ee = []; 
    options_ee.stability = 0; 
    options_ee.verbose = 0; 
    if lrm_i == 1
        options_ee.linear = 1; 
    else
        options_ee.linear = 0; 
    end

    res_tmp = ee_nonlin_fast(nlid_data_u2,options_ee);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_2(ud2,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (ee u2): %.5e\n',se_tmp)
    
    res_ee_u2{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['EE U2 (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])

end

end



%% Run LRSE (U2)

if run_lrse_u2

    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
    lrm = lrms_lree_u2{lrm_i};
    
    nlid_data_u2.e_monos = lrm.e_monos;
    nlid_data_u2.f_monos = lrm.f_monos;
    nlid_data_u2.g_monos = lrm.g_monos;        
    
    
% Fit model with LR:
    options_lree_ae = [];
    options_lree_ae.stability = 'none'; 
    options_lree_ae.verbose = 0; % Full
    %     options_lree_ae.solver = 'sedumi'; % Full

    res_tmp = lrse_nonlin_specialized(nlid_data_u2,options_lree_ae);

    [se_tmp,ys_tmp,xs_tmp] = se_nonlinear_2(ud2,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);
    res_tmp.se = se_tmp;
    res_tmp.ys = ys_tmp;
    res_tmp.xs = xs_tmp;

    fprintf('Training error (lrse U2): %.5e\n',se_tmp)
    
    res_lrse_u2{lrm_i} = res_tmp;

    % Check output on training data
    figure
    plot(1:length(y_raw),y_raw)
    hold on
    plot(1:length(ys_tmp),ys_tmp,'r--')
    title(['LRSE U2 (M ' num2str(lrm_i) '), Training error: ' num2str(se_tmp)])

end

end



%%

if run_lrse_ae
    
nlid_data_ae = nlid_data;
nlid_data_ae.e_monos = e_monos_ae;
nlid_data_ae.f_monos = f_monos_ae;
nlid_data_ae.g_monos = g_monos_ae;   

% Fit model with LR:
options_lree_ae = [];
options_lree_ae.stability = 'none'; 
options_lree_ae.verbose = 1; % Full
%     options_lree_ae.solver = 'sedumi'; % Full

res_lrse_ae = lr_nonlin_yalmip_ae(nlid_data_ae,options_lree_ae);

[se_lrse_ae,ys_lrse_ae,xs_lrse_ae] = se_nonlinear_ae(ud,yd,xd(:,1),lrm_ae.e,lrm_ae.f,lrm_ae.g,res_lrse_ae.ec,res_lrse_ae.fc,res_lrse_ae.gc);
res_lree_ae.se = se_lrse_ae;

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_lrse_ae),ys_lrse_ae,'r--')
title(['LRSE (AE) Training error: ' num2str(se_lrse_ae)])

end





%% Run LREE with a continuous time model

if run_lree_ct1

nlid_data_ct1 = nlid_data;
nlid_data_ct1.f_monos = res_prepro_ct1.f_monos;
nlid_data_ct1.g_monos = res_prepro_ct1.g_monos;  
nlid_data_ct1.dt = dt;

% Fit model with LR:
options_lree_ct1 = [];
options_lree_ct1.stability = 1; % Basic
options_lree_ct1.multi = 1;
    
res_lree_ct1 = lree_nonlinear_ct1(nlid_data_ct1,options_lree_ct1);

[se_lree_ct1,ys_lree_ct1,xs_lree_ct1] = se_nonlinear_ct1(ud,y_raw,xd(:,1),lrm_ct1.e,lrm_ct1.f,lrm_ct1.g,res_lree_ct1.ec,res_lree_ct1.fc,res_lree_ct1.gc,dt);
res_lree_ct1.se = se_lree_ct1;
res_lree_ct1.dt = nlid_data_ct1.dt;

fprintf('Training error (lree ct1): %.5e\n',se_lree_ct1)

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_lree_ct1),ys_lree_ct1,'r--')
title(['LREE CT1 Training error: ' num2str(se_lree_ct1)])

end

%% Run LREE with CT1 and SOS

if 0

nlid_data_ct1 = nlid_data;
nlid_data_ct1.f_monos = res_prepro_ct1.f_monos;
nlid_data_ct1.g_monos = res_prepro_ct1.g_monos;  
nlid_data_ct1.dt = dt;

% Fit model with LR:
options_lree_ct1 = [];
options_lree_ct1.stability = 1; % Basic
options_lree_ct1.multi = 1e-4;
    
res_lree_sos_ct1 = lree_sos_nonlinear_ct1(nlid_data_ct1,options_lree_ct1);

[se_lree_sos_ct1,ys_lree_sos_ct1,xs_lree_sos_ct1] = se_nonlinear_ct1(ud,y_raw,xd(:,1),lrm_ct1.e,lrm_ct1.f,lrm_ct1.g,res_lree_sos_ct1.ec,res_lree_sos_ct1.fc,res_lree_sos_ct1.gc,dt);
res_lree_sos_ct1.se = se_lree_ct1;
res_lree_sos_ct1.dt = nlid_data_ct1.dt;

fprintf('Training error (lree sos ct1): %.5e\n',se_lree_sos_ct1)

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_lree_sos_ct1),ys_lree_sos_ct1,'r--')
title('LREE (SOS) CT1 (training)')

end

%% Simulate CT model with ODE45?

if 0

dt = nlid_data_ct1.dt; 
tspan = 0:dt:(dt*(length(ud)-1));
x0 = xd(:,1);

[t_sim, x_sim] = ode45(@(t,x) ode_ct1(x,ud,t,lrm_ct1,res_lree_ct1), tspan, x0);
x_sim = x_sim';

% To plot the output, run the states through the output map:
ys_ode = zeros(ny,length(x_sim));

for t = 1:length(x_sim)
    
    ys_ode(:,t) = lrm_ct1.g(x_sim(:,t),ud(:,t),res_lree_ct1.gc);
    
end

se_ct1_ode45 = nse(y_raw,ys_ode);

fprintf('Training error (ct1 ode45): %.5e\n',se_ct1_ode45)

figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_ode),ys_ode,'r--')
title('LREE ode45 (training)')

end



%% Run RIE   

if run_rie

options_rie.verbose = 0;

res_rie = localrie(nlid_data,options_rie);

[se_rie,ys_rie,xs_rie] = se_nonlinear_2(ud,y_raw,xd(:,1),lrm.e,lrm.f,lrm.g,res_rie.ec,res_rie.fc,res_rie.gc);

res_rie.se = se_rie;
res_rie.ys = ys_rie;
fprintf('\nSimulation error: %.4e\n',se_rie);    
        
% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_rie),ys_rie,'r--')

end


%% Run stable EE

if run_stblee

options_ee = []; 
options_ee.stability = 1; 
options_ee.verbose = 0; 
res_stblee = ee_nonlin_fast(nlid_data,options_ee);

[se_ee,ys_ee,xs_ee] = se_nonlinear_2(ud,yd,xd(:,1),lrm.e,lrm.f,lrm.g,res_stblee.ec,res_stblee.fc,res_stblee.gc);

res_stblee.se = se_ee;
res_stblee.ys = ys_ee;
fprintf('\nSimulation error (stblee): %.4e\n',se_ee); 

% Check output on training data
figure
plot(1:length(y_raw),y_raw)
hold on
plot(1:length(ys_ee),ys_ee,'r--')
title(['EE Training error: ' num2str(se_ee)])

end

