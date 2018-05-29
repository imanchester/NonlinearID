%% Validation: lree

% If possible, just `pick up' where you left off in terms of the initial state 

if run_lree
    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
res_tmp = res_lree{lrm_i};

lrm = lrms_lree{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_2(uv,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (lree): %.5e\n',ve_tmp)  


res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_lree{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['LREE Validation error: ' num2str(ve_tmp)])

end

end


%% Validation: lrse

if run_lrse
    
for lrm_i = lr_model_sel 

    
% Update the id_data structure with the appropriate model,
res_tmp = res_lrse{lrm_i};

lrm = lrms_lrse{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_2(uv,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (lrse): %.5e\n',ve_tmp)  

res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_lrse{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['LRSE Validation error: ' num2str(ve_tmp)])

end

end



%% Validation: lree ae

if run_lree_ae
    
for lrm_i = lr_model_sel 

if lrm_i ~= 1    
% Update the id_data structure with the appropriate model,
res_tmp = res_lree_ae{lrm_i};

lrm = lrms_lree_ae{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_ae(uv,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (lree ae): %.5e\n',ve_tmp)  

res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_lree_ae{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['LREE AE Validation error: ' num2str(ve_tmp)])

end

end

end



%% Validation: lree u2

if run_lree_u2
    
for lrm_i = lr_model_sel 
    
% Update the id_data structure with the appropriate model,
res_tmp = res_lree_u2{lrm_i};

lrm = lrms_lree_u2{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_2(uv2,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (lree u2): %.5e\n',ve_tmp)  

res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_lree_u2{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['LREE U2 Validation error: ' num2str(ve_tmp)])

end

end



%% Validation: lree u2

if run_stblee_u2
    
for lrm_i = lr_model_sel 
    
% Update the id_data structure with the appropriate model,
res_tmp = res_stblee_u2{lrm_i};

lrm = lrms_lree_u2{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_2(uv2,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (stblee u2): %.5e\n',ve_tmp)  

res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_stblee_u2{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['STBLEE U2 Validation error: ' num2str(ve_tmp)])

end

end





%% Validation: lree u2

if run_ee_u2
    
for lrm_i = lr_model_sel 
    
% Update the id_data structure with the appropriate model,
res_tmp = res_ee_u2{lrm_i};

lrm = lrms_lree_u2{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_2(uv2,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (ee u2): %.5e\n',ve_tmp)  

res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_ee_u2{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['EE U2 Validation error: ' num2str(ve_tmp)])

end

end



%% Validation: lrse u2

if run_lrse_u2
    
for lrm_i = lr_model_sel 
    
% Update the id_data structure with the appropriate model,
res_tmp = res_lrse_u2{lrm_i};

lrm = lrms_lree_u2{lrm_i};

[ve_tmp,yv_tmp,xv_tmp] = se_nonlinear_2(uv2,yv,res_tmp.xs(:,end),lrm.e,lrm.f,lrm.g,res_tmp.ec,res_tmp.fc,res_tmp.gc);

fprintf('Validation error (lrse u2): %.5e\n',ve_tmp)  

res_tmp.ve = ve_tmp;
res_tmp.yv = yv_tmp;
res_tmp.xv = xv_tmp;

res_lrse_u2{lrm_i} = res_tmp;

figure
plot(1:Tval,yv)
hold on
plot(1:Tval,yv_tmp,'r--')
title(['LRSE U2 Validation error: ' num2str(ve_tmp)])

end

end





%% Validation: lree (u2)

% if run_lree_u2
% 
% [ve_lree_u2,yv_lree_u2,xv_lree_u2] = se_nonlinear_2(uv2,yv,xs_lree_u2(:,end),lrm_u2.e,lrm_u2.f,lrm_u2.g,res_lree_u2.ec,res_lree_u2.fc,res_lree_u2.gc);
% 
% fprintf('Validation error (lree_u2): %.5e\n',ve_lree_u2)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_lree_u2,'r--')
% title(['LREE (U2) Validation error: ' num2str(ve_lree_u2)])
% 
% end
% 
% % Validation: lrse ae
% 
% if run_lree_ae
% 
% [ve_lree_ae,yv_lree_ae,xv_lree_ae] = se_nonlinear_ae(uv,yv,xs_lree_ae(:,end),lrm_ae.e,lrm_ae.f,lrm_ae.g,res_lree_ae.ec,res_lree_ae.fc,res_lree_ae.gc);
% 
% fprintf('Validation error (lree ae): %.5e\n',ve_lree_ae)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_lree_ae,'r--')
% title(['LREE AE Validation error: ' num2str(ve_lree_ae)])
% 
% end
% 
% 
% % Validation: ct1
% 
% if run_lree_ct1
%     
% [ve_lree_ct1,yv_lree_ct1,xv_lree_ct1] = se_nonlinear_ct1(uv,yv,xs_lree_ct1(:,end),lrm_ct1.e,lrm_ct1.f,lrm_ct1.g,res_lree_ct1.ec,res_lree_ct1.fc,res_lree_ct1.gc,dt);
% 
% fprintf('Validation error (lree ct1): %.5e\n',ve_lree_ct1)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_lree_ct1,'r--')
% title(['LREE CT1 Validation error: ' num2str(ve_lree_ct1)])
% 
% end
% 
% % Validation: lrse
% 
% % If possible, just `pick up' where you left off in terms of the initial state 
% 
% if run_lrse
% 
% [ve_lr,yv_lr,xv_lr] = se_nonlinear_2(uv,yv,xs_lr(:,end),lrm_se.e,lrm_se.f,lrm_se.g,res_lr.ec,res_lr.fc,res_lr.gc);
% 
% fprintf('Validation error (lrse): %.5e\n',ve_lr)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_lr,'r--')
% title(['LRSE Validation error: ' num2str(ve_lr)])
% 
% end
% 
% 
% %
%  
% if run_lrse_u2
% 
% [ve_lr_u2,yv_lr_u2,xv_lr_u2] = se_nonlinear_2(uv2,yv,xs_lr_u2(:,end),lrm_se_u2.e,lrm_se_u2.f,lrm_se_u2.g,res_lr_u2.ec,res_lr_u2.fc,res_lr_u2.gc);
% 
% fprintf('Validation error (lrse_u2): %.5e\n',ve_lr_u2)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_lr_u2,'r--')
% title(['LRSE_u2 Validation error: ' num2str(ve_lr_u2)])
% 
% end
% 
% 
% % Validation: rie
% 
% if run_rie
% 
% % If possible, just `pick up' where you left off in terms of the initial state 
% 
% [ve_rie,yv_rie,xv_rie] = se_nonlinear_2(uv,yv,xs_rie(:,end),lrm.e,lrm.f,lrm.g,res_rie.ec,res_rie.fc,res_rie.gc);
% 
% fprintf('Validation error (rie): %.5e\n',ve_rie)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_rie,'r--')
% 
% end
% 
% % Validation: ee
% 
% if run_stblee
%     
% % If possible, just `pick up' where you left off in terms of the initial state 
% 
% [ve_ee,yv_ee,xv_ee] = se_nonlinear_2(uv,yv,xs_ee(:,end),lrm.e,lrm.f,lrm.g,res_stblee.ec,res_stblee.fc,res_stblee.gc);
% 
% fprintf('Validation error (ee): %.5e\n',ve_ee)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_ee,'r--')
% title(['EE (stbl) Validation error: ' num2str(ve_ee)])
% 
% end
% 
% % Validation: ee (u2)
% 
% % if run_stblee_u2
% %     
% % [ve_ee_u2,yv_ee_u2,xv_ee_u2] = se_nonlinear_2(uv2,yv,xs_ee_u2(:,end),lrm_u2.e,lrm_u2.f,lrm_u2.g,res_stblee_u2.ec,res_stblee_u2.fc,res_stblee_u2.gc);
% % 
% % fprintf('Validation error (ee_u2): %.5e\n',ve_ee_u2)  
% % 
% % figure
% % plot(1:Tval,yv)
% % hold on
% % plot(1:Tval,yv_ee_u2,'r--')
% % title(['EE (stbl_u2) Validation error: ' num2str(ve_ee_u2)])
% % 
% % end