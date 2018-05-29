% -------------------------------------------------------------------------
%% Models for LREE
% -------------------------------------------------------------------------

gx_deg = 1;

clear prg
clear spot_monos_e spot_monos_f spot_monos_g

nu = 1;
ny = 1;

prg = spotsosprog;
[prg,x] = prg.newIndeterminate('x',nx);
[prg,u] = prg.newIndeterminate('u',nu);

%% M (1,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 1)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:1)'; 
    spot_monos_e{2} = monomials(x,1:1)';
    spot_monos_e{3} = monomials(x,1:1)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree{1} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (3,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 2)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:3)'; 
    spot_monos_e{2} = monomials(x,1:3)';
    spot_monos_e{3} = monomials(x,1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree{2} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (3,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 3)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:3)'; 
    spot_monos_e{2} = monomials(x,1:3)';
    spot_monos_e{3} = monomials(x,1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree{3} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (5,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 4)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:5)'; 
    spot_monos_e{2} = monomials(x,1:5)';
    spot_monos_e{3} = monomials(x,1:5)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree{4} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end



% -------------------------------------------------------------------------
%% Models for LRSE
% -------------------------------------------------------------------------

clear prg
clear spot_monos_e spot_monos_f spot_monos_g

nu = 1;
ny = 1;

prg = spotsosprog;
[prg,x] = prg.newIndeterminate('x',nx);
[prg,u] = prg.newIndeterminate('u',nu);


%% M (1,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 1)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:1)'; 
    spot_monos_e{2} = monomials(x,1:1)';
    spot_monos_e{3} = monomials(x,1:1)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lrse{1} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end

%% M (3,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 2)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:3)'; 
    spot_monos_e{2} = monomials(x,1:3)';
    spot_monos_e{3} = monomials(x,1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lrse{2} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (3,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 3)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:3)'; 
    spot_monos_e{2} = monomials(x,1:3)';
    spot_monos_e{3} = monomials(x,1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lrse{3} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (5,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 4)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:5)'; 
    spot_monos_e{2} = monomials(x,1:5)';
    spot_monos_e{3} = monomials(x,1:5)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lrse{4} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end



% -------------------------------------------------------------------------
%% Models for LREE (alternate e)
% -------------------------------------------------------------------------

if run_lree_ae

clear prg
clear spot_monos_e spot_monos_f spot_monos_g

nu = 1;
ny = 1;

prg = spotsosprog;
[prg,x] = prg.newIndeterminate('x',nx);
[prg,u] = prg.newIndeterminate('u',nu);


%% M (1,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 1)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials([x;u],1:1)'; 
    spot_monos_e{2} = monomials([x;u],1:1)';
    spot_monos_e{3} = monomials([x;u],1:1)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_ae{1} = buildModel_ae(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (3,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 2)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials([x;u],1:3)'; 
    spot_monos_e{2} = monomials([x;u],1:3)';
    spot_monos_e{3} = monomials([x;u],1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_ae{2} = buildModel_ae(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (3,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 3)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials([x;u],1:3)'; 
    spot_monos_e{2} = monomials([x;u],1:3)';
    spot_monos_e{3} = monomials([x;u],1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_ae{3} = buildModel_ae(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (5,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 4)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials([x;u],1:5)'; 
    spot_monos_e{2} = monomials([x;u],1:5)';
    spot_monos_e{3} = monomials([x;u],1:5)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_ae{4} = buildModel_ae(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


end


% -------------------------------------------------------------------------
%% Models for LREE U2
% -------------------------------------------------------------------------

clear prg
clear spot_monos_e spot_monos_f spot_monos_g

% nx = 3;
nu = 2;
% ny = 1;

prg = spotsosprog;
[prg,x] = prg.newIndeterminate('x',nx);
[prg,u] = prg.newIndeterminate('u',nu);


%% M (1,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 1)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:1)'; 
    spot_monos_e{2} = monomials(x,1:1)';
    spot_monos_e{3} = monomials(x,1:1)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_u2{1} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end

%% M (3,1,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 2)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:3)'; 
    spot_monos_e{2} = monomials(x,1:3)';
    spot_monos_e{3} = monomials(x,1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:1)'];
    spot_monos_f{2} = [monomials([x;u],1:1)'];
    spot_monos_f{3} = [monomials([x;u],1:1)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_u2{2} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (3,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 3)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:3)'; 
    spot_monos_e{2} = monomials(x,1:3)';
    spot_monos_e{3} = monomials(x,1:3)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_u2{3} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end


%% M (5,3,1)
% -------------------------------------------------------------------------

if sum(lr_model_sel == 4)


if narx_structure   
    
    
    spot_monos_g{1} = [x(1)];

    
else
    
    spot_monos_e{1} = monomials(x,1:5)'; 
    spot_monos_e{2} = monomials(x,1:5)';
    spot_monos_e{3} = monomials(x,1:5)';    
    
    spot_monos_f{1} = [monomials([x;u],1:3)'];
    spot_monos_f{2} = [monomials([x;u],1:3)'];
    spot_monos_f{3} = [monomials([x;u],1:3)'];    
    
    spot_monos_g{1} = [monomials(u,1)'  monomials(x,1:gx_deg)'];

    
end

lrms_lree_u2{4} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);

end




% %% Model: alternate E
% % -------------------------------------------------------------------------
% 
% % Monomials in e
% % -------------------------------------------------------------------------
% % spot_monos_e_ae{1} = [monomials([x;u],1:3)' monomials(x,4:5)']; 
% 
% if narx_structure
% 
%     spot_monos_e_ae{1} = [monomials([x;u],1:e_degree)']; 
%     spot_monos_e_ae{2} = [monomials([x;u],1)']; 
%     spot_monos_e_ae{3} = [monomials([x;u],1)'];
% 
% else
%     
%     spot_monos_e_ae{1} = [monomials([x;u],1:e_degree)']; 
%     spot_monos_e_ae{2} = [monomials([x;u],1:e_degree)']; 
%     spot_monos_e_ae{3} = [monomials([x;u],1:e_degree)'];   
%     
% end
% % spot_monos_e_ae{1} = [monomials([x;u],1:5)']; 
% % spot_monos_e_ae{2} = [monomials([x;u],1:5)']; 
% % spot_monos_e_ae{3} = [monomials([x;u],1:5)'];
% 
% % Monomials in f
% % -------------------------------------------------------------------------
% if narx_structure
% 
%     spot_monos_f_ae{1} = [monomials([x;u],1:f_degree)']; 
%     spot_monos_f_ae{2} = [monomials([x],1)'];
%     spot_monos_f_ae{3} = [monomials([x],1)'];
%     
% else
%      
%     spot_monos_f_ae{1} = [monomials([x;u],1:f_degree)']; 
%     spot_monos_f_ae{2} = [monomials([x;u],1:f_degree)'];
%     spot_monos_f_ae{3} = [monomials([x;u],1:f_degree)'];
%     
% end
% 
% % Monomials in g
% % -------------------------------------------------------------------------
% if narx_structure
%     spot_monos_g_ae{1} = [x(1)];
% else
%     spot_monos_g_ae{1} = [monomials(u,1:3)'  monomials(x,1:g_degree)'];
% end
% 
% res_prepro_ae = buildModel_ae(prg,x,u,spot_monos_e_ae,spot_monos_f_ae,spot_monos_g_ae);
% 
% e_monos_ae = res_prepro_ae.e_monos;
% f_monos_ae = res_prepro_ae.f_monos;
% g_monos_ae = res_prepro_ae.g_monos;
% 
% e_id_ae = res_prepro_ae.e;
% f_id_ae = res_prepro_ae.f;
% g_id_ae = res_prepro_ae.g;
% 
% lrm_ae.e = e_id_ae;
% lrm_ae.f = f_id_ae;
% lrm_ae.g = g_id_ae;


% %% Model: CT1
% 
% spot_monos_f_ct1{1} = [monomials([x;u],1:5)']; 
% spot_monos_f_ct1{2} = [monomials([x],1)'];
% spot_monos_f_ct1{3} = [monomials([x],1)'];
% 
% res_prepro_ct1 = buildModel_ct1(prg,x,u,spot_monos_f_ct1,spot_monos_g_ae);
% 
% lrm_ct1.e = res_prepro_ct1.e;
% lrm_ct1.f = res_prepro_ct1.f;
% lrm_ct1.g = res_prepro_ct1.g;
% 
% 
% %% What happens if you add an additional input 
% 
% clear prg
% clear spot_monos_e spot_monos_f spot_monos_g
% 
% nx = 3;
% nu = 2;
% % ny = 1;
% 
% prg = spotsosprog;
% [prg,x] = prg.newIndeterminate('x',nx);
% [prg,u] = prg.newIndeterminate('u',nu);
% 
% % Monomials in e
% % -------------------------------------------------------------------------
% % spot_monos_e{1} = monomials(x,1:5)'; 
% % spot_monos_e{2} = monomials(x,1)';
% % spot_monos_e{3} = monomials(x,1)';
% 
% if narx_structure
%     
% else
%     spot_monos_e{1} = monomials(x,1:e_degree)'; 
%     spot_monos_e{2} = monomials(x,1:e_degree)';
%     spot_monos_e{3} = monomials(x,1:e_degree)';
% end
% 
% % Monomials in f
% % -------------------------------------------------------------------------
% % spot_monos_f{1} = [monomials(u,1:3)' monomials([x],1:3)'];
% % spot_monos_f{2} = [monomials([x],1)'];
% % spot_monos_f{3} = [monomials([x],1)'];
% 
% if narx_structure
% 
% else
% % spot_monos_f{1} = [monomials(u,1:3)' monomials([x],1:3)'];
% % spot_monos_f{2} = [monomials(u,1:3)' monomials([x],1:3)'];
% % spot_monos_f{3} = [monomials(u,1:3)' monomials([x],1:3)'];
% 
%     spot_monos_f{1} = [monomials([x;u],1:f_degree)'];
%     spot_monos_f{2} = [monomials([x;u],1:f_degree)'];
%     spot_monos_f{3} = [monomials([x;u],1:f_degree)'];
% 
% end
% 
% % Monomials in g
% % -------------------------------------------------------------------------
% if narx_structure
%     spot_monos_g{1} = [x(1)];
% else
%     spot_monos_g{1} = [monomials(u,1:3)'  monomials(x,1:g_degree)'];
% end
% 
% % Preprocessing of models
% 
% res_prepro_u2 = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);
% 
% e_monos_u2 = res_prepro_u2.e_monos;
% f_monos_u2 = res_prepro_u2.f_monos;
% g_monos_u2 = res_prepro_u2.g_monos;
% 
% lrm_u2.e = res_prepro_u2.e;
% lrm_u2.f = res_prepro_u2.f;
% lrm_u2.g = res_prepro_u2.g;
% 
% nlid_data_u2.e_monos = e_monos_u2;
% nlid_data_u2.f_monos = f_monos_u2;
% nlid_data_u2.g_monos = g_monos_u2;
% 
% ud2 = [ud;ud(2:end),0];
% 
% nlid_data_u2.x = xd;
% nlid_data_u2.u = ud2;
% nlid_data_u2.y = yd; 
% 
% 
% %% additional model for lrse u2
% 
% % Monomials in e
% % -------------------------------------------------------------------------
% % spot_monos_e{1} = monomials(x,1:5)'; 
% % spot_monos_e{2} = monomials(x,1)';
% % spot_monos_e{3} = monomials(x,1)';
% 
% if narx_structure
%     
% else
%     spot_monos_e{1} = monomials(x,1:e_degree)'; 
%     spot_monos_e{2} = monomials(x,1:e_degree)';
%     spot_monos_e{3} = monomials(x,1:e_degree)';
% end
% 
% % Monomials in f
% % -------------------------------------------------------------------------
% % spot_monos_f{1} = [monomials(u,1:3)' monomials([x],1:3)'];
% % spot_monos_f{2} = [monomials([x],1)'];
% % spot_monos_f{3} = [monomials([x],1)'];
% 
% if narx_structure
% 
% else
% % spot_monos_f{1} = [monomials(u,1:3)' monomials([x],1:3)'];
% % spot_monos_f{2} = [monomials(u,1:3)' monomials([x],1:3)'];
% % spot_monos_f{3} = [monomials(u,1:3)' monomials([x],1:3)'];
% 
%     spot_monos_f{1} = [monomials([x;u],1:f_degree)'];
%     spot_monos_f{2} = [monomials([x;u],1:f_degree)'];
%     spot_monos_f{3} = [monomials([x;u],1:f_degree)'];
% 
% end
% 
% % Monomials in g
% % -------------------------------------------------------------------------
% if narx_structure
%     spot_monos_g{1} = [x(1)];
% else
%     spot_monos_g{1} = [monomials(u,1:3)'  monomials(x,1:g_degree)'];
% end
% 
% % Preprocessing of models
% 
% res_prepro_se_u2 = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);
% 
% e_monos_u2 = res_prepro_se_u2.e_monos;
% f_monos_u2 = res_prepro_se_u2.f_monos;
% g_monos_u2 = res_prepro_se_u2.g_monos;
% 
% lrm_se_u2.e = res_prepro_se_u2.e;
% lrm_se_u2.f = res_prepro_se_u2.f;
% lrm_se_u2.g = res_prepro_se_u2.g;
% 
% nlid_data_se_u2 = nlid_data_u2;
% 
% nlid_data_se_u2.e_monos = e_monos_u2;
% nlid_data_se_u2.f_monos = f_monos_u2;
% nlid_data_se_u2.g_monos = g_monos_u2;