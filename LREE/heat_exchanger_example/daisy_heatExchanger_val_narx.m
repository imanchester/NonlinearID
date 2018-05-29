%%

if run_poly_p1
    
    res_tmp = res_narx_poly_p{1};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
    fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  

    res_tmp.ve = ve_narx_poly;
    res_tmp.yv = yv_narx_poly;
    
    res_narx_poly_p{1} = res_tmp;
    
    figure
    plot(1:Tval,yv)
    hold on
    plot(1:Tval,yv_narx_poly,'r--')
    ax=gca;
    set(ax,'YLim',[-1,1])
    title(['narx (poly) Validation error: ' num2str(ve_narx_poly)])

end


if run_poly_p2
    
    res_tmp = res_narx_poly_p{2};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
    fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  

    res_tmp.ve = ve_narx_poly;
    res_tmp.yv = yv_narx_poly;
    
    res_narx_poly_p{2} = res_tmp;
    
    figure
    plot(1:Tval,yv)
    hold on
    plot(1:Tval,yv_narx_poly,'r--')
    ax=gca;
    set(ax,'YLim',[-1,1])
    title(['narx (poly) Validation error: ' num2str(ve_narx_poly)])

end



if run_poly_p3
    
    res_tmp = res_narx_poly_p{3};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
    fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  

    res_tmp.ve = ve_narx_poly;
    res_tmp.yv = yv_narx_poly;
    
    res_narx_poly_p{3} = res_tmp;
    
    figure
    plot(1:Tval,yv)
    hold on
    plot(1:Tval,yv_narx_poly,'r--')
    ax=gca;
    set(ax,'YLim',[-1,1])
    title(['narx (poly) Validation error: ' num2str(ve_narx_poly)])

end

%%

if run_poly_s1
    
    res_tmp = res_narx_poly_s{1};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
    fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  

    res_tmp.ve = ve_narx_poly;
    res_tmp.yv = yv_narx_poly;
    
    res_narx_poly_s{1} = res_tmp;
    
    figure
    plot(1:Tval,yv)
    hold on
    plot(1:Tval,yv_narx_poly,'r--')
    ax=gca;
    set(ax,'YLim',[-1,1])
    title(['narx (poly) Validation error: ' num2str(ve_narx_poly)])

end




if run_poly_s2
    
    res_tmp = res_narx_poly_s{2};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
    fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  

    res_tmp.ve = ve_narx_poly;
    res_tmp.yv = yv_narx_poly;
    
    res_narx_poly_s{2} = res_tmp;
    
    figure
    plot(1:Tval,yv)
    hold on
    plot(1:Tval,yv_narx_poly,'r--')
    ax=gca;
    set(ax,'YLim',[-1,1])
    title(['narx (poly) Validation error: ' num2str(ve_narx_poly)])

end



if run_poly_s3
    
    res_tmp = res_narx_poly_s{3};
    
    nlarx_poly = res_tmp.nlarx;

    figure(101)
    nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
    close(101)

    yv_narx_poly = get(nlarx_poly_val,'outputData')';

    ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
    fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  

    res_tmp.ve = ve_narx_poly;
    res_tmp.yv = yv_narx_poly;
    
    res_narx_poly_s{3} = res_tmp;
    
    figure
    plot(1:Tval,yv)
    hold on
    plot(1:Tval,yv_narx_poly,'r--')
    ax=gca;
    set(ax,'YLim',[-1,1])
    title(['narx (poly) Validation error: ' num2str(ve_narx_poly)])

end









% %% Validation: sigmoid p
% 
% figure(101)
% nlarx_val = compare(narx_dataset_val,nlarx_sig_p);
% close(101)
% 
% yv_narx = get(nlarx_val,'outputData')';
% 
% ve_narx_sig_p = nse(yv(:,:),yv_narx(1,:));
% fprintf('Validation error (sig): %.5e\n',ve_narx_sig_p)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx,'r--')
% title(['narx (sig) Validation error: ' num2str(ve_narx_sig_p)])
% 
% % Validation: sigmoid s
% 
% figure(101)
% nlarx_val = compare(narx_dataset_val,nlarx_sig_s);
% close(101)
% 
% yv_narx = get(nlarx_val,'outputData')';
% 
% ve_narx_sig_s = nse(yv(:,:),yv_narx(1,:));
% fprintf('Validation error (sig s): %.5e\n',ve_narx_sig_s)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx,'r--')
% title(['narx (sig s) Validation error: ' num2str(ve_narx_sig_s)])
% 
% 
% % Validation: wave p
% 
% figure(101)
% nlarx_wav_val = compare(narx_dataset_val,nlarx_wav_p);
% close(101)
% 
% yv_narx_wav_p = get(nlarx_wav_val,'outputData')';
% 
% ve_narx_wav_p = nse(yv(:,:),yv_narx_wav_p(1,:));
% fprintf('Validation error (wav): %.5e\n',ve_narx_wav_p)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx_wav,'r--')
% title(['narx (wav) Validation error: ' num2str(ve_narx_wav_p)])
% 
% % Validation: wave s
% 
% figure(101)
% nlarx_wav_val = compare(narx_dataset_val,nlarx_wav_s);
% close(101)
% 
% yv_narx_wav_s = get(nlarx_wav_val,'outputData')';
% 
% ve_narx_wav_s = nse(yv(:,:),yv_narx_wav_s(1,:));
% fprintf('Validation error (wav s): %.5e\n',ve_narx_wav_s)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx_wav,'r--')
% title(['narx (wav s) Validation error: ' num2str(ve_narx_wav_s)])
% 
% % Validation: sigmoid ps
% 
% if run_sig_ps
% 
% figure(101)
% nlarx_val = compare(narx_dataset_val,nlarx_sig_ps);
% close(101)
% 
% yv_narx = get(nlarx_val,'outputData')';
% 
% ve_narx_sig_ps = nse(yv(:,:),yv_narx(1,:));
% fprintf('Validation error (sig ps): %.5e\n',ve_narx_sig_ps)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx,'r--')
% title(['narx (sig ps) Validation error: ' num2str(ve_narx_sig_ps)])
% 
% end
% 
% % Validation: wave ps
% 
% if run_wav_ps
%     
% figure(101)
% nlarx_wav_val = compare(narx_dataset_val,nlarx_wav_ps);
% close(101)
% 
% yv_narx_wav = get(nlarx_wav_val,'outputData')';
% 
% ve_narx_wav_ps = nse(yv(:,:),yv_narx_wav(1,:));
% fprintf('Validation error (wav ps): %.5e\n',ve_narx_wav_ps)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx_wav,'r--')
% title(['narx (wav ps) Validation error: ' num2str(ve_narx_wav_ps)])
% 
% end
% 
% % Validation: sigmoid ss
% 
% if run_sig_ss
%     
% figure(101)
% nlarx_val = compare(narx_dataset_val,nlarx_sig_ss);
% close(101)
% 
% yv_narx = get(nlarx_val,'outputData')';
% 
% yv_narx_sig_ss = yv_narx;
% 
% ve_narx_sig_ss = nse(yv(:,:),yv_narx(1,:));
% fprintf('Validation error (sig ss): %.5e\n',ve_narx_sig_ss)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx,'r--')
% title(['narx (sig ss) Validation error: ' num2str(ve_narx_sig_ss)])
% 
% end
% 
% % Validation: wave ss
% 
% if run_wav_ss
%        
% figure(101)
% nlarx_wav_val = compare(narx_dataset_val,nlarx_wav_ss);
% close(101)
% 
% yv_narx_wav = get(nlarx_wav_val,'outputData')';
% 
% ve_narx_wav_ss = nse(yv(:,:),yv_narx_wav(1,:));
% fprintf('Validation error (wav ss): %.5e\n',ve_narx_wav_ss)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx_wav,'r--')
% title(['narx (wav ss) Validation error: ' num2str(ve_narx_wav_ss)])
% 
% end
% 
% % Validation: poly
% 
% figure(101)
% nlarx_poly_val = compare(narx_dataset_val,nlarx_poly);
% close(101)
% 
% yv_narx_poly = get(nlarx_poly_val,'outputData')';
% 
% ve_narx_poly = nse(yv(:,:),yv_narx_poly(1,:));
% fprintf('Validation error (poly): %.5e\n',ve_narx_poly)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx_poly,'r--')
% ax=gca;
% set(ax,'YLim',[-1,1])
% title(['narx (wav) Validation error: ' num2str(ve_narx_poly)])
% 
% % Validation: poly p
% 
% figure(101)
% nlarx_poly_val = compare(narx_dataset_val,nlarx_poly_p);
% close(101)
% 
% yv_narx_poly_p = get(nlarx_poly_val,'outputData')';
% 
% ve_narx_poly_p = nse(yv(:,:),yv_narx_poly_p(1,:));
% fprintf('Validation error (poly p): %.5e\n',ve_narx_poly_p)  
% 
% figure
% plot(1:Tval,yv)
% hold on
% plot(1:Tval,yv_narx_poly_p,'r--')
% ax=gca;
% set(ax,'YLim',[-1,1])
% title(['narx (poly) Validation error: ' num2str(ve_narx_poly_p)])





