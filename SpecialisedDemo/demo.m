%% Demonstrate specialized algorithm

% Fit a stable nonlinear dynamical system to data from a simulated
% mass-spring-damper system.

% Please refer to readme.txt for dependencies.

close all
clear variables
clc

%% True system

% Mass - nonlinear spring - damper

% Dimensions
nx = 4;
ny = 1;
nu = 1;
nw = 1;

param.idt = 0.1; 

param.m = [0.5,0.1];        % Mass
param.c = [0.5,0.75];       % Damping
param.k = [2,1];            % Spring 
param.nl = [1,1];           % Spring nonlinearity type
param.slim = [1.25,1.25];   % Spring travel limits

nlsprng_models{1} = param;

%% Experiment parameters

% Number of data points used for fitting,
datapoints = [500]; 

% No. of repetitions per model
numreps = 1;

% Methods to run
run_lrse = 1; % Specialized method for Lagrangian relaxation of simulation error


%% Polynomial models

% Each model structure goes into a cell-array called 'lrms'

clear prg

nx = 4;
nu = 1;
ny = 1;

prg = spotsosprog;
[prg,x] = prg.newIndeterminate('x',nx);
[prg,u] = prg.newIndeterminate('u',nu);

% Model 1 (3rd order)
% -------------------------------------------------------------------------

% Monomials in e
spot_monos_e{1} = monomials(x,1:3)'; 
spot_monos_e{2} = monomials(x,1:3)'; 
spot_monos_e{3} = monomials(x,1:3)'; 
spot_monos_e{4} = monomials(x,1:3)'; 

% Monomials in f
spot_monos_f{1} = [u(1) monomials(x,1:3)'];
spot_monos_f{2} = [u(1) monomials(x,1:3)'];
spot_monos_f{3} = [u(1) monomials(x,1:3)'];
spot_monos_f{4} = [u(1) monomials(x,1:3)'];


% Monomials in g
spot_monos_g{1} = [u(1) x(1) x(2) x(3) x(4)];

% Process model
lrms{1} = buildModel_2(prg,x,u,spot_monos_e,spot_monos_f,spot_monos_g);


%% Store preamble

% Select models, from lrms, to fit.
lr_model_sel = [1];

preamble.datapoints = datapoints;
preamble.numreps = numreps;
preamble.lr_models_sel = lr_model_sel;
preamble.lrms = lrms;

rslts.preamble = preamble;
rslts.results = [];

% Generate fresh file name:
s_date = datestr(clock,'ddmmmyyyy_HHMM');
% s_name = ['timing_res_' s_date];
s_name = ['demo_results_' s_date];
s_name_t = [s_name '.mat'];

save(s_name_t,'rslts','-v7.3')

clear preamble

%% Run the external loop 


for tindex = 1:length(datapoints)
    
       
    for rep = 1:numreps
        
        fprintf('\n--------------------------------------------\n')
        fprintf('--------------------------------------------\n')
        fprintf('*** Rep = %d, Time: %s ***\n',rep,datestr(clock,'HH:MM'))
        fprintf('--------------------------------------------\n')
        fprintf('--------------------------------------------\n\n')
        

 
% Generate problem data
% -------------------------------------------------------------------------
% Randomly generate an input signal

        ndp = datapoints(tindex); % No. of data points
        dt = param.idt;  % Sampling time
        nseg = 5; 

        ud = [];

        for i = 1:nseg

            tt = 0:dt:(ndp*dt/nseg); % Time interval

            ts = (round(1.0*rand)*(2 + 0.75*rand))*sin(2*pi/(12+5*rand)*tt + pi/2+pi/4*rand) + ...
                 (round(0.8*rand)*(1 + 0.5*rand))*sin(2*pi/(19 + 8*rand)*tt + pi/2+pi/2*rand) + ...
                 (1 + 0.75*rand)*sin(2*pi/(17+10*rand)*tt + pi/2+pi/4*rand) + ...
                 (round(0.95*rand)*(1 + 0.5*rand))*sin(2*pi/(8+5*rand)*tt + pi/2+pi/2*rand) + ... 
                 (0.75 + 0.2*rand)*sin(2*pi/(5+2*rand)*tt + pi/2+pi/4*rand) + ...
                 (round(rand)*(1 + 0.5*rand))*sin(2*pi/(4+1.5*rand)*tt + pi/2+pi/4*rand) + ...
                 (1*round(rand))*sin(2*pi/(4.1+3*rand)*tt + pi/2+pi/4*rand);     

            ud = [ud,ts];
        end


% Simulation
% -------------------------------------------------------------------------

        dt = param.idt; 
        tspan = 0:dt:(ndp*dt);
        x0 = zeros(nx,1);

        [t_sim, x_sim] = ode45(@(t,x) msd_nonlinearSpring_gen(x,ud,t,param), tspan, x0);

        x_sim = x_sim';
        y_sim = x_sim(2,:);

% Add some noise to the simulated quantities:
        xn = x_sim + 1*mvnrnd(zeros(size(x_sim,1),1),1e-4*eye(size(x_sim,1)),length(x_sim))';
        yn = xn(2,:); 

% Scale the quantities for better numerical properties
        ssf = diag(max(abs(xn),[],2));    % State
        osf = diag(max(abs(yn),[],2));    % Output
        xsc = ssf\xn;
        ysc = osf\yn;
        
% Filter the outputs for state estimates        
        [xd,yd,ud] = processStates4nlmsd(xsc,ysc,ud,dt);


% Set and save the problem data
% -------------------------------------------------------------------------

% Structure to pass problem data to identification algorithms
        nlid_data.x = xd;
        nlid_data.u = ud;
        nlid_data.y = yd;

        % Signal to noise ratio:
        snr_emp = var(y_sim)/var(y_sim-yn);
        fprintf('SNR: %.5e (%.5edB)\n',snr_emp,10*log10(snr_emp))

        % Save
        load(s_name_t)
        rslts.results{rep}.data = nlid_data;
        rslts.results{rep}.ssf = ssf;
        rslts.results{rep}.osf = osf;
        save(s_name_t,'rslts','-v7.3')

        clear rslts 
        clear x_sim y_sim xn yn 
    

% Loop over different models
% -------------------------------------------------------------------------
 
        for lrm_i = lr_model_sel 


        % Update the id_data structure with the appropriate model,
            lrm = lrms{lrm_i};

            nlid_data.e_monos = lrm.e_monos;
            nlid_data.f_monos = lrm.f_monos;
            nlid_data.g_monos = lrm.g_monos;    


%% Run Specialized
% -------------------------------------------------------------------------

            if run_lrse

                fprintf('\n\n--------------------------------------------\n')            
                fprintf('\n\nRunning Specialized: model = %d, Time: %s ***\n',lrm_i,datestr(clock,'HH:MM'))
                fprintf('--------------------------------------------\n\n')

% Fit model with LR:
                options_specialized = [];
                options_specialized.verbose = 0;   % Minimal output to terminal
                options_specialized.trcnst = 0;    % Impose limit on the trace? No
                options_specialized.earlyterm = 1; % Terminate inner loop early? Yes
                
% Run specialized algorithm
                res_lr = lrse_nonlin_specialized(nlid_data,options_specialized);

% Compute simulation error of identified model (on scaled problem data)                 
                [se_lr,ys_lr] = se_nonlinear_2(ud,yd,xd(:,1),lrm.e,lrm.f,lrm.g,res_lr.ec,res_lr.fc,res_lr.gc);
                res_lr.se = se_lr;

% Report anything interesting:
                fprintf('\nSimulation error: %.4e\n',se_lr); 
                fprintf('Time (LR): %.5e\n',res_lr.solvertime)
                fprintf('Bound (LR): %.9e\n',res_lr.Jf)

                load(s_name_t)
                rslts.results{rep}.lr{lrm_i} = res_lr;
                save(s_name_t,'rslts','-v7.3')

                clear rslts res_lr_red ys_lr options_specialized res_lr

            end


        end
       

    end
    
    
end

 


