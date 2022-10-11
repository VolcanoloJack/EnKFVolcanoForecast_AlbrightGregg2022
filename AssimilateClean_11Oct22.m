function [] = AssimilateClean_11Oct22(savename,N,assim,infl,init,...
    data,model,Dload,iniskip)
% Assimilates synthetic GPS and InSAR data

%% Setup

% Inputs
% savename - file name for saving results (string)
% N - number of ensemble members
% assim - assimilation algorithm used (1 = pseudoinverse, 2 = low rank,
%           3 = rank reduction, 4 = slow convergence rank reduction,
%           5 = square root, 6 = square root random roatation, 7 = DEnKF)
% infl - inflation algorithm used (0 = none, 1 = flat, 2 = individual, 
%           3 = all w/ limit, 4 = all, 5 = constant)
% init - method for initiating ensemble (1 = flat, 2 = gaussian, 3 =
%           supersample)
% data - data set used (1 = P, 2 = R)
% model - predictive model (1 = Yang, 2 = 2.5D FEM)
% Dload - scaling of diagonal loading for some assimilation methods
% iniskip - [0,1] skip initial time steps with near-0 data

% Assimilation Options
p = 6; % Number of model parameters: dV, gamma, psi, delta, dX, dY
drv = 7; % Number of derived values: V, theta, r1, r2, dP, dZ, ERR
gps_reps = 1; % Number of iterations per GPS step
insar_reps = 1; % Number of iterations per InSAR step
saveloc = 'Outputs'; % Subfolder to save results
flattol = 2.5; % threshold for mean(abs(d)./e) to end initial skipping [2.5]
isfilttag = 0; % toggle for filtering near-0 insar data
kinfl = 1.05; % constant inflation factor (when applicable)

% Scaling Options
wth = 0.05; wset = 0.1; wmax = 0.25; % Threshold, goal, and limit stdevs
qtol = 5; % Tolerance on interquartile outlier filter

% Misc Options
rng('shuffle'); RNGsus = rng;

addpath(genpath('YangEllipsoid'))

%% Define Constants

% Material Constants
E = 50e9; % Young's Modulus [Pa]
nu = 0.25; % Poisson's Ratio
G = E/(2*(1+nu)); % Shear Modulus [Pa]
K = E/(3*(1-2*nu)); % Bulk Modulus [Pa]
lambda = (E*nu)/((1+nu)*(1-2*nu)); % Lame Parameter [Pa]

% Limits and Tuning Factors
eps = [300, 300, 500]; % Minima for r1, r2, and depth to crest
sig = [7750 7750 11000]; % eps+sig = maxima for r1, r2, and dZ
k = [5e-2,5e-2,5e-2]; % Slope factors for logistic functions
o = [50,120,30]; % offset factors for logistic functions
Vr = (((4*pi())/3)*1000^3*1e6)/K; % Reference volume change (1 km at 1 MPa)

pwiden = [5 5 5 5 150 150]; % Alternative thresholds for inflation

% Time Constants
ysec = 60*60*24*365.25; % seconds in a year

% Geometric Constants
theta = 89.99 * (pi/180); % Plunge of reservoir for Yang (cannot be 90 deg)

%% Load Data
if data == 1
    load('Synth_P_06Mar21.mat')
else
    load('Synth_R_06Mar21.mat')
end
t0 = INSAR{1,1}(1);

synths0 = [r1s;r2s;dPs;dXs;dYs;dZs];
[synths,synthsB] = EnKFParamShift(synths0,eps,sig,Vr,k,o,K,0);

synths = [synths;synthsB;synths0([1,2,3,6],:)];

%% Define Time Steps

insar = size(INSAR,1); gps = size(GPS,1);
gps_steps = []; insar_steps = [];

for i = 1:gps
    gps_steps = union(gps_steps,GPS{i,2});
end

for i = 1:insar
    insar_steps = union(insar_steps,(INSAR{i,1}(1)-t0)*ysec);
end

steps = union(gps_steps,insar_steps);
reps = length(gps_steps)*gps_reps + length(insar_steps)*insar_reps;
tend = round(((steps(end)/ysec)+t0),1);

%% Initialize Ensemble

% Setup A tracking
Alog = zeros(p,N,reps*5+2); % Keeps full matrix of parameters
Astep = zeros(4,reps*5+2); % Records which step and data match each A
Amean = zeros(p,reps*5+2); % Mean of each parameter at each time step
Astd = zeros(p,reps*5+2); % Standard deviations
Asynths = zeros(p+drv-1,reps*5+2); % Synthetic values at each time step
datlog = cell(3,length(gps_steps) + length(insar_steps)); % Stores data
Acount = 1; datcount = 0;

% For Astep: 0 = initial, 1 = EnKF, 2 = outliers, 3 = widen, 4 = impossible,      
%   5 = shuffle

% Initial Parameter Estimates and Standard Deviations
est = zeros(1,6); devs = zeros(1,6);
est(1) = 0;  devs(1) = 70;    % dV - Normalized Elastic Volume Change
est(2) = 0;  devs(2) = 40;    % gamma - Aspect Ratio Scale
est(3) = 0;  devs(3) = 40;    % psi - Volume/Pressure partition
est(4) = 0;  devs(4) = 40;    % delta - Depth partition
est(5) = 0;  devs(5) = 2000;  % dX - E-W Offset (m)
est(6) = 0;  devs(6) = 2000;  % dY - N-S Offset (m)

% Build A matrix
A = ones(6,N);

if init ~= 3 
    for i = 1:p
        if init == 1 % Flat Distribution
            Q = rand(1,N);
        elseif init == 2 % Normal Distribution
            Q = randn(1,N);
        end
        Q = Q - mean(Q); Q = (devs(i)/std(Q)) * Q;
        A(i,:) = A(i,:)*est(i) + Q;
    end
else % Supersample from Evansen (2004)
    % Generate Larger Ensemble
    beta = 3; A = ones(6,N*beta);
    for i = 1:p
        Q = randn(1,N*beta); Q = Q - mean(Q);
        A(i,:) = A(i,:)*est(i) + devs(i)*Q;
    end

    % Forward Calculation to Sample Model State
    xy = zeros(gps,2);
    for i = 1:gps
        xy(i,:) = GPS{i,1};
    end
    xy = [xy ; INSAR{end,2}]; A = [A ; ones(length(xy)+2*gps,N*beta)];
    look = INSAR{end,1}(2); track = INSAR{end,1}(3);
    
    [Ap,~] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);

    for i = 1:N*beta
        [Ux,Uy,Uz] = yangdisp(Ap(4,i),Ap(5,i),Ap(6,i),...
            Ap(1,i),Ap(2,i),lambda,G,nu,Ap(3,i),theta,0,...
            xy(:,1),xy(:,2),xy(:,1)*0);
        Ux = real(Ux); Uy = real(Uy); Uz = real(Uz);

        A(p+1:3:p+3*gps,i) = Ux(1:gps);
        A(p+2:3:p+3*gps,i) = Uy(1:gps);
        A(p+3:3:p+3*gps,i) = Uz(1:gps);

        U = (Uy(gps+1:end)*sind(track) - Ux(gps+1:end)*cosd(track))...
                    *sind(look) + Uz(gps+1:end)*cosd(look);
        A(p+3*gps+1:end,i) = U;
    end

    % Take SVD of Model Perturbations
    N1 = ones(N*beta,N*beta) * (1/(N*beta));
    A_ = A* N1; Ap = A - A_;

    [U,Sig,~] = svd(Ap,'econ');
    Sig = Sig(1:N,1:N); U = U(:,1:N);

    % Generate random orthogonal matrix
    Y = randn(N,N); [~,~,V] = svd(Y,'econ');
    Ap = U * (1/sqrt(beta)) * Sig * V';
    
    % Build A Matrix
    A = zeros(6,N);
    for i = 1:p
        Ap(i,:) = Ap(i,:) - mean(Ap(i,:));
        Ap(i,:) = (devs(i)/std(Ap(i,:))) * Ap(i,:);
        A(i,:) = Ap(i,:) + A_(i,1);
    end
end

Alog(:,:,Acount) = A(1:p,:);
Astep(:,Acount) = [0;0;0;0];
Amean(:,Acount) = mean(A(1:p,:),2);
Astd(:,Acount) = std(A(1:p,:),0,2);
Asynths(:,Acount) = synths(:,1);
Acount = Acount + 1;

%% Prepare Data Matricies for Extracting Results

params = zeros(p+drv,N,reps+1);
iterations_total = zeros(1,reps+1);
count = 2; gpscount = 1; insarcount = 1;

ERR_gps = zeros(4,length(gps_steps)*gps_reps + 1);
t_gps = zeros(1,length(gps_steps)*gps_reps + 1);
model_gps = cell(gps,7);
runtime_gps = zeros(2,length(gps_steps)*gps_reps);

ERR_insar = zeros(4,length(insar_steps)*insar_reps + 1);
t_insar = zeros(1,length(insar_steps)*insar_reps + 1);
model_insar = cell(length(insar_steps)*insar_reps + 1,8);
runtime_insar = zeros(2,length(insar_steps)*insar_reps);

% Tensile and Mohr-Coulomb Failure for FEM
if model == 2
    TF = cell(N,reps+1);
    MC = cell(N,reps+1);
end

%% EnKF Setup

% Set initial values
params(1:p,:,1) = A(1:p,:);
iterations_total(1) = steps(1);

% Initialize redistributions counters
redist_count = 1;                  % Tracks how many times the model widens
iterations_redist = zeros(1,reps); % Lists time steps of widenings
params_redist = zeros(p,N,reps);   % Lists parameter sets at widening
count_narrow = zeros(p,reps);      % Number of widenings triggered by each parameter
OLcount = zeros(1,reps);           % Counts number of outliers filtered

% Initialize Loop Timer
year = 0;
tsteps = zeros(1,length(steps));
cpusteps = zeros(1,length(steps));
timer = 0;
tic
cpu0 = cputime;

spinup = 1;

%% EnKF Loop
for t = steps'
    %% Extract Measurement Vectors
    timer = timer + 1;
    
    % GPS
    xy_gps = [];  d_gps = [];
    e_gps =  [];  idx_gps = [];
    
    for i = 1:gps
        for j = 1:length(GPS{i,2})
            if GPS{i,2}(j) == t
                xy_gps = [xy_gps ; GPS{i,1}];
                d_gps = [d_gps ; GPS{i,3}(j) ; GPS{i,4}(j) ; GPS{i,5}(j)];
                e_gps = [e_gps ; GPS{i,6}(j) ; GPS{i,7}(j) ; GPS{i,8}(j)];
                idx_gps = [idx_gps ; i];
            end
        end
    end
    
    % InSAR
    xy_insar = [];  d_insar = [];
    e_insar = [];   idx_insar = [];
    
    for i = 1:insar
        if ((INSAR{i,1}(1) - t0)*ysec) == t
            d_insar = INSAR{i,3};
            e_insar = INSAR{i,4};
            xy_insar = INSAR{i,2};
            idx_insar = i;
            look = INSAR{i,1}(2); track = INSAR{i,1}(3);
        end
    end
    
    % Filter near-0 InSAR measurements
    if ~isempty(d_insar) && isfilttag
        isfilt = d_insar > 1.5 * e_insar;
        d_insar = d_insar(isfilt);
        e_insar = e_insar(isfilt);
        xy_insar = xy_insar(isfilt,:);
    end
    
    if isempty(d_insar) && isempty(d_gps)
        continue
    end
    
    %% Detect when flat ends and match step to synthetic values
    if iniskip
        if ~isempty(d_gps)
            if mean(abs(d_gps)./e_gps) > flattol
                spinup = 0;
            end
        end

        if ~isempty(d_insar)
            if mean(abs(d_insar)./e_insar) > flattol
                spinup = 0;
            end
        end

        if spinup
            continue
        end
    end
    [~,synthst] = min(abs(steps_s-t));
    synthst = synths(:,synthst);
    %% Assimilate GPS Data
    if ~isempty(d_gps)
        % Track data for matching to corresponding A matricies
        datcount = datcount + 1; 
        datlog{1,datcount} = d_gps; datlog{2,datcount} = e_gps;
        datlog{3,datcount} = xy_gps;
        
        for q = 1:gps_reps
            %% Predictive Model
            
            % Convert to standard parameters
            A = [A(1:p,:) ; zeros(3*size(xy_gps,1),N)];
            [Ap,B] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);
            
            if model == 1
                for i = 1:N
                    [Ux,Uy,Uz] = yangdisp(Ap(4,i),Ap(5,i),Ap(6,i),...
                        Ap(1,i),Ap(2,i),lambda,G,nu,Ap(3,i),theta,0,...
                        xy_gps(:,1),xy_gps(:,2),xy_gps(:,1)*0);
                    Ux = real(Ux); Uy = real(Uy); Uz = real(Uz);

                    A(p+1:3:end,i) = Ux;
                    A(p+2:3:end,i) = Uy;
                    A(p+3:3:end,i) = Uz;
                end
            else
                [Ux,Uy,Uz,MCi,TFi] = Elastic2D_17Jun19(Ap,xy_gps,E,nu);
                
                A(p+1:3:end,:) = Ux;
                A(p+2:3:end,:) = Uy;
                A(p+3:3:end,:) = Uz;
            end
            
            %% Extract predicted deformation for comparison
            for j = 1:length(idx_gps)
                i = idx_gps(j);
                
                model_gps{i,1} = [model_gps{i,1} t];
                model_gps{i,2} = [model_gps{i,2} mean(A(p+3*j-2,:))];
                model_gps{i,3} = [model_gps{i,3} mean(A(p+3*j-1,:))];
                model_gps{i,4} = [model_gps{i,4} mean(A(p+3*j,:))];
                model_gps{i,5} = [model_gps{i,5} std(A(p+3*j-2,:))];
                model_gps{i,6} = [model_gps{i,6} std(A(p+3*j-1,:))];
                model_gps{i,7} = [model_gps{i,7} std(A(p+3*j,:))];
            end
            
            %% Calculate Prediction Error
            
            ERR = zeros(1,N);
            
            for i = 1:N
                ERR(i)= sqrt(sum((A(p+1:end,i)-d_gps).^2)/length(d_gps));
            end
            
            lo = prctile(ERR,10); hi = prctile(ERR,90); emean = mean(ERR);
            
            ERR_gps(:,gpscount) = [emean ; std(ERR) ; emean-lo ; hi-emean]; 
            t_gps(gpscount) = t; 
            
            %% Retroactively Record Derived Terms
            % These are calculated as part of the current time step but are
            % associated with the old set of parameters pre-EnKF update
            
            params(p+1:end,:,count-1) = [B ; Ap([1,2,3,6],:) ; ERR];
            
            if model == 2
                MC(:,count-1) = MCi;
                TF(:,count-1) = TFi;
            end
            
            %% Apply EnKF Algorithm
            
            [A,runtime] = EnKF_14Jun21(A,d_gps,e_gps,assim,Dload);
            
            Alog(:,:,Acount) = A(1:p,:);
            Astep(:,Acount) = [1 ; datcount ; 1 ; t];
            Amean(:,Acount) = mean(A(1:p,:),2);
            Astd(:,Acount) = std(A(1:p,:),0,2);
            Asynths(:,Acount) = synthst;
            Acount = Acount + 1;
            
            runtime_gps(:,gpscount) = runtime; gpscount = gpscount + 1;
            
            %% Remove outliers with a quartile filter and resample
            
            A = A(1:p,:);
            A = qfilt_rows(A,1:p,qtol);
            
            % Determine number of ensembles to resample
            [~,j] = size(A); j = N - j;
            OLcount(count) = j;
            
            if j > 0
                ACov = cov(A'); [U,S,V] = svd(ACov); ACov = U*sqrt(S)*V';
                new_ens = randn(j,p) * ACov; new_ens = new_ens';
                for i = 1:p
                    new_ens(i,:) = new_ens(i,:) + mean(A(i,:));
                end
                A = [A new_ens];
                
                Alog(:,:,Acount) = A(1:p,:);
                Astep(:,Acount) = [2 ; datcount ; 1 ; t];
                Amean(:,Acount) = mean(A(1:p,:),2);
                Astd(:,Acount) = std(A(1:p,:),0,2);
                Asynths(:,Acount) = synthst;
                Acount = Acount + 1;
            end
            
            %% Widen ensemble if base parameters are too narrow
            
            [Ap,~] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);
            
            % Set threshold standard deviations
            widen_std = zeros(p,3);
            widen_std(:,1) = abs(mean(A(1:p,:),2)) * wth;
            
            % Widen X and Y based on depth
            widen_std(5,1) = abs(mean(Ap(6,:))) * wth;
            widen_std(6,1) = abs(mean(Ap(6,:))) * wth;
            
            for i = 1:p
                widen_std(i,1) = max([widen_std(i,1),pwiden(i)]);
            end
            
            widen_std(:,2) = widen_std(:,1) * (wset/wth);
            widen_std(:,3) = widen_std(:,1) * (wmax/wth);
            
            red_cent = mean(A(1:p,:),2);
            redist = 0; c = ones(1,p);
            
            if infl
                for i = 1:p
                    if abs(std(A(i,:))) < widen_std(i,1)
                        if ~redist
                            redist = 1; % Flags step for widening
                            params_redist(:,:,redist_count) = A(1:p,:);
                            iterations_redist(redist_count) = t;
                            redist_count = redist_count + 1;
                        end
                        c(i) = widen_std(i,2) / abs(std(A(i,:)));
                        count_narrow(i,count) = c(i);
                    end
                end
                cmax = max(c);
                
                % Widen
                for i = 1:p
                    if infl == 1 && redist % Flat redistribution
                        A(i,:) = red_cent(i) + widen_std(i,2)*(rand(1,N)-0.5);
                        
                    elseif infl == 2 && redist % Widen individually
                        A(i,:) = red_cent(i) + c(i)*(A(i,:)-mean(A(i,:)));
                        
                    elseif infl == 3 && redist % Widen most
                        if abs(std(A(i,:))) < widen_std(i,3)
                            A(i,:) = red_cent(i) + cmax*(A(i,:)-mean(A(i,:)));
                        end
                        
                    elseif infl == 4 && redist % Widen all
                        A(i,:) = red_cent(i) + cmax*(A(i,:)-mean(A(i,:)));
                        
                    elseif infl == 5 % Constant inflation factor
                        A(i,:) = red_cent(i) + kinfl*(A(i,:)-mean(A(i,:)));

                    end
                end

            end

            if redist
                Alog(:,:,Acount) = A(1:p,:);
                Astep(:,Acount) = [3 ; datcount ; 1 ; t];
                Amean(:,Acount) = mean(A(1:p,:),2);
                Astd(:,Acount) = std(A(1:p,:),0,2);
                Asynths(:,Acount) = synthst;
                Acount = Acount + 1;
            end
            
            %% Record Values
            
            params(1:p,:,count) = A(1:p,:);
            iterations_total(count) = t; count = count+1;
            
        end
        
        %% Run final model prediction if last step of GPS dataset
        if t == gps_steps(end)
            %% Predictive Model
            
            % Convert to standard parameters
            A = [A(1:p,:) ; zeros(3*size(xy_gps,1),N)];
            [Ap,B] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);
            
            if model == 1
                for i = 1:N
                    [Ux,Uy,Uz] = yangdisp(Ap(4,i),Ap(5,i),Ap(6,i),...
                        Ap(1,i),Ap(2,i),lambda,G,nu,Ap(3,i),theta,0,...
                        xy_gps(:,1),xy_gps(:,2),xy_gps(:,1)*0);
                    Ux = real(Ux); Uy = real(Uy); Uz = real(Uz);

                    A(p+1:3:end,i) = Ux;
                    A(p+2:3:end,i) = Uy;
                    A(p+3:3:end,i) = Uz;
                end
            else
                [Ux,Uy,Uz,MCi,TFi] = Elastic2D_17Jun19(Ap,xy_gps,E,nu);
                
                A(p+1:3:end,:) = Ux;
                A(p+2:3:end,:) = Uy;
                A(p+3:3:end,:) = Uz;
            end
            
            
            %% Extract predicted deformation for comparison
            
            for j = 1:length(idx_gps)
                i = idx_gps(j);
                
                model_gps{i,1} = [model_gps{i,1} t];
                model_gps{i,2} = [model_gps{i,2} mean(A(p+3*j-2,:))];
                model_gps{i,3} = [model_gps{i,3} mean(A(p+3*j-1,:))];
                model_gps{i,4} = [model_gps{i,4} mean(A(p+3*j,:))];
                model_gps{i,5} = [model_gps{i,5} std(A(p+3*j-2,:))];
                model_gps{i,6} = [model_gps{i,6} std(A(p+3*j-1,:))];
                model_gps{i,7} = [model_gps{i,7} std(A(p+3*j,:))];
            end
            
            %% Calculate Prediction Error
            
            ERR = zeros(1,N);
            
            for i = 1:N
                ERR(i)= sqrt(sum((A(p+1:end,i)-d_gps).^2)/length(d_gps));
            end
            
            lo = prctile(ERR,10); hi = prctile(ERR,90); emean = mean(ERR);
            
            ERR_gps(:,gpscount) = [emean ; std(ERR) ; emean-lo ; hi-emean]; 
            t_gps(gpscount) = t; gpscount = gpscount+1;
            
            %% Retroactively Record Derived Terms
            % These are calculated as part of the current time step but are
            % associated with the old set of parameters pre-EnKF update
            
            if t == steps(end) && isempty(d_insar)
                params(p+1:end,:,count-1) = [B ; Ap([1,2,3,6],:) ; ERR];

                if model == 2
                    MC(:,count-1) = MCi;
                    TF(:,count-1) = TFi;
                end
            end
        end
    end
    
    %% Assimilate InSAR Data
    if ~isempty(d_insar)
        % Track data for matching to corresponding A matricies
        datcount = datcount + 1; 
        datlog{1,datcount} = d_insar; datlog{2,datcount} = e_insar;
        datlog{3,datcount} = xy_insar;
        
        for q = 1:insar_reps
            %% Predictive Model
            
            % Convert to standard parameters
            A = [A(1:p,:) ; zeros(size(xy_insar,1),N)];
            [Ap,B] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);
            
            if model == 1
                for i = 1:N
                    [Ux,Uy,Uz] = yangdisp(Ap(4,i),Ap(5,i),Ap(6,i),...
                        Ap(1,i),Ap(2,i),lambda,G,nu,Ap(3,i),theta,0,...
                        xy_insar(:,1),xy_insar(:,2),xy_insar(:,1)*0);
                    Ux = real(Ux); Uy = real(Uy); Uz = real(Uz);

                    A(p+1:end,i) = (Uy*sind(track) - Ux*cosd(track))...
                        *sind(look) + Uz*cosd(look);
                end
            else
                [Ux,Uy,Uz,MCi,TFi] = Elastic2D_17Jun19(Ap,xy_insar,E,nu);
                
                A(p+1:end,:) = (Uy*sind(track) - Ux*cosd(track))...
                    *sind(look) + Uz*cosd(look);
            end
            
            %% Extract predicted deformation for comparison
            
            model_insar{insarcount,1} = INSAR{idx_insar,1};
            model_insar{insarcount,2} = xy_insar;
            model_insar{insarcount,3} = mean(A(p+1:end,:),2);
            model_insar{insarcount,4} = std(A(p+1:end,:),0,2);
            
            %% Calculate Prediction Error
            
            ERR = zeros(1,N);
            
            for i = 1:N
                ERR(i)= sqrt(sum((A(p+1:end,i)-d_insar).^2)/length(d_insar));
            end
            
            lo = prctile(ERR,10); hi = prctile(ERR,90); emean = mean(ERR);
            
            ERR_insar(:,insarcount) = [emean ; std(ERR) ; emean-lo ; hi-emean]; 
            t_insar(insarcount) = t;
            
            %% Retroactively Record Derived Terms
            % These are calculated as part of the current time step but are
            % associated with the old set of parameters pre-EnKF update
            
            params(p+1:end,:,count-1) = [B ; Ap([1,2,3,6],:) ; ERR];
            
            if model == 2
                MC(:,count-1) = MCi;
                TF(:,count-1) = TFi;
            end
            
            %% Apply EnKF Algorithm
            
            [A,runtime] = EnKF_14Jun21(A,d_insar,e_insar,assim,Dload);
            
            Alog(:,:,Acount) = A(1:p,:);
            Astep(:,Acount) = [1 ; datcount ; 2 ; t];
            Amean(:,Acount) = mean(A(1:p,:),2);
            Astd(:,Acount) = std(A(1:p,:),0,2);
            Asynths(:,Acount) = synthst;
            Acount = Acount + 1;
            
            runtime_insar(:,insarcount) = runtime; 
            insarcount = insarcount + 1;
            
            %% Remove outliers with a quartile filter and resample
            
            A = A(1:p,:);
            A = qfilt_rows(A,1:p,qtol);
            
            % Determine number of ensembles to resample
            [~,j] = size(A); j = N - j;
            OLcount(count) = j;
            
            if j > 0
                ACov = cov(A'); [U,S,V] = svd(ACov); ACov = U*sqrt(S)*V';
                new_ens = randn(j,p) * ACov; new_ens = new_ens';
                for i = 1:p
                    new_ens(i,:) = new_ens(i,:) + mean(A(i,:));
                end
                A = [A new_ens];
                
                Alog(:,:,Acount) = A(1:p,:);
                Astep(:,Acount) = [2 ; datcount ; 2 ; t];
                Amean(:,Acount) = mean(A(1:p,:),2);
                Astd(:,Acount) = std(A(1:p,:),0,2);
                Asynths(:,Acount) = synthst;
                Acount = Acount + 1;
            end
            
            %% Widen ensemble if base parameters are too narrow
            
            [Ap,~] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);
            
            % Set threshold standard deviations
            widen_std = zeros(p,3);
            widen_std(:,1) = abs(mean(A(1:p,:),2)) * wth;
            
            % Widen X and Y based on depth
            widen_std(5,1) = abs(mean(Ap(6,:))) * wth;
            widen_std(6,1) = abs(mean(Ap(6,:))) * wth;
            
            for i = 1:p
                widen_std(i,1) = max([widen_std(i,1),pwiden(i)]);
            end
            
            widen_std(:,2) = widen_std(:,1) * (wset/wth);
            widen_std(:,3) = widen_std(:,1) * (wmax/wth);
            
            red_cent = mean(A(1:p,:),2);
            redist = 0; c = ones(1,p);
            
            if infl
                for i = 1:p
                    if abs(std(A(i,:))) < widen_std(i,1)
                        if ~redist
                            redist = 1; % Flags step for widening
                            params_redist(:,:,redist_count) = A(1:p,:);
                            iterations_redist(redist_count) = t;
                            redist_count = redist_count + 1;
                        end
                        c(i) = widen_std(i,2) / abs(std(A(i,:)));
                        count_narrow(i,count) = c(i);
                    end
                end
                cmax = max(c);
                
                % Widen
                for i = 1:p
                    if infl == 1 && redist % Flat redistribution
                        A(i,:) = red_cent(i) + widen_std(i,2)*(rand(1,N)-0.5);
                        
                    elseif infl == 2 && redist % Widen individually
                        A(i,:) = red_cent(i) + c(i)*(A(i,:)-mean(A(i,:)));
                        
                    elseif infl == 3 && redist % Widen most
                        if abs(std(A(i,:))) < widen_std(i,3)
                            A(i,:) = red_cent(i) + cmax*(A(i,:)-mean(A(i,:)));
                        end
                        
                    elseif infl == 4 && redist % Widen all
                        A(i,:) = red_cent(i) + cmax*(A(i,:)-mean(A(i,:)));
                        
                    elseif infl == 5 % Constant inflation factor
                        A(i,:) = red_cent(i) + kinfl*(A(i,:)-mean(A(i,:)));
                        
                    end
                end

            end

            if redist
                Alog(:,:,Acount) = A(1:p,:);
                Astep(:,Acount) = [3 ; datcount ; 2 ; t];
                Amean(:,Acount) = mean(A(1:p,:),2);
                Astd(:,Acount) = std(A(1:p,:),0,2);
                Asynths(:,Acount) = synthst;
                Acount = Acount + 1;
            end

            %% Record Values
            
            params(1:p,:,count) = A(1:p,:);
            iterations_total(count) = t; count = count+1;
            
        end
        
        %% Run final model prediction if last step of InSAR dataset
        if t == insar_steps(end)
            %% Predictive Model
            
            % Convert to standard parameters
            A = [A(1:p,:) ; zeros(size(xy_insar,1),N)];
            [Ap,B] = EnKFParamShift(A,eps,sig,Vr,k,o,K,1);
            
            if model == 1
                for i = 1:N
                    [Ux,Uy,Uz] = yangdisp(Ap(4,i),Ap(5,i),Ap(6,i),...
                        Ap(1,i),Ap(2,i),lambda,G,nu,Ap(3,i),theta,0,...
                        xy_insar(:,1),xy_insar(:,2),xy_insar(:,1)*0);
                    Ux = real(Ux); Uy = real(Uy); Uz = real(Uz);

                    A(p+1:end,i) = (Uy*sind(track) - Ux*cosd(track))...
                        *sind(look) + Uz*cosd(look);
                end
            else
                [Ux,Uy,Uz,MCi,TFi] = Elastic2D_17Jun19(Ap,xy_insar,E,nu);
                    
                A(p+1:end,:) = (Uy*sind(track) - Ux*cosd(track))...
                    *sind(look) + Uz*cosd(look);
            end
            
            %% Extract predicted deformation for comparison
            
            model_insar{insarcount,1} = INSAR{idx_insar,1};
            model_insar{insarcount,2} = xy_insar;

            model_insar{insarcount,3} = mean(A(p+1:end,:),2);
            model_insar{insarcount,4} = std(A(p+1:end,:),0,2);
            
            %% Calculate Prediction Error
            
            ERR = zeros(1,N);
            
            for i = 1:N
                ERR(i)= sqrt(sum((A(p+1:end,i)-d_insar).^2)/length(d_insar));
            end
            
            lo = prctile(ERR,10); hi = prctile(ERR,90); emean = mean(ERR);
            
            ERR_insar(:,insarcount) = [emean ; std(ERR) ; emean-lo ; hi-emean]; 
            t_insar(insarcount) = t; insarcount = insarcount + 1;
            
            %% Retroactively Record Derived Terms

            if t == steps(end)
                params(p+1:end,:,count-1) = [B ; Ap([1,2,3,6],:) ; ERR];

                if model == 2
                    MC(:,count-1) = MCi;
                    TF(:,count-1) = TFi;
                end
            end
        end
    end
    
    %% Update timer and report progress
    
    tsteps(timer) = toc;
    cpusteps(timer) = cputime - cpu0;
    
    if (t/ysec) >= year + 0.1
        year = year + 0.1;
    end
    
    clc
    fprintf(['Completed Time Step: ' num2str(timer-1) '/'...
        num2str(length(steps)) '\n' 'Year: ' num2str(year) '/'...
        num2str(tend) '\n'])
    
end

%% Cleanup
Alog(:,:,Acount:end) = [];
Astep(:,Acount:end) = [];
Amean(:,Acount:end) = [];
Astd(:,Acount:end) = [];
Asynths(:,Acount:end) = [];
datlog(:,datcount+1:end) = [];

params(:,:,count:end) = [];
iterations_total(count:end) = [];
if model == 2
    MC(:,count:end) = [];
    TF(:,count:end) = [];
end
count_narrow(:,count:end) = [];
OLcount(count:end) = [];

ERR_gps(:,gpscount:end) = []; 
t_gps(gpscount:end) = [];
runtime_gps(:,gpscount:end) = [];

ERR_insar(:,insarcount:end) = []; 
t_insar(insarcount:end) = [];
runtime_insar(:,insarcount:end) = [];
model_insar(insarcount:end,:) = [];

params_redist(:,:,redist_count:end) = [];
iterations_redist(redist_count:end) = [];

%% Calculate means at each time step and median values across assimilation
% Calculate means at each time step and median value across assimilation
means = zeros(p+drv,count-1); std_devs = zeros(p+drv,count-1);
Pmeans = mean(params,2); Pstd = std(params,0,2);
for i = 1:count-1
    means(:,i) = Pmeans(:,:,i);
    std_devs(:,i) = Pstd(:,:,i);
end

medians = zeros(1,p+drv);
for i = 1:p+drv
    medians(i) = median(means(i,:));
end

%% Convert units of output times to years
itr_years = (iterations_total / ysec);
itr_redistyr = (iterations_redist / ysec);

t_gpsyr = (t_gps / ysec);
modelgps_yr = cell(gps,1);
gps_yr = cell(gps,1);
for i = 1:gps
    modelgps_yr{i,1} = (model_gps{i,1}/ysec);
    gps_yr{i,1} = (GPS{i,2}/ysec);
end

t_insaryr = (t_insar/ysec);

%% Save Outputs

root = cd(saveloc);

save([savename '_' datestr(now,'ddmmmyy') '.mat'])

cd(root)

end

