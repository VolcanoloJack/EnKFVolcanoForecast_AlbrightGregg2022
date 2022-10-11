function [A2,t] = EnKF_14Jun21(A,d,e,mode,load)
% Unified function for various EnKF implementations. Computes EnKF update
% of ensemble matrix A -> A2 given measurements d and error vector e
% derived from original data source. 'mode' determines which implementation
% is used: 1 = stochastic pseudoinverse, 2 = stochastic low rank Re,
% 3 = stochastic rank reduction, 4 = stochastic rank reduction (bugged), 
% 5 = square root pseudoinverse, 6 = square root random rotation, 7 = DEnKF

%% Changelog

% 03 Mar 21 - Original version synthesized from sEnKF_24Feb21,
%  DEnKF_15Dec20, and SqRtEnKF_01Mar21
%   - changed inversion step in DEnKF to pseudoinverse

% 07 Jun 21 - Removal of small eigenvalues was distorting inversion, so
%   rolled back

% 10 Jun 21 - Error was not in eigenvalues but in rank of C, reverting and
%   instead changing the diagonal entries of Re to match variance

% 14 Jun 21 - Adding in diagonal loading to stabilize stochastic inversion

%% Extract ensemble and measurement sizes
t = zeros(2,1);
t(2) = cputime;
tic

[p,N] = size(A); m = length(d);

%% Build necessary matricies

H = zeros(m,p); % Translation matrix
H(:,p-m+1:p) = eye(m);
HA = H*A;

N1 = ones(N,N) * (1/N);
A_ = A * N1; % Ensemble mean matrix
Ap = A - A_; % Ensemble perturbation matrix

S = H*Ap;

%% Stochastic Models
if mode <= 4
    D = zeros(m,N);
    Y = zeros(m,N);

    for i = 1:m
        Y(i,:) = e(i) * randn(1,N);
        D(i,:) = Y(i,:) + d(i);
    end
    
    Dp = D - H*A;
    Re = (1/(N-1))*(Y*Y');
    for i = 1:m
        Re(i,i) = e(i)^2;
    end
    
    if mode == 1 % Stochastic pseudoinverse
        C = (S*S')+(N-1)*Re;
        if load
            diagload = load*mean(diag(C))/10;
        else
            diagload = mean(diag((N-1)*Re))/10;
        end
        C = C + diagload*eye(m);
        [Zc,Delc] = eig(C); Delc = diag(Delc);
        Delc = Delc.^-1;
        q = length(Delc);
        Delc = diag(Delc);
        Zc = Zc(:,1:q);

        Cinv = Zc * Delc * Zc';
        X = eye(N) + S'*Cinv*Dp;
        A2 = A*X;
        
    elseif mode == 2 % Stochastic Low Rank Re
        [U0,Sig0,~] = svd(S,'econ');
        sig = diag(Sig0); sig = sig.^-1; q = length(sig);
        Sig0t = diag(sig);
        U0 = U0(:,1:q);

        X0 = Sig0t*U0'*Y;
        [U1,Sig1,~] = svd(X0,'econ');

        X1 = U0*Sig0t'*U1;
        Cinv = X1*((eye(q)+Sig1^2)^(-1))*X1';
        X = eye(N) + S'*Cinv*Dp;
        A2 = A*X;
        
    elseif mode == 3 % Stochastic Low Rank
        HA_ = HA*N1;
        HAp = HA - HA_;
        HApE = HAp + Y; 

        [U, sig, ~] = svd(HApE,'econ'); % Decompose eigenvalues
        sig = diag(sig);
        q = length(sig);
        sig = diag(sig);
        U = U(:,1:q);

        sig2 = sig*sig';
        sigt = sig2^(-1);

        X = sigt * U';
        X = X * Dp;
        X = U*X;
        X = (HAp)' * X;

        A2 = A + Ap*X;
        
    elseif mode == 4 % Stochastic Low Rank (from Okmok paper)
        HA_ = HA*N1;
        HAp = HA - HA_;
        HApE = HAp + Y; 

        [U, sig, ~] = svd(HApE); % Decompose eigenvalues
        sig2 = (sig(find(sig > 0)))';
        sigN = sig2*sig2';
        sigN2 = sigN^(-1);

        X = sigN2 * U';
        X = X * Dp;
        X = U*X;
        X = (HAp)' * X;

        A2 = A + Ap*X;
    end

%% Square Root Models
elseif mode <= 6
    D = d*ones(1,N); 
    Re = diag(e.^2);
    C = S*S' + (N-1)*Re;
%     diagload = load*mean(diag((N-1)*Re))/10;
    diagload = load*mean(diag(C))/10;
    C = C + diagload*eye(m);
    
    if mode == 6
        w = ones(N,1) / sqrt(N);
        [w1,~] = svd(w);
        [w2,~] = qr(randn(N-1,N-1));
        RR = w1 * blkdiag(1,w2) * w1';
    else
        RR = eye(N);
    end
    
    [Z,Del] = eig(C);
    Del = diag(Del); Del = Del.^-1;
    q = length(Del); Del = diag(Del); Z = Z(:,1:q);
    
    Cinv = Z * Del * Z';
    
    
    A_a = A_ + Ap*S'*Cinv*(D-H*A_);
    X2 = Del^(1/2)*Z'*S;
    
    [~,Sig,V] = svd(X2);
    Apa = Ap*V*((eye(N)-Sig'*Sig)^(1/2))*V'*RR';
    
    A2 = A_a + Apa;
    
%% DEnKF
elseif mode == 7
    n1 = ones(N,1);
    x = A_(:,1); % Mean vector
    
    % Assuming errors uncorrelated, generate observation covariance
    Rp = eye(m); % Inverse square root of Re
    for i = 1:m
        Rp(i,i) = 1/e(i);
    end
    
    % Standardize observation anomalies and innovations
    s = Rp * (d - H * x) / sqrt(N-1);
    S = Rp * (H * Ap) / sqrt(N-1);

    % Calculate update factors
    if N < m
        G = (eye(N) + S'*S);
        [U,sig,V] = svd(G,'econ');
        sig = diag(sig); sig = sig.^-1; 
        q = length (sig); sig = diag(sig);
        Ginv = U(:,1:q) * sig * V(:,1:q)';
        
        G = Ginv * S';
    else
        G = (eye(m) + S*S');
        [U,sig,V] = svd(G,'econ');
        sig = diag(sig); sig = sig.^-1; 
        q = length (sig); sig = diag(sig);
        Ginv = U(:,1:q) * sig * V(:,1:q)';
        
        G = S' * Ginv;
    end
    TR = (eye(N) - 0.5*(G*S));
    w = G * s;
    X5 = w * n1' + TR;

    % Update Ensemble
    A2 = A * X5;
end

t(1) = toc;
t(2) = cputime - t(2);
end

