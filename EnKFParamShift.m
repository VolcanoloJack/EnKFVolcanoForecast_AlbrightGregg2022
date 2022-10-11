function [ A2, B ] = EnKFParamShift( A, eps, sig, Vr, k, o, K, inv )
% Converts back and forth between standard model parameters
%   (r1, r2, dP, dX, dY, dZ) and reformulated parameters
%   (dV, gamma, psi, delta, dX, dY), along with intermediate parameters
%   in B
%% Changelog

% 18 Jun 21 - Original Version

%% Reformulated -> Standard
if inv == 1 
    
    dV = A(1,:); gamma = A(2,:); psi = A(3,:);
    delta = A(4,:); dX = A(5,:); dY = A(6,:);
    
    % Aspect Ratio
    thmin = eps(2)/(eps(1)+sig(1));
    thmax = (eps(2)+sig(2))/eps(1);
    
    theta = thmin + (thmax-thmin)*logist_18Jun21(gamma,k(1),o(1),0);
    
    % Volume
    [Vmin,Vmax] = VlimsTheta(theta,eps,sig);
    V = Vmin + (Vmax - Vmin).*logist_18Jun21(psi,k(2),o(2),0);
    
    % Pressure
    dP = (Vr*dV*K)./V;
    
    % Reservoir Dimensions
    kV = (4*pi())/3;
    r1 = (V./(kV*theta.^2)).^(1/3); r2 = (V.*theta/kV).^(1/3);
    
    % Depth
    dZ = eps(3) + r1 + (sig(3)-r1).*logist_18Jun21(delta,k(3),o(3),0);
    
    % Collect
    A2 = [r1;r2;dP;dX;dY;dZ];
    B = [V;theta];
    
%% Standard -> Reformulated
else 
    
    r1 = A(1,:); r2 = A(2,:); dP = A(3,:);
    dX = A(4,:); dY = A(5,:); dZ = A(6,:);
    
    % Aspect Ratio, Volume, and Volume Change
    theta = r2./r1;
    V = (4*pi()/3)*r1.*r2.^2;
    dV = ((V.*dP)/K) / Vr;
    
    % Delta
    delta = (dZ - r1 - eps(3))./(sig(3)-r1);
    delta = logist_18Jun21(delta,k(3),o(3),1);
    
    % Gamma
    thmin = eps(2)/(eps(1)+sig(1));
    thmax = (eps(2)+sig(2))/eps(1);
    
    gamma = (theta - thmin)/(thmax - thmin);
    gamma = logist_18Jun21(gamma,k(1),o(1),1);
    
    % Psi
    [Vmin,Vmax] = VlimsTheta(theta,eps,sig);
    psi = (V - Vmin)./(Vmax - Vmin);
    psi = logist_18Jun21(psi,k(2),o(2),1);
    
    % Collect
    A2 = [dV;gamma;psi;delta;dX;dY];
    B = [V;theta];
    
end

end

