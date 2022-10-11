function [ Ux, Uy, Uz, MC, TF ] = Elastic2D_17Jun19( A, xy, E, nu )
% Given a set of parameters A, a list of observation coordinates xy,
% material properties E and eta, and cohesions coh, outputs 3-component
% surface deformation as well as measures of tensile and Mohr-Coulomb
% failure

%% Changelog

% 20 May 2019 - Original version based on FEM_Ellip_30Jul18.m

% 04 June - Exceeded memory capacity, attempting version that only extracts
%   TF and MC when factors exceed 0

% 05 June - Memory problem persisted, reducing radius of MC search cylinder
%   to 2x lateral extent, altered both here and in original COMSOL file

% 06 June - Continued memory problems, extracting MC failure as
%   cross-sections along and pependicular to plunge axis

% 10 June - Iterating outside of function avoids crash but multiple model
%   calls seem to cause progressive slowdown. Reverting June 6 changes and
%   altering function to perform multiple runs on same model.

% 17 Jun - Altering for use with 2.5D axisymmetric model

%% Setup

% Extract model parameters
N = size(A,2);
points = size(xy,1);

r1    = A(1,:); % Half height
r2    = A(2,:); % Half width
dP    = A(3,:); % Pressure
dX    = A(4,:); % Lateral Offsets
dY    = A(5,:);
dZ    = A(6,:); % Depth to center


% Setup Outputs
MC = cell(N,1); TF = cell(N,1);
Ux = zeros(points,N); Uy = zeros(points,N); Uz = zeros(points,N);

%% Prepare COMSOL model

model = mphload('2DElasticEllip_17Jun19.mph');
model.param.set('E',E);
model.param.set('nu',nu);

%% Loop over ensemble members to extract data

for i = 1:N
    
    % Set parameters and run model
    model.param.set('r1',r1(i));       model.param.set('r2',r2(i));
    model.param.set('dP',dP(i));       model.param.set('dZ',dZ(i));
    model.study('std1').run
    
    % Extract offsets from model
    R = sqrt((xy(:,1)-dX(i)).^2 + (xy(:,2)-dY(i)).^2);
    theta = atan2(xy(:,2)-dY(i),xy(:,1)-dX(i));
    Z = zeros(points,1);
    p = [R';Z'];
    
    [Ur,Uzi] = mphinterp(model,{'u','w'},'coord',p);
    Ux(:,i) = Ur' .* cos(theta); Uy(:,i) = Ur' .* sin(theta); Uz(:,i) = Uzi;
    
    % Extract location and magnitude of maximum tensile stress on wall 
    pd = mpheval(model,{'solid.sp1','alpha'},'edim',1,'selection',[11 12]);
    TF{i,1} = [pd.d1 ; pd.p ; pd.d2];
    
    % Extract MC failure in surrounding rock
    pd = mpheval(model,{'MCFail','meshvol'},'edim','domain',...
        'selection',2);
    id = pd.d1 > 0;
    MC{i,1} = [pd.d1(id) ; pd.d2(id) ; pd.p(:,id)];
    
    if i > 1
        fprintf(repmat('\b',1,strlen))
    end
    strlen = fprintf(['Extracted ensemble member ' num2str(i) '/' num2str(N) '\n']);
    
end

end
