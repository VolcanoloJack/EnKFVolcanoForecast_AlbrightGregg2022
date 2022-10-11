function    X2v = X2v_spheroidVertical(VAR,a,mu,nu,x,y,Wd,data)
% chi square function for a prolate spheroid
%
% VAR   vector with unknown parameters x0, y0, z0, A, P_G, theta and phi
%
% SOURCE PARAMETERS
% a         semimajor axis
% A         geometric aspect ratio [dimensionless]
% P_G       dimennsionless excess pressure (pressure/shear modulus) 
% x0,y0     surface coordinates of the center of the prolate spheroid
% z0        depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% theta     plunge (dip) angle [deg] [90 = vertical spheroid]
% phi       trend (strike) angle [deg] [0 = aligned to North]
%
% CRUST PARAMETERS
% mu        shear modulus
% nu        Poisson's ratio 
%
% BENCHMARKS
% x,y       benchmark location
% z         depth within the crust (z=0 is the free surface)
%
% x0,y0     coordinates of the center of the sphere 
% z0        depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% P_G       dimensionless excess pressure (pressure/shear modulus)
% b     radius of the sphere
% nu    Poisson's ratio
% x,y   data point location
% Wd    weight matrix
% data  deformation data vector

    x0 = VAR(1); y0 = VAR(2); z0 = VAR(3); A = VAR(4); dP = VAR(5); 
    theta = 89.99; phi = 0; 
    [u, v, w] = yang(x0,y0,z0,a,A,dP,mu,nu,theta,phi,x,y,0);
    p = length(VAR);                        % number of parameters


model = [u v w]; model = model(:);
r = data - model;                       % residual
X2 = r'*Wd*r;                           % Chi Square - full covariance
N = length(data);                       % number of data points
X2v = X2/(N-p);                         % Chi Square per degrees of freedom