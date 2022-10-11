function [ Vmin, Vmax ] = VlimsTheta( theta, eps, sig )
% Given a set of aspect ratios theta and bounds on half-width and height
% (contained in eps and sig), calculate the smallest and largest possible
% volumes

%% Changelog

% 18 Jun 21 - Original version

%%

k = (4*pi())/3;

Vmin = zeros(size(theta)); Vmax = zeros(size(theta));

for i = 1:numel(theta)
    % Calculate minimum
    if theta(i) < eps(2)/eps(1)
        Vmin(i) = k*eps(2)^3/theta(i);
    else
        Vmin(i) = k*eps(1)^3 * theta(i)^2;
    end
    
    % Calculate maximum
    if theta(i) < (eps(2)+sig(2))/(eps(1)+sig(1))
        Vmax(i) = k*(eps(1)+sig(1))^3 * theta(i)^2;
    else
        Vmax(i) = k*(eps(2)+sig(2))^3 / theta(i);
    end
end


end

