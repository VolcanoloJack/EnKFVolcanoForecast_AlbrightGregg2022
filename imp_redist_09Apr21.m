function [ A2, ext, breach, far ] = imp_redist_09Apr21( A, lims )
% Redistributes physically impossible/infeasible values as part of EnKF
% workflow. Outputs redistributed parameter matrix A and counters for
% number of redistributions implemented

%% Changelog

% 05 Sep 19 - Filters now check against current parameter values instead of
%   original values, prevents early filters creating impossible values that
%   aren't caught by later filters

% 08 Mar 21 - Rewrite to improve functionality in more extreme cases
%   - Never resets to a single value like mean; if entire ensemble gets
%      near an extreme a single value would cause collapse
%   - For extreme cases, reflects across the boundary. The more anomalous
%      the value, the stronger the correction
%   - Not using the parameter scaling term because the only parameter
%      changed in typical workflow is Pressure which doesn't appear here

% 12 Mar 21 - Option for square root scaling

% 09 Apr 21 - Changed to adjustable root scaling

%% Setup

rt = 0.9; % Exponent of offset
brdiv = 0.75; % Amount of adjustment partioned to depth increase for breach

N = size(A,2);
ext = zeros(3,2); breach = 0; far = 0;
loops = 0; cont = 0;

A2 = A;

while cont == 0
    ext0 = ext; breach0 = breach; far0 = far; loops = loops + 1;
    
    A3 = A2;
    
    for i = 1:N
        
        %% Half Height
        
        if A2(1,i) < lims(1,1) % Height too small
            dist = lims(1,1) - A2(1,i);
            A2(1,i) = lims(1,1) + (dist^rt);
            ext(1,1) = ext(1,1) + 1;
        elseif A2(1,i) > lims(1,2) % Height too large
            dist = A2(1,i) - lims(1,2);
            A2(1,i) = lims(1,2) - (dist^rt);
            ext(1,2) = ext(1,2) + 1;
        end
        
        %% Half Width
        
        if A2(2,i) < lims(2,1) % Width too small
            dist = lims(2,1) - A2(2,i);
            A2(2,i) = lims(2,1) + (dist^rt);
            ext(2,1) = ext(2,1) + 1;
        elseif A2(2,i) > lims(2,2) % Width too large
            dist = A2(2,i) - lims(2,2);
            A2(2,i) = lims(2,2) - (dist^rt);
            ext(2,2) = ext(2,2) + 1;
        end
        
        %% Depth
        
        if A2(6,i) < lims(3,1) % Too Shallow
            dist = lims(3,1) - A2(6,i);
            A2(6,i) = lims(3,1) + (dist^rt);
            ext(3,1) = ext(3,1) + 1;
        elseif A2(6,i) > lims(3,2) % Too Deep
            dist = A2(6,i) - lims(3,2);
            A2(6,i) = lims(3,2) - (dist^rt);
            ext(3,2) = ext(3,2) + 1;
        end
        
        %% Breach at Surface
        if A2(1,i) >= A2(6,i) - lims(4,1)
            dist = lims(4,1) + A2(1,i) - A2(6,i);
            offset = dist + (dist^rt);
            A2(1,i) = A2(1,i) - (offset * (1-brdiv));
            A2(6,i) = A2(6,i) + (offset * brdiv);
            
            breach = breach + 1;
        end
        
        %% Too Distant
        if sqrt(A2(4,i)^2 + A2(5,i)^2) > lims(4,2)
            dist = sqrt(A2(4,i)^2 + A2(5,i)^2) - lims(4,2);
            rad = lims(4,2) - (dist^rt); 
            ratio = rad/sqrt(A2(4,i)^2 + A2(5,i)^2);
            A2(4,i) = A2(4,i) * ratio; A2(5,i) = A2(5,i) * ratio; 
            far = far + 1;
        end
    end
        
        %% Check for End of Loop
    if isequal(ext0,ext) && isequal(breach0,breach) && isequal(far0,far)
        cont = 1;
    elseif loops >= 5
        fprintf('Impossible Redist Failed \n')
        cont = 1;
    end
end

end

