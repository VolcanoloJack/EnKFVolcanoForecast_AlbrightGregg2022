function [Y] = logist_18Jun21(X,k,o,inv)
% Element-wise calculation of X (all real numbers) to Y (0 to 1) according
% to a logistic function, inverts if inv == 1

% 18 Jun 21 - Adding shift factor o

Y = zeros(size(X));

for i = 1:numel(X)
    if inv
        Y(i) = (-1/k)*log((1-X(i))/X(i)) + o;
    else
        Y(i) = 1/(1+exp(-k*(X(i)-o)));
    end
end

end

