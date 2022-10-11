function [ B ] = qfilt_rows( A, rows, mult )
% Removes colums from matrix A that contain outliers in row 'row' defined as
%  points farther than 1.5 interquartile ranges from the upper 
%  and lower quartiles

%% Changelog

% 10 Sep 19 - Function returning empty set when given 1-row matricies,
%   adjusting to return original matrix if containing less then 5 points

% 26 Feb 21 - Changing from columns to rows

%% Function
if size(A,1) < 5
    B = A;
else
    IN = ones(1,size(A,2));
    for row = rows
        
        Arow = A(row,:);

        M = median(Arow);
        upper = Arow(Arow>=M); lower = Arow(Arow<=M);

        one = median(lower); three = median(upper);
        IQR = (three - one)*mult;

        in = Arow<(three+IQR) & Arow>(one-IQR);
        IN = IN & in;
    end
    
    B = A(:,IN);
end

end

