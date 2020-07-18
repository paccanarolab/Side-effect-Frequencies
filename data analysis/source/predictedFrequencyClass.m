function [y] = predictedFrequencyClass(x)
% PREDICTEDFREQUENCYCLASS returns the frequency class given the predicted 
% score in x.
% For more details see methods in Galeano et al. Predicting the Frequency
% of Drug Side effects.
% Copyright (C) 2017 by Diego Galeano.
boundaryThresholds = [0.42, 1.26,2.43, 3.25, 3.93];
y = zeros(length(x), 1);

for i = 1:length(x)
    if x(i) < boundaryThresholds(1)
        y(i) = 0;
    elseif x(i) <= boundaryThresholds(2)
        y(i) = 1;
    elseif x(i) <= boundaryThresholds(3)
        y(i) = 2;

    elseif x(i) <= boundaryThresholds(4)
        y(i) = 3;

    elseif  x(i) <= boundaryThresholds(5)
        y(i) = 4;

    else
        y(i) = 5;
    end
end



