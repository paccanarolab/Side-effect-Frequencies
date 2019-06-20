function [AUC, XROC, YROC, TROC] = getAUROC(Res,YTest, YTrain)
%AUROC returns area under the receiver operating curve (ROC)
%  input: 
%       + Res: matrix of predicted scores
%       + YTest: matrix of test set.
%       + YTrain: matrix of training set.
%  output:
%       + AUC: area under the ROC curve.
%       + XROC: X coordinate of the ROC curve.
%       + YROC: Y coordinate of the ROC curve.
%       + TROC: thresholds of the ROC curve.

        poslabels = YTest > 0;
        neglabels = YTrain == 0 & YTest == 0;             
        
        no_poslabels = sum(poslabels(:));
        no_neglabels = sum(neglabels(:));
        fprintf('\n number of positive labels %d negative %d\n', no_poslabels, no_neglabels);
        
        % Define the labels
        labels = zeros(no_poslabels + no_neglabels,1);
        labels(1:no_poslabels) = ones(no_poslabels,1);

        % Get the predicted scores.
        scores = zeros(no_poslabels + no_neglabels,1);
        scores(1:no_poslabels) = Res(poslabels);   
        scores(no_poslabels+1:end) = Res(neglabels);   

        [XROC,YROC,TROC, AUC] = perfcurve(labels,scores,1);
        

end

