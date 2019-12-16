function [faRate, hitRate, Acc] = ROCcurve(score, trueLabel)
%% Plot a receiver operating curve
% function [faRate, hitRate, Acc] = plotROCcurve(score, trueLabel)
%
% score(i) = confidence in i'th detection (bigger means more confident)
% trueLabel(i) = 0 if background or 1 if target
%
% faRate(t) = false alarm rate at t'th threshold
% hitRate(t) = detection rate at t'th threshold 
% Acc = Accuracy

class1 = find(trueLabel==1);
class0 = find(trueLabel==0);

thresh = sort(score);
Nthresh = length(thresh);
hitRate = zeros(1, Nthresh);
faRate = zeros(1, Nthresh);
rec = zeros(1, Nthresh);
for thi=1:length(thresh)
    th = thresh(thi);
    % hit rate = TP/P
    hitRate(thi) = sum(score(class1) >= th) / length(class1);
    % fa rate = FP/N
    faRate(thi) = sum(score(class0) >= th) / length(class0);
    rec(thi) = (sum(score(class1) >= th) + sum(score(class0) < th)) / Nthresh;
end

% Accuracy
Acc = max(rec);

plot(faRate,hitRate)
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression')

end