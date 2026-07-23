function s = scoreBadge(metric, goodThr, fairThr)
% SCOREBADGE  Compact ASCII pass/marginal/off badge for the accuracy scorecard.
%   metric, goodThr and fairThr must share the same units; smaller is better.
%     metric <= goodThr -> [ + ]  (within tolerance)
%     metric <= fairThr -> [ ~ ]  (marginal)
%     otherwise         -> [ ! ]  (off)
%     non-finite metric -> [ ? ]  (not available, e.g. estimator returned NaN)
    if ~isfinite(metric)
        s = '[ ? ]';
    elseif metric <= goodThr
        s = '[ + ]';
    elseif metric <= fairThr
        s = '[ ~ ]';
    else
        s = '[ ! ]';
    end
end
