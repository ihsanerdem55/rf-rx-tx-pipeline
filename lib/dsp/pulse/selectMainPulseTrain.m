    function [rise,fall] = selectMainPulseTrain(riseAll,fallAll,env,P)
        % Multipath often creates two detections per PRI. Rank runs by peak
        % envelope and greedily retain at most one strong run per minimum-PRI
        % neighbourhood. This remains an estimator (no truth timing is used).
        rise = []; fall = [];
        if isempty(riseAll), return; end
        amp = zeros(size(riseAll));
        for k = 1:numel(riseAll)
            amp(k) = max(env(riseAll(k):fallAll(k)));
        end
        [~,ord] = sort(amp,'descend');
        minSep = max(1,round(0.45*P.PRI*(1-P.pridev)*P.Fs));
        keep = [];
        for ii = 1:numel(ord)
            q = ord(ii);
            if isempty(keep) || all(abs(riseAll(q)-riseAll(keep))>=minSep)
                keep(end+1) = q; %#ok<AGROW>
                if numel(keep)>=P.Np, break; end
            end
        end
        keep = sort(keep);
        rise = riseAll(keep); fall = fallAll(keep);
    end
