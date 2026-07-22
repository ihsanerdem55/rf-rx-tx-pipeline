    function n=countFreqLevels(f,tol)
        if isempty(f), n=0; return; end
        fs=sort(f); n=1; c=fs(1);
        for k=2:numel(fs)
            if abs(fs(k)-c)>tol, n=n+1; c=fs(k); else, c=0.7*c+0.3*fs(k); end
        end
    end
