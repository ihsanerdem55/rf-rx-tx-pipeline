    function y = filtfiltLite(b, a, x)
        % Zero-phase IIR filtering: run forward, then backward. The phase shifts
        % cancel exactly; the magnitude is applied twice, so the effective
        % response is |H|^2 (the tab-2 plot accounts for this).
        % Reflection padding suppresses the start/end transients.
        x = x(:).'; N = numel(x);
        npad = min(N-1, 3*max(numel(a), numel(b)));
        if npad > 0
            pre  = 2*x(1)   - x(npad+1:-1:2);
            post = 2*x(end) - x(end-1:-1:end-npad);
            xp = [pre, x, post];
        else
            xp = x;
        end
        y = filter(b, a, xp);
        y = filter(b, a, fliplr(y));
        y = fliplr(y);
        y = y(npad+1 : npad+N);
    end
