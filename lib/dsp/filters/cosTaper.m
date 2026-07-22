    function m = cosTaper(a, f1, f2, tw)
        % Mask = 1 inside [f1,f2], 0 outside [f1-tw, f2+tw], with raised-cosine
        % transitions in between. tw <= 0 degenerates to a brick-wall.
        m = zeros(size(a));
        m(a >= f1 & a <= f2) = 1;
        if tw > 0
            lo = (a >= f1-tw) & (a < f1);
            m(lo) = 0.5 * (1 + cos(pi * (f1 - a(lo)) / tw));
            hi = (a > f2) & (a <= f2+tw);
            m(hi) = 0.5 * (1 + cos(pi * (a(hi) - f2) / tw));
        end
    end
