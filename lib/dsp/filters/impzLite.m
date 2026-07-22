    function hh = impzLite(b, a, N)
        % Impulse response = filter a unit impulse (impz() is toolbox-only).
        imp = [1, zeros(1, N-1)];
        hh = filter(b, a, imp);
    end
