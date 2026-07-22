    function m = gaussMask(a, f1, f2)
        % Smooth Gaussian band-pass shape, -3 dB exactly at f1 and f2 (so the
        % arithmetic band centre is the peak). Because the transition is
        % infinitely smooth, the time-domain response has essentially no ringing
        % -- the opposite extreme from the brick-wall.
        fc = 0.5*(f1 + f2);
        sigma = (f2 - f1) / (2*sqrt(log(2)));            % |H(f1)|=|H(f2)|=1/sqrt(2)
        m = exp(-0.5 * ((abs(a) - fc)/sigma).^2);
    end
