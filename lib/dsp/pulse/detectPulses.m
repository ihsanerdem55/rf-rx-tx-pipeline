    function [rise, fall] = detectPulses(x, PW, Fs)
        % Reusable pulse detector: envelope -> half-max threshold -> contiguous
        % above-threshold runs. Returns rising/falling sample indices. Used both
        % for the received signal and for the TRANSMITTED signal (so DRFM false
        % targets are read off the real waveform, overlaps included).
        if isreal(x), z = analytic(x); else, z = x; end
        env = abs(z);
        w = max(1, round(0.15 * PW * Fs));
        env = movmean(env, w);
        thr = median(env) + 0.5*(max(env) - median(env));
        b = env > thr;
        d = diff([0 b 0]);
        rise = find(d ==  1);
        fall = find(d == -1) - 1;
        n = min(numel(rise), numel(fall));
        rise = rise(1:n); fall = fall(1:n);
    end
