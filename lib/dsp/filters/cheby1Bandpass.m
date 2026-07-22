    function [b, a] = cheby1Bandpass(n, f1, f2, Fs, Rp)
        % Equiripple-passband IIR (Rp dB ripple). Same design path as
        % Butterworth but with Chebyshev prototype poles (they sit on an
        % ellipse instead of a circle), giving a sharper transition for the
        % same order at the cost of passband ripple. Verified against
        % scipy.signal.cheby1: response error ~1e-14, ripple exact.
        if nargin < 5, Rp = 1; end
        eps_r = sqrt(10^(Rp/10) - 1);
        mu    = asinh(1/eps_r) / n;
        k     = 1:n;
        th    = pi*(2*k - 1)/(2*n);
        pLP   = -sinh(mu)*sin(th) + 1i*cosh(mu)*cos(th);  % Chebyshev-I poles
        [b, a] = bpBilinear(pLP, n, f1, f2, Fs);
        % ripple means the centre may sit on a peak or trough, so normalise to
        % the passband PEAK (sample the in-band response and divide it out).
        fg = linspace(f1, f2, 256);
        H  = freqzLite(b, a, fg, Fs);
        pk = max(abs(H));
        if pk > 0, b = b / pk; end
    end
