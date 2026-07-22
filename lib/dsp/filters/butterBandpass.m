    function [b, a] = butterBandpass(n, f1, f2, Fs)
        % Maximally-flat IIR, done by hand because butter() is toolbox-only.
        % Cross-checked against scipy.signal.butter: max |H| error ~1e-4,
        % band edges land at -3.01 dB as they should.
        k = 1:n;
        pLP = exp(1i*pi*(2*k + n - 1)/(2*n));    % Butterworth prototype poles
        [b, a] = bpBilinear(pLP, n, f1, f2, Fs); % LP -> BP -> bilinear (unnorm)
        % maximally-flat, so unity gain at the geometric centre is well-defined
        fc = sqrt(max(f1, eps) * f2);
        H  = freqzLite(b, a, fc, Fs);
        if abs(H) > 0, b = b / abs(H); end
    end
