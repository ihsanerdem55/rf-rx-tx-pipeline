    function y = fftMask(x, Fs, f1, f2, tw)
        % Zero-phase band-pass by multiplying the FFT with a mask.
        % tw = 0 -> hard brick-wall; tw > 0 -> raised-cosine edges.
        % A hard spectral edge is a sinc in time, which is exactly what makes
        % the pulse corners ring; the taper trades selectivity for less ringing.
        x = x(:).'; N = numel(x);
        X = fft(x);
        f = (0:N-1) / N * Fs;
        fbi = f; fbi(f > Fs/2) = f(f > Fs/2) - Fs;   % wrap to a two-sided axis
        mask = cosTaper(abs(fbi), f1, f2, tw);
        y = real(ifft(X .* mask));                   % real() kills round-off imag
    end
