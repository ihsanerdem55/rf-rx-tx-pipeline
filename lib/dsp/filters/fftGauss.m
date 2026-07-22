    function y = fftGauss(x, Fs, f1, f2)
        % Zero-phase Gaussian band-pass applied by multiplying the FFT.
        x = x(:).'; N = numel(x);
        X = fft(x);
        f = (0:N-1) / N * Fs;
        fbi = f; fbi(f > Fs/2) = f(f > Fs/2) - Fs;       % two-sided axis
        y = real(ifft(X .* gaussMask(fbi, f1, f2)));
    end
