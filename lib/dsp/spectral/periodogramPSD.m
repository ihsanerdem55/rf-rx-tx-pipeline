    function [f, pxx] = periodogramPSD(x, Fs, Nfft, twosided)
        % Periodogram PSD from a single Hann-windowed, zero-padded FFT.
        if nargin < 4, twosided = false; end
        x = x(:).'; N = numel(x);
        win = 0.5 - 0.5*cos(2*pi*(0:N-1)/max(1,N-1));    % Hann
        U = sum(win.^2);
        X = fft(x .* win, Nfft);
        if twosided
            Xs = fftshift(X);
            f = ((0:Nfft-1) - Nfft/2) / Nfft * Fs;
            pxx = abs(Xs).^2 / (Fs * U);
        else
            half = 1:floor(Nfft/2)+1;
            f = (half-1) / Nfft * Fs;
            pxx = abs(X(half)).^2 / (Fs * U);
            pxx(2:end-1) = 2 * pxx(2:end-1);              % one-sided fold
        end
    end
