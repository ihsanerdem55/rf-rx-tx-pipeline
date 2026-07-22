    function [f, pxx] = welchPSD(x, Fs, nseg, overlap, Nfft, twosided)
        % Welch's method (hand-rolled): overlapping Hann segments, FFT each,
        % average magnitude-squared. One-sided (real) or two-sided (complex).
        if nargin < 6, twosided = false; end
        x = x(:).'; N = numel(x);
        ov = min(0.95, max(0, overlap));
        L = floor(N / (1 + (nseg-1)*(1-ov)));
        L = max(16, min(N, L));
        step = max(1, round(L*(1-ov)));
        win = 0.5 - 0.5*cos(2*pi*(0:L-1)/max(1,L-1));     % Hann
        U = sum(win.^2);

        if twosided
            f = ((0:Nfft-1) - Nfft/2) / Nfft * Fs;
            acc = zeros(1, Nfft);
        else
            half = 1:floor(Nfft/2)+1;
            f = (half-1) / Nfft * Fs;
            acc = zeros(1, numel(half));
        end
        cnt = 0; start = 1;
        while start + L - 1 <= N
            seg = x(start:start+L-1) .* win;
            X = fft(seg, Nfft);
            if twosided
                acc = acc + abs(fftshift(X)).^2;
            else
                acc = acc + abs(X(half)).^2;
            end
            cnt = cnt + 1;
            start = start + step;
        end
        if cnt == 0
            [f, pxx] = periodogramPSD(x, Fs, Nfft, twosided);
            return
        end
        pxx = acc / cnt / (Fs * U);
        if ~twosided, pxx(2:end-1) = 2 * pxx(2:end-1); end
    end
