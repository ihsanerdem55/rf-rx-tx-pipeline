    function [f, magdB] = spectrum(x, Fs)
        % Magnitude spectrum in dB (peak-normalised). One-sided for a real (IF)
        % signal, two-sided (fftshift) for a complex baseband signal.
        x = x(:).'; N = numel(x);
        w = 0.5 - 0.5*cos(2*pi*(0:N-1)/max(1,N-1));   % Hann
        Nfft = 2^nextpow2(max(1024, N));
        if isreal(x)
            X = abs(fft(x .* w, Nfft));
            half = 1:floor(Nfft/2)+1;
            f = (half-1) / Nfft * Fs;
            mag = X(half);
        else
            X = abs(fftshift(fft(x .* w, Nfft)));
            f = ((0:Nfft-1) - Nfft/2) / Nfft * Fs;
            mag = X;
        end
        magdB = 20*log10(mag / max(mag + eps) + eps);
    end
