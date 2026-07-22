    function [f, magdB] = spectrumTwoSided(x, Fs)
        % Always-two-sided magnitude spectrum in dB (peak-normalised). Used for
        % the transmit MIXER view: a real IF signal has symmetric +/- content,
        % so shifting this two-sided spectrum up by the LO reveals BOTH the
        % wanted sideband and its IMAGE. (A complex/IQ signal is one-sided, so
        % its image is naturally absent -- exactly what an IQ upconverter does.)
        x = x(:).'; N = numel(x);
        w = 0.5 - 0.5*cos(2*pi*(0:N-1)/max(1,N-1));   % Hann
        Nfft = 2^nextpow2(max(1024, N));
        X = abs(fftshift(fft(x .* w, Nfft)));
        f = ((0:Nfft-1) - Nfft/2) / Nfft * Fs;
        magdB = 20*log10(X / max(X + eps) + eps);
    end
