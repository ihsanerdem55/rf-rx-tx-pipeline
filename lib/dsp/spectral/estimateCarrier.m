    function [fc, fBand, magBandDB] = estimateCarrier(x, P)
        % Carrier = strongest in-band spectral peak, refined to sub-bin
        % resolution by 3-point parabolic interpolation. Two spectral estimators
        % are available (both hand-rolled -- pwelch is Signal Processing Toolbox):
        %
        %   'Periodogram (FFT)' : one windowed, zero-padded FFT of the whole
        %       record. Highest frequency resolution (~Fs/N).
        %   'Welch PSD'         : average the periodograms of overlapping Hann
        %       segments -- lower variance, coarser resolution.
        %
        % For BASEBAND the PSD is two-sided (the carrier sits near 0 Hz and can
        % be either sign); for IF it is one-sided over [Blo, Bhi].
        Fs = P.Fs;
        N = numel(x);
        Nfft = 2^nextpow2(max(4096, 4*N));   % common interpolation grid
        ts = P.iq;                           % two-sided for complex baseband

        switch P.cest
            case 'Welch PSD'
                [f, pxx] = welchPSD(x, Fs, P.wseg, P.wov, Nfft, ts);
            otherwise    % 'Periodogram (FFT)'
                [f, pxx] = periodogramPSD(x, Fs, Nfft, ts);
        end

        if ts
            Wc = 0.5*(P.Bhi - P.Blo);
            inband = abs(f) <= Wc;
        else
            inband = f >= P.Blo & f <= P.Bhi;
        end
        if ~any(inband)                                  % degenerate band
            fc = NaN; fBand = f;
            magBandDB = 10*log10(pxx / max(pxx + eps) + eps);
            return
        end
        pm = pxx; pm(~inband) = -inf;
        [~, im] = max(pm);
        fc = f(im);
        if im > 1 && im < numel(pxx)
            aL = 10*log10(pxx(im-1) + eps);
            bC = 10*log10(pxx(im)   + eps);
            cR = 10*log10(pxx(im+1) + eps);
            den = aL - 2*bC + cR;
            if den ~= 0
                delta = 0.5 * (aL - cR) / den;           % peak offset, in bins
                delta = max(-0.5, min(0.5, delta));
                fc = fc + delta * (Fs / Nfft);           % bin spacing = Fs/Nfft
            end
        end

        fBand = f(inband);
        pb = pxx(inband);
        magBandDB = 10*log10(pb / max(pb + eps) + eps);  % peak-normalised PSD (dB)
    end
