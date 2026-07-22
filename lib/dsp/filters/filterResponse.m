    function [Hmag, himp] = filterResponse(P, f)
        % Effective magnitude response at frequencies f (Hz), plus the impulse
        % response where one exists (empty for the FFT-mask filters).
        himp = [];
        if P.iq
            % baseband low-pass response (|f| axis), half-width Wc
            Wc = 0.5*(P.Bhi - P.Blo); a = abs(f);
            switch P.filt
                case 'FFT Gaussian'
                    sigma = Wc/sqrt(log(2)); Hmag = exp(-0.5*(a/sigma).^2);
                case 'FFT cosine-taper'
                    tw = P.tap*2*Wc; Hmag = double(a<=Wc);
                    e = (a>Wc)&(a<=Wc+tw); Hmag(e) = 0.5*(1+cos(pi*(a(e)-Wc)/tw));
                case 'FFT brick-wall'
                    Hmag = double(a <= Wc);
                otherwise    % FIR/IIR fall back to cosine low-pass
                    tw = 0.2*2*Wc; Hmag = double(a<=Wc);
                    e = (a>Wc)&(a<=Wc+tw); Hmag(e) = 0.5*(1+cos(pi*(a(e)-Wc)/tw));
            end
            return
        end
        switch P.filt
            case 'FFT cosine-taper'
                tw   = P.tap * (P.Bhi - P.Blo);
                Hmag = cosTaper(f, P.Blo, P.Bhi, tw);
            case 'FFT Gaussian'
                Hmag = gaussMask(f, P.Blo, P.Bhi);
            case 'FIR windowed-sinc'
                himp = firBandpass(P.ord, P.Blo, P.Bhi, P.Fs, 'hamming');
                Hmag = abs(freqzLite(himp, 1, f, P.Fs));   % linear phase: |H|
            case 'FIR Kaiser'
                himp = firBandpass(P.ord, P.Blo, P.Bhi, P.Fs, 'kaiser', P.kbeta);
                Hmag = abs(freqzLite(himp, 1, f, P.Fs));
            case 'Butterworth IIR'
                [b, a] = butterBandpass(P.ord, P.Blo, P.Bhi, P.Fs);
                himp = impzLite(b, a, 200);
                % forward-backward pass => the signal actually sees |H|^2
                Hmag = abs(freqzLite(b, a, f, P.Fs)).^2;
            case 'Chebyshev I IIR'
                [b, a] = cheby1Bandpass(P.ord, P.Blo, P.Bhi, P.Fs, P.crip);
                himp = impzLite(b, a, 200);
                Hmag = abs(freqzLite(b, a, f, P.Fs)).^2;   % |H|^2 (zero-phase)
            otherwise    % 'FFT brick-wall'
                Hmag = double(f >= P.Blo & f <= P.Bhi);
        end
    end
