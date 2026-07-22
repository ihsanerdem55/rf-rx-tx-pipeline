    function y = applyFilter(x, P)
        % Dispatch to the selected filter. In BASEBAND the signal is complex and
        % centred at 0, so the receiver filter is an FFT LOW-PASS of half-width
        % Wc = (Bhi-Blo)/2 (the FFT-shape choice is honoured; FIR/IIR band-pass
        % designs have no baseband meaning and fall back to a cosine low-pass).
        if P.iq
            Wc = 0.5*(P.Bhi - P.Blo);
            switch P.filt
                case 'FFT Gaussian',    y = fftLowpass(x, P.Fs, Wc, 'gauss', 0);
                case 'FFT cosine-taper', y = fftLowpass(x, P.Fs, Wc, 'cos', P.tap*2*Wc);
                case 'FFT brick-wall',  y = fftLowpass(x, P.Fs, Wc, 'brick', 0);
                otherwise,              y = fftLowpass(x, P.Fs, Wc, 'cos', 0.2*2*Wc);
            end
            return
        end
        switch P.filt
            case 'FFT cosine-taper'
                tw = P.tap * (P.Bhi - P.Blo);           % taper width in Hz
                y  = fftMask(x, P.Fs, P.Blo, P.Bhi, tw);
            case 'FFT Gaussian'
                y  = fftGauss(x, P.Fs, P.Blo, P.Bhi);
            case 'FIR windowed-sinc'
                hf = firBandpass(P.ord, P.Blo, P.Bhi, P.Fs, 'hamming');
                y  = firZeroPhase(x, hf);
            case 'FIR Kaiser'
                hf = firBandpass(P.ord, P.Blo, P.Bhi, P.Fs, 'kaiser', P.kbeta);
                y  = firZeroPhase(x, hf);
            case 'Butterworth IIR'
                [b, a] = butterBandpass(P.ord, P.Blo, P.Bhi, P.Fs);
                y  = filtfiltLite(b, a, x);
            case 'Chebyshev I IIR'
                [b, a] = cheby1Bandpass(P.ord, P.Blo, P.Bhi, P.Fs, P.crip);
                y  = filtfiltLite(b, a, x);
            otherwise    % 'FFT brick-wall'
                y  = fftMask(x, P.Fs, P.Blo, P.Bhi, 0);
        end
    end
