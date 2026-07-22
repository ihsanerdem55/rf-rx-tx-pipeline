    function y = fftLowpass(x, Fs, Wc, shape, tw)
        % Zero-phase complex low-pass by FFT masking (keeps complex I/Q).
        x = x(:).'; N = numel(x);
        X = fft(x);
        f = (0:N-1) / N * Fs;
        fbi = f; fbi(f > Fs/2) = f(f > Fs/2) - Fs;   % two-sided axis
        a = abs(fbi);
        switch shape
            case 'gauss'
                sigma = Wc / sqrt(log(2));            % -3 dB at Wc
                m = exp(-0.5 * (a/sigma).^2);
            case 'cos'
                m = double(a <= Wc);
                if tw > 0
                    edge = (a > Wc) & (a <= Wc+tw);
                    m(edge) = 0.5*(1 + cos(pi*(a(edge)-Wc)/tw));
                end
            otherwise  % brick
                m = double(a <= Wc);
        end
        y = ifft(X .* m);                            % complex in, complex out
    end
