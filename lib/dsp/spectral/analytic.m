    function z = analytic(x)
        % Analytic signal (Hilbert transform) built with the FFT -- no toolbox.
        % Keep DC & Nyquist, double the positive frequencies, zero the negative
        % ones. |z| is then the amplitude envelope of x.
        x = x(:).'; N = numel(x);
        X = fft(x); Hf = zeros(1, N);
        if mod(N, 2) == 0
            Hf([1 N/2+1]) = 1; Hf(2:N/2) = 2;
        else
            Hf(1) = 1; Hf(2:(N+1)/2) = 2;
        end
        z = ifft(X .* Hf);
    end
