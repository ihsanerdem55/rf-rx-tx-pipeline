    function H = freqzLite(b, a, f, Fs)
        % Frequency response of b/a evaluated at frequencies f (Hz).
        % (freqz() is toolbox-only, so evaluate the polynomials by hand.)
        f = f(:).';
        w = 2*pi*f/Fs;
        B = zeros(size(w)); A = zeros(size(w));
        for k = 1:numel(b), B = B + b(k)*exp(-1i*w*(k-1)); end
        for k = 1:numel(a), A = A + a(k)*exp(-1i*w*(k-1)); end
        H = B ./ A;
    end
