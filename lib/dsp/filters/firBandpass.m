    function hf = firBandpass(Ntaps, f1, f2, Fs, wintype, beta)
        % Ideal band-pass = difference of two ideal low-passes. The ideal
        % impulse response is infinite, so truncate it to Ntaps and apply a
        % window to tame the truncation ripple. Two windows are offered:
        %   'hamming' : classic, ~ -53 dB sidelobes.
        %   'kaiser'  : adjustable beta; larger beta -> lower sidelobes but a
        %               wider transition band (beta~8 gives ~ -84 dB).
        if nargin < 5, wintype = 'hamming'; end
        if nargin < 6, beta = 8; end
        if mod(Ntaps,2) == 0, Ntaps = Ntaps + 1; end     % odd -> Type-I
        M = (Ntaps-1)/2;
        n = (0:Ntaps-1) - M;                             % symmetric index
        hf = (2*f2/Fs)*sincLocal(2*f2*n/Fs) - (2*f1/Fs)*sincLocal(2*f1*n/Fs);
        switch wintype
            case 'kaiser'
                r = ((0:Ntaps-1) - M) / M;               % -1..+1
                w = besseli(0, beta*sqrt(max(0,1-r.^2))) / besseli(0, beta);
            otherwise    % 'hamming'
                w = 0.54 - 0.46*cos(2*pi*(0:Ntaps-1)/(Ntaps-1));
        end
        hf = hf .* w;
        % normalise to unity gain at the geometric centre of the band
        fc = sqrt(max(f1, eps) * f2);
        H  = sum(hf .* exp(-2j*pi*fc*(0:Ntaps-1)/Fs));
        if abs(H) > 0, hf = hf / abs(H); end
    end
