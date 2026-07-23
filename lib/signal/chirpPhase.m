function phi = chirpPhase(tl, P)
% CHIRPPHASE  Intra-pulse phase (radians) of the transmitted waveform,
%   evaluated at the pulse-relative sample times tl (tl = 0 at the pulse
%   centre).
%
%   The phase is returned as a function of tl, so every caller -- the signal
%   generator and each matched-filter reference -- stays mutually consistent
%   regardless of its own sample grid. Base MATLAB only (fft-free).
%
%     'LFM chirp'      linear sweep, phi = pi*(B/PW)*tl.^2
%     'NLFM chirp'     Hamming-weighted stationary-phase design: the frequency
%                      dwells longer near band centre, so the compressed
%                      spectrum is tapered and range sidelobes fall to about
%                      -24 dB (TB ~ 30) without the SNR loss of a receive
%                      window.
%     'Barker-13 code' bi-phase code: phi is 0 or pi per chip; compresses to a
%                      chip-wide peak with ~ -22 dB sidelobes (peak 1/13).
%     other (CW)       phi = 0
    switch P.modtype
        case 'LFM chirp'
            if P.bchirp > 0
                phi = pi*(P.bchirp/P.PW)*tl.^2;
            else
                phi = zeros(size(tl));
            end
        case 'NLFM chirp'
            if P.bchirp > 0
                phi = nlfmPhase(tl, P.bchirp, P.PW);
            else
                phi = zeros(size(tl));
            end
        case 'Barker-13 code'
            phi = barkerPhase(tl, P.PW);
        otherwise
            phi = zeros(size(tl));
    end
end

% =====================================================================
function phi = nlfmPhase(tl, B, PW)
% Numerical stationary-phase NLFM: choose the group delay t(f) so the time
% spent near each frequency follows a Hamming power weighting, giving a
% tapered spectrum and low matched-filter sidelobes with no receive window.
    Ng = 2048;
    fg = linspace(-B/2, B/2, Ng);
    W  = 0.54 + 0.46*cos(2*pi*fg/B);          % Hamming weighting across the band
    C  = cumtrapz(fg, W.^2);
    Tf = PW*(C/C(end)) - PW/2;                % group delay t(f) in [-PW/2, PW/2]
    tg = linspace(-PW/2, PW/2, Ng);
    fT = interp1(Tf, fg, tg, 'linear', 'extrap');   % invert to instantaneous f(t)
    phiG = 2*pi*cumtrapz(tg, fT);             % integrate frequency -> phase(t)
    phi = interp1(tg, phiG, tl, 'linear', 'extrap');
    phi = reshape(phi, size(tl));
end

% =====================================================================
function phi = barkerPhase(tl, PW)
% Bi-phase Barker-13 code: 0 or pi per chip across the pulse width.
    code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    N = numel(code);
    idx = floor((tl + PW/2)/PW * N) + 1;      % chip index for each sample
    idx = min(N, max(1, idx));
    phi = double(code(idx) < 0) * pi;
    phi = reshape(phi, size(tl));
end
