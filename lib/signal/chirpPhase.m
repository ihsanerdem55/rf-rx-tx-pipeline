function phi = chirpPhase(tl, P)
% CHIRPPHASE  Intra-pulse phase (radians) of the transmitted chirp, evaluated
%   at the pulse-relative sample times tl (tl = 0 at the pulse centre).
%
%   The phase is returned as a function of tl, so every caller -- the signal
%   generator and each matched-filter reference -- stays mutually consistent
%   regardless of its own sample grid. Base MATLAB only (fft-free: cumtrapz,
%   interp1, linspace).
%
%     'LFM chirp'   linear sweep,  phi = pi*(B/PW)*tl.^2
%     'NLFM chirp'  Hamming-weighted stationary-phase design: the instantaneous
%                   frequency dwells longer near band centre, so the compressed
%                   spectrum is tapered and range sidelobes fall to about
%                   -24 dB (TB ~ 30) WITHOUT the SNR loss of a receive window.
%     other (CW)    phi = 0
    B = P.bchirp; PW = P.PW;
    if ~isChirp(P) || ~(B > 0)
        phi = zeros(size(tl));
        return
    end
    if strcmp(P.modtype, 'LFM chirp')
        phi = pi*(B/PW)*tl.^2;
        return
    end

    % --- NLFM: numerical stationary-phase design -----------------------------
    % Choose the group delay t(f) so the time spent near each frequency follows
    % a Hamming power weighting W(f)^2; the resulting spectrum is tapered and
    % the matched-filter sidelobes drop without any receive weighting.
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
