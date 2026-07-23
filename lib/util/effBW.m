function B = effBW(P)
% EFFBW  Effective occupied bandwidth (Hz) used for range resolution, the
%   time-bandwidth product and compressed-pulse diagnostics. For a swept chirp
%   this is the sweep bandwidth; for a phase code it is the chip rate (code
%   length / pulse width), which is what actually sets the compressed range
%   resolution. CW (and a zero-bandwidth chirp) returns 0.
    switch P.modtype
        case {'LFM chirp', 'NLFM chirp'}
            B = P.bchirp;
        case 'Barker-13 code'
            B = 13 / P.PW;              % chip rate = code length / pulse width
        otherwise
            B = 0;
    end
end
