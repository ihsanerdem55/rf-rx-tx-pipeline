function tf = isCompressible(P)
% ISCOMPRESSIBLE  True when the transmitted waveform carries an intra-pulse
%   modulation with a matched-filter (pulse-compression) reference -- a swept
%   chirp (LFM/NLFM) OR a phase code (Barker). CW returns false. Use this to
%   gate pulse compression; use ISCHIRP only where a *frequency sweep* is
%   specifically meant (e.g. clamping the sweep bandwidth field).
    tf = any(strcmp(P.modtype, {'LFM chirp', 'NLFM chirp', 'Barker-13 code'}));
end
