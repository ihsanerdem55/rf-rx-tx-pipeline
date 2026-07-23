function tf = isChirp(P)
% ISCHIRP  True when the transmitted waveform is a swept chirp (LFM or NLFM),
%   i.e. it carries an intra-pulse frequency modulation that supports pulse
%   compression. CW returns false. Centralises the waveform-family test so the
%   generator, the matched filters and the validation share one definition.
    tf = any(strcmp(P.modtype, {'LFM chirp', 'NLFM chirp'}));
end
