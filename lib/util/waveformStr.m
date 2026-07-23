    function s = waveformStr(P)
        % One-line description of the transmitted ground-truth waveform.
        if strcmp(P.modtype,'LFM chirp')
            s=sprintf('LFM %.2f MHz | %s PRI | %s',P.bchirp/1e6,P.primode,P.fhmode);
        elseif strcmp(P.modtype,'NLFM chirp')
            s=sprintf('NLFM %.2f MHz | %s PRI | %s',P.bchirp/1e6,P.primode,P.fhmode);
        else
            s=sprintf('CW pulses | %s PRI | %s',P.primode,P.fhmode);
        end
    end
