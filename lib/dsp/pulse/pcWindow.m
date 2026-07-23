function w = pcWindow(L, kind)
% PCWINDOW  Real amplitude taper (length L) for a matched-filter reference.
%   Tapering the pulse-compression reference trades compressed-pulse SNR and
%   mainlobe width for lower range sidelobes -- the classic weighting
%   trade-off. 'None' returns a flat window, i.e. the ordinary (matched)
%   filter. Base MATLAB only.
%
%       'None'      flat            -> matched filter, ~ -13 dB LFM sidelobes
%       'Hamming'   0.54/0.46       -> ~ -42 dB sidelobes, slight SNR loss
%       'Blackman'  0.42/0.5/0.08   -> ~ -58 dB sidelobes, more SNR/width cost
    if nargin < 2 || isempty(kind), kind = 'None'; end
    L = max(1, round(L));
    if L < 2
        w = ones(1, L);
        return
    end
    n = (0:L-1) / (L-1);
    switch kind
        case 'Hamming'
            w = 0.54 - 0.46*cos(2*pi*n);
        case 'Blackman'
            w = 0.42 - 0.5*cos(2*pi*n) + 0.08*cos(4*pi*n);
        otherwise                         % 'None' (or unknown) -> flat
            w = ones(1, L);
    end
end
