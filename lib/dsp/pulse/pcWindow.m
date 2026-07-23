function w = pcWindow(L, kind)
% PCWINDOW  Real amplitude taper (length L) for a matched-filter reference.
%   Tapering the pulse-compression reference trades compressed-pulse SNR and
%   mainlobe width for lower range sidelobes -- the classic weighting
%   trade-off. 'None' returns a flat window, i.e. the ordinary (matched)
%   filter. Base MATLAB only (no Signal Processing Toolbox).
%
%   Applied to the RECEIVE reference only (a mismatched filter -- the signal
%   itself is still plain LFM), so the achievable sidelobe floor is limited to
%   about -30 dB, unlike transmit+receive weighting. Compressed-LFM peak
%   sidelobe level / mainlobe broadening (TB ~ 30):
%       'None'      flat               -14 dB,  x1.0   (matched filter)
%       'Taylor'    n-bar = 4, -30 dB  -26 dB,  x1.4   (least broadening)
%       'Hann'      0.5/0.5            -31 dB,  x1.7
%       'Hamming'   0.54/0.46          -32 dB,  x1.6
%       'Blackman'  0.42/0.5/0.08      -31 dB,  x1.9   (widest lobe / most loss)
    if nargin < 2 || isempty(kind), kind = 'None'; end
    L = max(1, round(L));
    if L < 2
        w = ones(1, L);
        return
    end
    n = (0:L-1) / (L-1);
    switch kind
        case 'Hann'
            w = 0.5 - 0.5*cos(2*pi*n);
        case 'Hamming'
            w = 0.54 - 0.46*cos(2*pi*n);
        case 'Blackman'
            w = 0.42 - 0.5*cos(2*pi*n) + 0.08*cos(4*pi*n);
        case 'Taylor'
            w = taylorWin(L, 4, -30);
        otherwise                         % 'None' (or unknown) -> flat
            w = ones(1, L);
    end
end

function w = taylorWin(L, nbar, SLL)
% Taylor n-bar window: designed for a flat near-in sidelobe floor at SLL dB
% (SLL negative). Standard closed form, implemented in base MATLAB.
    R   = 10^(-SLL/20);
    A   = acosh(R)/pi;
    sp2 = nbar^2 / (A^2 + (nbar-0.5)^2);   % scaling factor sigma^2
    Fm  = zeros(1, nbar-1);
    for m = 1:nbar-1
        num = 1; den = 1;
        for i = 1:nbar-1
            num = num * (1 - m^2/(sp2*(A^2 + (i-0.5)^2)));
            if i ~= m
                den = den * (1 - m^2/i^2);
            end
        end
        Fm(m) = ((-1)^(m+1) * num) / (2*den);
    end
    x = ((0:L-1) - (L-1)/2) / L;
    w = ones(1, L);
    for m = 1:nbar-1
        w = w + 2*Fm(m)*cos(2*pi*m*x);
    end
end
