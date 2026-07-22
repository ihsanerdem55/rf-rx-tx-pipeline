    function jam=victimJammer(N,Fs,P,hop,rs)
        w=(randn(rs,N,1)+1i*randn(rs,N,1))/sqrt(2);
        if strcmp(P.jammode,'Spot')
            W=fft(w); f=((0:N-1)-floor(N/2)).'/N*Fs; W=fftshift(W);
            W(abs(f)>P.spotbw/2)=0; w=ifft(ifftshift(W));
            t=(0:N-1).'/Fs;
            jam=w.*exp(-1i*2*pi*hop*t); % victim hopping moves captured spot off channel
        else
            jam=w;
        end
    end
