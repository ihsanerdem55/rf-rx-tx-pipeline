    function y = addNoiseJammer(y,P,seed)
        % Spot noise occupies a narrow band about the captured Tx IF; barrage
        % noise fills the sampled band. J/S is referenced to average replay power.
        if strcmp(P.jammode,'Off') || isempty(y), return, end
        N=numel(y); rs=RandStream('mt19937ar','Seed',seed);
        w=(randn(rs,1,N)+1i*randn(rs,1,N))/sqrt(2);
        if strcmp(P.jammode,'Spot')
            W=fft(w); f=((0:N-1)-floor(N/2))/N*P.Fs;
            W=fftshift(W); W(abs(f)>P.spotbw/2)=0; w=ifft(ifftshift(W));
        end
        t=(0:N-1)/P.Fs;
        jam=w.*exp(1i*2*pi*P.Ft*t);
        jam=jam/sqrt(mean(abs(jam).^2)+eps);
        Ps=max(mean(abs(y).^2),eps);
        jam=jam*sqrt(Ps*10^(P.jsdb/10));
        if isreal(y), y=y+real(jam); else, y=y+jam; end
    end
