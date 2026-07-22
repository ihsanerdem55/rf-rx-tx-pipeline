    function [Tc, F, Sdb] = stftLocal(x, Fs)
        % Hand-rolled STFT. Real IF data is displayed one-sided; complex IQ is
        % displayed two-sided so negative baseband hops are not hidden.
        x=x(:).'; N=numel(x);
        L=max(64,2^round(log2(max(N/64,1)))); L=min(L,N);
        H=max(1,round(L/2)); Nfft=max(256,2^nextpow2(L));
        w=0.5-0.5*cos(2*pi*(0:L-1)/max(1,L-1));
        starts=1:H:(N-L+1); if isempty(starts), starts=1; end
        if isreal(x)
            bins=1:floor(Nfft/2)+1; F=(bins-1)/Nfft*Fs;
        else
            bins=1:Nfft; F=((0:Nfft-1)-Nfft/2)/Nfft*Fs;
        end
        Sdb=zeros(numel(bins),numel(starts));
        for c=1:numel(starts)
            seg=x(starts(c):min(N,starts(c)+L-1));
            if numel(seg)<L, seg(end+1:L)=0; end
            X=fft(seg.*w,Nfft);
            if isreal(x), v=abs(X(bins)); else, v=abs(fftshift(X)); end
            Sdb(:,c)=v(:);
        end
        Sdb=20*log10(Sdb/(max(Sdb(:))+eps)+eps); Sdb=max(Sdb,-60);
        Tc=(starts-1+L/2)/Fs;
    end
