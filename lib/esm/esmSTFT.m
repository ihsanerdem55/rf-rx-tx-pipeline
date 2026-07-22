    function [Tc,F,Sdb]=esmSTFT(x,Fs,iq)
        x=x(:).'; N=numel(x); L=min(384,max(128,2^round(log2(max(64,N/180)))));
        L=min(L,N); H=max(1,round(L/2)); Nfft=max(512,2^nextpow2(L));
        w=0.5-0.5*cos(2*pi*(0:L-1)/max(1,L-1)); starts=1:H:max(1,N-L+1);
        if iq
            bins=1:Nfft; F=((-Nfft/2):(Nfft/2-1))/Nfft*Fs;
        else
            bins=1:floor(Nfft/2)+1; F=(bins-1)/Nfft*Fs;
        end
        Sdb=zeros(numel(bins),numel(starts));
        for k=1:numel(starts)
            ix=starts(k):min(N,starts(k)+L-1); seg=x(ix);
            if numel(seg)<L, seg(end+1:L)=0; end
            Z=abs(fft(seg.*w,Nfft));
            if iq, Z=fftshift(Z); end
            Sdb(:,k)=Z(bins).';
        end
        Sdb=20*log10(Sdb/(max(Sdb(:))+eps)+eps); Sdb=max(Sdb,-55);
        Tc=(starts-1+L/2)/Fs;
    end
