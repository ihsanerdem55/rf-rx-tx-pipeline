    function [RD,fdAxis]=slowTimeDFT(rc,tSlow,Fc,PRI,Ndop)
        Np=size(rc,2);
        if nargin<5, Ndop=max(128,2^nextpow2(8*max(Np,1))); end
        PRF=1/PRI; fdAxis=((0:Ndop-1)-Ndop/2)/Ndop*PRF;
        if Np<4, sw=ones(1,Np); else, sw=0.5-0.5*cos(2*pi*(0:Np-1)/(Np-1)); end
        E=exp(-1i*2*pi*tSlow(:)*fdAxis);
        RD=(rc.*sw)*E;
    end
