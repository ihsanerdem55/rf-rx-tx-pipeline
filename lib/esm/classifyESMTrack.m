    function c=classifyESMTrack(T)
        pri=T.PRI; pw=T.pw;
        if pri>0 && pri<70e-6 && pw<4e-6
            c='Missile seeker';
        elseif pri>135e-6 && pw>=4e-6
            c='Search radar';
        elseif T.priAgile || T.freqAgile
            c='Agile multifunction';
        elseif pri>0 && pri<=140e-6
            c='Track radar';
        else
            c='Unknown emitter';
        end
    end
