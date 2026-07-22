    function [w,nrej]=priConsistencyWeights(rc)
        % Odd/even pulses correspond to the two stagger values. The harmonic
        % agreement metric is one for equal power and approaches zero when a
        % cell exists in only one stagger group.
        Nr=size(rc,1); Np=size(rc,2); w=ones(Nr,1); nrej=0;
        ia=1:2:Np; ib=2:2:Np;
        if isempty(ia) || isempty(ib), return, end
        pa=mean(abs(rc(:,ia)).^2,2); pb=mean(abs(rc(:,ib)).^2,2);
        agree=2*min(pa,pb)./(pa+pb+eps);
        energetic=max(pa,pb)>median(max(pa,pb))+4*median(abs(max(pa,pb)-median(max(pa,pb))));
        nrej=sum(energetic & agree<0.25);
        w=max(0,min(1,agree)).^1.5;
    end
