    function [mask,threshold]=caCfar2D(Pwr,pfa)
        [Nr,Nd]=size(Pwr); Tr=min(10,max(2,floor(Nr/20))); Gr=2;
        Td=min(8,max(2,floor(Nd/24))); Gd=1;
        outer=ones(2*(Tr+Gr)+1,2*(Td+Gd)+1);
        guard=zeros(size(outer));
        guard(Tr+1:Tr+2*Gr+1,Td+1:Td+2*Gd+1)=1;
        train=outer-guard; Ntrain=sum(train(:));
        noise=conv2(Pwr,train,'same')/max(Ntrain,1);
        alpha=Ntrain*(pfa^(-1/max(Ntrain,1))-1);
        threshold=alpha*noise; mask=Pwr>threshold;
        er=Tr+Gr; ed=Td+Gd;
        mask(1:min(er,Nr),:)=false; mask(max(1,Nr-er+1):Nr,:)=false;
        mask(:,1:min(ed,Nd))=false; mask(:,max(1,Nd-ed+1):Nd)=false;
    end
