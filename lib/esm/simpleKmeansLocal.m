    function labels = simpleKmeansLocal(F,K)
        n=size(F,1); C=zeros(K,size(F,2)); C(1,:)=F(1,:);
        for k=2:K
            d=inf(n,1);
            for q=1:k-1, d=min(d,sum((F-C(q,:)).^2,2)); end
            [~,ix]=max(d); C(k,:)=F(ix,:);
        end
        labels=ones(n,1);
        for it=1:50
            D=zeros(n,K);
            for k=1:K, D(:,k)=sum((F-C(k,:)).^2,2); end
            [~,new]=min(D,[],2);
            if it>1 && all(new==labels), break; end
            labels=new;
            for k=1:K
                if any(labels==k), C(k,:)=mean(F(labels==k,:),1);
                else
                    [~,ix]=max(min(D,[],2)); C(k,:)=F(ix,:); labels(ix)=k;
                end
            end
        end
    end
