    function m=localMax2(A)
        m=true(size(A));
        for dr=-1:1
            for dc=-1:1
                if dr==0 && dc==0, continue, end
                m=m & A>=circshift(A,[dr dc]);
            end
        end
        if size(A,1)>1, m([1 end],:)=false; end
        if size(A,2)>1, m(:,[1 end])=false; end
    end
