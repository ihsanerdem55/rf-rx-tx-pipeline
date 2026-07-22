    function peaks=findMusicPeaks(pdb,grid,n)
        cand=find(pdb(2:end-1)>=pdb(1:end-2) & pdb(2:end-1)>pdb(3:end))+1;
        if isempty(cand), [~,cand]=max(pdb); end
        [~,ord]=sort(pdb(cand),'descend'); cand=cand(ord); keep=[];
        for k=1:numel(cand)
            if isempty(keep) || all(abs(grid(cand(k))-grid(keep))>3)
                keep(end+1)=cand(k); %#ok<AGROW>
                if numel(keep)>=n, break; end
            end
        end
        peaks=grid(keep);
    end
