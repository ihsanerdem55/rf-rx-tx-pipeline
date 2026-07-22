    function M = runESMMUSIC(P,X,t,pdw,labels,tracks)
        % Narrowband MUSIC is evaluated independently for each deinterleaved
        % frequency group. A three-sensor array leaves a two-dimensional signal
        % subspace at most, hence the explicit two co-channel source cap.
        M=struct('groups',[],'trackDoa',NaN(1,numel(tracks)),'grid',-90:P.musicStep:90);
        if isempty(tracks) || isempty(pdw), return; end
        K=numel(tracks); used=false(1,K); groups={};
        for k=1:K
            if used(k), continue; end
            g=k; used(k)=true;
            for j=k+1:K
                tol=max(0.30e6,0.18*(tracks(k).freqSpan+tracks(j).freqSpan)+0.12e6);
                if abs(tracks(j).freqMed-tracks(k).freqMed)<tol
                    g(end+1)=j; used(j)=true; %#ok<AGROW>
                end
            end
            groups{end+1}=g; %#ok<AGROW>
        end
        pos=P.esmPos(:);
        G=repmat(struct('members',[],'freq',0,'pdb',[],'peaks',[]),0,1);
        for gi=1:numel(groups)
            members=groups{gi}; ids=find(ismember(labels,members)); V=[];
            for q=1:numel(ids)
                d=pdw(ids(q)); ic=round((d.toa+0.5*d.pw)*P.Fs)+1;
                half=max(3,round(min(0.22*d.pw,0.5e-6)*P.Fs/2));
                ix=max(1,ic-half):min(numel(t),ic+half);
                if isempty(ix), continue; end
                v=mean(X(:,ix).*exp(-1i*2*pi*d.freq*t(ix)),2);
                if norm(v)>0, V(:,end+1)=v/norm(v); end %#ok<AGROW>
            end
            if size(V,2)<2, continue; end
            R=V*V'/size(V,2); R=(R+R')/2+1e-6*trace(R)/3*eye(3);
            [U,D]=eig(R); [~,ord]=sort(real(diag(D)),'descend'); U=U(:,ord);
            ns=min(2,numel(members)); En=U(:,ns+1:end);
            fref=median([tracks(members).freqMed]); rf=P.Flo+fref;
            ps=zeros(size(M.grid));
            for a=1:numel(M.grid)
                sv=exp(1i*2*pi*pos*sind(M.grid(a))/(P.c/max(rf,eps)));
                ps(a)=1/(sum(abs(En'*sv).^2)+eps);
            end
            pdb=10*log10(ps/max(ps)+eps);
            peaks=findMusicPeaks(pdb,M.grid,ns);
            z=struct('members',members,'freq',fref,'pdb',pdb,'peaks',peaks);
            G(end+1)=z; %#ok<AGROW>
            % Assign MUSIC peaks to tracks, but only as a REFINEMENT: a co-channel
            % group with more emitters than a 3-sensor array can resolve produces
            % fewer/blended peaks, so a peak is accepted only when it corroborates
            % the robust single-source beamformer DOA (median over the track's own
            % pulses). Otherwise the beamformer angle is kept -- which is correct,
            % because the deinterleaver has already isolated each emitter.
            avail=peaks; refineTol=6;   % deg
            for m=members
                if isempty(avail), break; end
                [dmin,ii]=min(abs(avail-tracks(m).doa));
                if dmin<=refineTol, M.trackDoa(m)=avail(ii); avail(ii)=[]; end
            end
        end
        M.groups=G;
    end
