classdef Plot < handle
properties(Constant,Hidden)
    colorsRE='[rgbcmykw]?';
    markersRE='[o+*.x_|sd^v><ph]?';
    linesRE='(-\.|-{1,2}|:)';
    colorkey=struct('r',[1 0 0], ...
                    'g',[0 1 0], ...
                    'b',[0 0 1], ...
                    'c',[0 1 0], ...
                    'm',[1 0 1], ...
                    'y',[1 1 0], ...
                    'k',[0 0 0], ...
                    'w',[1 1 1])
end
methods(Static)
    function arrowTest_()
        x=[0 1];
        y=[0 10];
        figure(1)
        hold off;
        len=[1 90];
        nArrow=2;
        Plot.arrow([0  10], [0  10],len,nArrow); hold on
        %Plot.arrow([0 -10], [0  10],len,nArrow);
        %Plot.arrow([0 -10], [0 -10],len,nArrow);
        %Plot.arrow([0  10], [0 -10],len,nArrow);

        Plot.arrow([0   0], [0  10],len,nArrow);
        Plot.arrow([0  10], [0   0],len,nArrow);
        Plot.arrow([0 -10], [0   0],len,nArrow);
        Plot.arrow([0   0], [0 -10],len,nArrow);

        % rev A
        %Plot.arrow([-10 0], [-10 0],len,nArrow);
        %Plot.arrow([ 10 0], [ 10 0],len,nArrow);
        %Plot.arrow([ 10 0], [-10 0],len,nArrow);
        %Plot.arrow([-10 0], [ 10 0],len,nArrow);

        xl=[-11 11];
        yl=[-11 11];
        xlim(xl);
        ylim(yl);
        axis square;


    end
    function [hline]=arrow(x,y,varargin)
    % XXX NOTE NOT FINISHED see arrowTest_()
    %function h=arrow(x,y,len,varargin)
    %function h=arrow(x,y,varargin)
        %x1 y1 xd yd

        nvargs=numel(varargin);
        nArrow=1;

        %% parse numeric input
        bNums=nvargs > 0 && isnumeric(varargin{1});
        if bNums
            nums=varargin{1};
            varargin=varargin(2:end);
            nvargs=nvargs-1;

            % radius
            r=nums(1);

            % angle
            if nvargs > 0 && isnumeric(varargin{1}) && numel(nums)==1
                ang=varargin{1};
                varargin=varargin(2:end);
                nvargs=nvargs-1;
            elseif numel(nums)==1
                ang=30;
            elseif numel(nums) > 1
                ang=nums(2);
                nums=nums(2:end);
            end

            % bBoth
            if nvargs > 0 && isnumeric(varargin{1})
                nArrow=varargin{1};
                varargin=varargin(2:end);
                nvargs=nvargs-1;

                bAlt=nArrow < 0;
                nArrow=abs(nArrow);
            end
        end

        %% Get Spec
        if nvargs > 0 && mod(nvargs,2)~=0
            spec=varargin(1);
            varargin=varargin(2:end);
            nvargs=nvargs-1;
        else
            spec={};
        end

        %% Get Pairs
        P_quiver={...
           'AutoScale','on';
           'AutoScaleFactor',1;
           'MaxHeadSize',0.2;
           'ShowArrowHead','on';
        };
        [vargsQuiv,vargsLine]=Args.simpleLoose('args',P_quiver,varargin{:});


        %% Plot no len
        if ~bNums
            % startx, starty, lenx, leny, scale
            hq=quiver(x(1),y(1),x(end)-x(1),y(end)-y(1),0,spec{:},vargsQuiv{:});
            if nargout >= 1
                hline=hq;
                hquiv=hq;
            end
            return
        end


        %% Plot line
        if ~ishold();
            cl=onCleanup(@() hold('off'));
        end
        hl=[];

        hl=plot(x,y,spec{:},vargsLine{:}); hold on;

        %% GET ARROW DIMENSIONS
        X=x(end)-x(end-1);
        Y=y(end)-y(end-1);
        R=sqrt(X.^2+Y.^2);

        Ang=acosd(X/R);
        if Y < 0
            Ang=Ang+180;
            C=-1;
        else
            C=1;
        end

        len(1)=cosd(Ang)*r*C;
        len(2)=sqrt(r^2 - len(1)^2)*C;
        if ang==90
            h=r;
        else
            h=tand(ang)*r;
        end
        yh=sind(90-Ang)*h;
        xh=sqrt(h^2-yh^2);


        %% HANDLE nArrow
        inds=round(linspace(1,numel(x),nArrow)); % TODO INTERP
        Xs=x(inds);
        Ys=y(inds);
        if ang~=90
            X0s=Xs-len(1);
            Y0s=Ys-len(2);
        else
            X0s=Xs;
            Y0s=Ys;
        end

        %% PLOT ARROWS
        for i = 1:nArrow
            xa(1)=X0s(i)-xh;
            xa(2)=X0s(i)+xh;
            ya(1)=Y0s(i)+yh;
            ya(2)=Y0s(i)-yh;

            args={'Color',hl.Color,'LineWidth',hl.LineWidth};
            plot([xa(1) Xs(i)],[ya(1) Ys(i)],args{:});
            plot([xa(2) Xs(i)],[ya(2) Ys(i)],args{:});
        end

        if nargout > 0
            hline=hl;
        end

    end
    function [counts,edges,h]=histline(data,edges,varargin)
        % histc for counts only
        if nargin < 2
            edges=[];
        end
        if isempty(edges)
            [counts,edges]=hist(data);
        else
            [counts,edges]=hist(data,edges);
        end
        hh=plot(X(:),Y(:),varargin{:});
        if nargout > 0
            h=hh;
        end
    end
    function [counts,ctrs,h]=hist(data,edges,varargin)
        if nargin < 2
            edges=[];
        end
        bEdges=isempty(edges);
        if bEdges
            [~,ctrs]=hist(data);
            edges=Plot.ctrs2edges(ctrs,'linear');
        end
        [counts,edges,h]=Plot.bar(data,edges,varargin{:});
        if ~bEdges
            ctrs=Plot.edges2ctrs(edges,'linear');
        end
    end
    function [counts,edges,h]=histogram(data,edges,varargin)
        if nargin < 2
            edges=[];
        end
        if numel(varargin) > 1 && isnumeric(varargin{1})
            counts=varargin{1};
        end
        % edges NOT centers
        if isempty(edges)
            [counts,edges]=histcounts(data);
        else
            [counts,edges]=histcounts(data,edges);
        end
        h=Plot.bar([],edges,counts,varargin{:});
    end
    function [counts,edges,h]=bar(data,edges,counts,varargin)
        % LEAVE DATA EMPTY IF HAVE COUNTS

        if ~isempty(data) && ~isempty(counts)
            varargin=[counts varargin];
        end
        bCounts=~isempty(counts) && (isempty(data) || isnumeric(counts));
        if nargin < 2
            edges=[];
        end
        if bCounts
            counts=Vec.row(counts);
        elseif isempty(edges)
            [counts,edges]=hist(data);
        else
            [counts,edges]=hist(data,edges);
        end
        d=diff(edges);
        dh=d(1)/2;
        Y=[counts; counts];
        X=[edges(1:end-1); edges(2:end)];
        X=X(:);
        Y=Y(:);

        P={'FaceAlpha',0.2;
           ...
           'Color',[];
           'EdgeColor',[]; % pref
            ...
           'FaceColor',[];
           'FillColor',[]; % pref
           'LineWidth',2;
           'bDiv',[];
        };
        [S,varargin]=Args.simpleLoose([],P,varargin{:});

        % Pref
        if ~isempty(S.Color) && isempty(S.EdgeColor)
            S.EdgeColor=S.Color;
        end
        if ~isempty(S.FaceColor) && isempty(S.FillColor);
            S.FillColor=S.FillColor;
        end

        if isempty(S.EdgeColor)
            color={};
        else
            color={'Color',S.EdgeColor};
        end

        if ~ishold()
            cl=onCleanup(@() hold('off'));
        end


        bLine=S.LineWidth>0 || (~isempty(S.bDiv) && S.bDiv);
        if isempty(S.bDiv)
            S.bDiv=bLine;
        end
        if S.bDiv
            XL=repelem(X,3,1);
            YL=repelem(Y,3,1);
            YL(1:3:end)=0;
            YL(3:3:end)=0;
        else
            XL=X;
            YL=Y;
        end

        if ~bLine
            S.LineWidth=1;
            if isempty(S.FillColor) && isempty(S.EdgeColor)
                S.FillColor = nextColor;
            elseif isempty(S.FillColor)
                S.FillColor=S.EdgeColor;
            end
            S.EdgeColor='none';
        else
            if isempty(S.FillColor) && isempty(S.EdgeColor)
                S.FillColor = nextColor;
            elseif isempty(S.FillColor)
                S.FillColor=S.EdgeColor;
            end
            %S.EdgeColor='none';
        end
        if bLine
            h1=plot(XL,YL,'Color',S.EdgeColor,'LineWidth',S.LineWidth,varargin{:}); hold on
        else
            h1=[];
        end

        if isempty(S.FillColor)
            S.FillColor=lastColor;
        end
        if Y(1)~=0
            Y=[0; Y];
            X=[X(1); X];
        end
        if Y(end)~=0
            Y=[Y; 0];
            X=[X; X(end)];
        end
        if S.bDiv
            h1=plot(X,Y,'Color',S.EdgeColor,'LineWidth',S.LineWidth,varargin{:}); hold on
        end
        h2=fill(X,Y,S.FillColor,'EdgeColor','none');
        h2.FaceAlpha=S.FaceAlpha;

        if nargout > 0
            h=[h1 h2];
        end
    end
    function h=poles(x,y,varargin)
        if nargin >= 3 && ~isnumeric(y)
            varargin=[{y} varargin];
            y=x;
            x=1:numel(y);
        elseif (nargin < 2 || isempty(y))
            y=x;
            x=1:numel(y);
        end
        [color,varargin,bSuccess]=Args.getPair('color',varargin{:});
        if ~bSuccess && length(varargin{1})==1 && ischar(varargin{1})
            color=varargin{1};
            varargin=varargin(2:end);
            ca={'Color',color,'MarkerFaceColor',color};
            cc={'Color',color};
        elseif ~bSuccess
            ca={};
            cc={};
        else
            ca={'Color',color,'MarkerFaceColor',color};
            cc={'Color',color};
        end

        X=Vec.row(x);
        Y=Vec.row(y);

        val=std(Y(~isinf(Y)))*5;
        ind=(isinf(Y) | Y > val);

        Y(ind)=val;
        % PLOT Markers
        a={'Marker','o','LineStyle','none'};
        args=[a ca varargin] ;
        plot(X(~ind),Y(~ind),args{:});
        hold on;

        a={'Marker','x','LineStyle','none'};
        args=[a ca varargin] ;
        plot(X(ind),Y(ind),args{:});

        % PLOT LINES
        n=nan(size(x));
        x=[X; X; X; n];

        z=zeros(1,numel(y));
        n=nan(size(z));
        y=[z; Y; z; n];

        args=[varargin cc];
        out=plot(x(:),y(:),args{:});
        if nargout > 0; h=out; end;
    end
    function fill(x,y,varargin)
        vargs=Plot.parseArgs('line',varargin{:});
        fill(x,y,'r',vargs{:});
    end
    function h=plot(x,y,varargin)
        vargs=Plot.parseArgs('line',varargin{:});
        out=plot(x,y,vargs{:});
        if nargout > 0; h=out; end;
    end
    function h=RC(RC,varargin)
        vargs=Plot.parseArgs('line',varargin{:});
        out=plot(RC(:,2),RC(:,1),vargs{:});
        if nargout > 0; h=out; end;
    end
    function h= XY(RC,varargin)
        vargs=Plot.parseArgs('line',varargin{:});
        out=plot(RC(:,1),RC(:,2),vargs{:});
        if nargout > 0; h=out; end;
    end
    function h = RC3(xyzRC,varargin)
        %function plot3sane(xyzRC,varargin)

        out=plot3(xyzRC(:,1),xyzRC(:,3),xyzRC(:,2),varargin{:});
        if nargout > 0; h=out; end;
        ylabel('Z');
        zlabel('Y');
        xlabel('X');
    end
    function h = XYZ(xyzRC,varargin)
        %function plot3sane(xyzRC,varargin)

        out=plot3(xyzRC(:,:,1),xyzRC(:,:,3),xyzRC(:,:,2),varargin{:});
        if nargout > 0; h=out; end;
        ylabel('Z');
        zlabel('Y');
        xlabel('X');
    end
    function h = imagescRC3(xyzRC,xMinMax,yMinMax,dnk,varargin)
        if nargin < 4
            dnk=[];
            if nargin < 3
                yMinMax=[];
                if nargin < 2
                    xMinMax=[];
                end
            end
        end

        [xq,yq,zq]=Plot.scatter_interp_fun(xyzRC,xMinMax,yMinMax,dnk);

        imagesc(xq(:),yq(:),zq);
        set(gca,'YDir','normal');
        Fig.formatIm();

        %ylabel('Z');
        %zlabel('Y');
        %xlabel('X');
    end
    function h = surfRC3(xyzRC,xMinMax,yMinMax,dnk,varargin)

        [xq,yq,zq]=scatter_interp_fun(xyzRC,xMinMax,yMinMax,dnk);

        f=surf(xq,zq,yq);
        f.EdgeColor='none';

        ylabel('Z');
        zlabel('Y');
        xlabel('X');

        %axis vis3d;
        %axis off;
        %l = light('Position',[-50 -15 29]);
        %set(gca,'CameraPosition',[208 -50 7687]);
        %lighting phong;
        %shading interp;
        %colorbar EastOutside;

        if nargout >= 1
            h=f;
        end
    end

    function error(mx,my,dx,dy,line,varargin)
        if size(mx,1)==1
            mx=transpose(mx);
        end
        if size(my,1)==1
            my=transpose(my);
        end
        if size(dx,1)==1
            dx=transpose(dx);
        end
        if size(dy,1)==1
            dy=transpose(dy);
        end
        if size(dx,2)==1
            dx=[dx dx];
        end
        if size(dy,2)==1
            dy=[dy dy];
        end
        if ~exist('line','var') || isempty(line)
            line='k-';
        end

        ly=my;
        ry=my;
        bx=mx;
        tx=mx;

        bAdd=1;
        if bAdd
            lx=mx-dx(:,1);
            rx=mx+dx(:,2);
            by=my-dy(:,1);
            ty=my+dy(:,2);
        else
            lx=dx(:,1);
            rx=dx(:,2);
            by=dy(:,1);
            ty=dy(:,2);
        end

        for i = 1:size(bx,1)
            plot([bx(i),tx(i)],[by(i),ty(i)],line,varargin{:}); hold on;
            plot([lx(i),rx(i)],[ly(i),ry(i)],line,varargin{:});
        end

    end
    function intervStd(X,m,std,C,varargin)

        if ~exist('C','var')
            C=[];
        end
        U=m+std;
        L=m-std;

        Plot.interv(X,U,L,C,varargin{:});
        hold on
        plot(X,m,'ko','MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
        hold off
    end
    function h=interv(X,U,L,varargin)

        bNum1=nargin >= 4 && isnumeric(varargin{1});
        bNum2=nargin >= 5 && isnumeric(varargin{2});
        bC=bNum1;
        bZ=bNum1 && bNum2;
        if Args.hasPair('Color',varargin{:});
            [C,varargin]=Args.getPair('Color',varargin{:});
        elseif bC
            C=varargin{1};
            varargin(1)=[];
        else
            C=[.5 .5 .5];
        end
        if bZ
            error('TODO')
        end
        X=Vec.row(X);
        X=[X fliplr(X)];
        U=Vec.row(U);
        L=Vec.row(L);
        Y=[U fliplr(L)];
        h=patch(X,Y,C,varargin{:});
    end
    function h=interRC(X,LU,varargin)
        h=Plot.interv(X,L(:,1),U(:,2),varargin{:});
    end

    function []=rect(ImCtrRC,h,w,color,LineWidth,bPlotCtr)
    %function []=Plot.rect(ImCtrRC,h,w)
    % example call:
    %    ImCtrRC=[500 1000; 1000 1000];
    %    h=100;
    %    w=100;
    %    pht=dbImg.getImg('LRSI','img','phtGamma',8,'L');
    %    IszRC=size(pht);
    %    imagesc(pht);hold on
    %    Plot.rect(ImCtrRC,h,w)
    %    Fig.formatIm();
        if ~exist('color','var') || isempty(color)
            color='g';
        end
        if ~exist('LineWidth','var') || isempty(LineWidth)
            LineWidth=1;
        end
        if ~exist('bPlotCtr','var') || isempty(bPlotCtr)
            bPlotCtr=0;
        end



        [x,y]=Rec.rect(ImCtrRC,h,w);
        x=x+0.5;
        y=y+0.5;
        x(:,5)=x(:,1);
        y(:,5)=y(:,1);
        hold on;
        for i = 1:size(x,1)
            plot(x(i,:),y(i,:),'color',color,'LineWidth',LineWidth);
        end
        if bPlotCtr;
            p=ImCtrRC+0.5;
            plot(p(2),p(1),'o','color',color,'MarkerFaceColor',color);
        end
        hold off

    end
    function hiddenData(varargin);
        [x,y]=Axis.getHiddenData();
        plot(x,y,varargin{:});
    end
    function opts=parse(opts,varargin)
        if nargin < 1 || isempty(opts)
            opts=struct();
        end
        bP=nargin >= 2 && ~isempty(varargin{1});
        if nargin ==2 && isstruct(varargin)
            P=varargin{1};
        else
            P=struct(varargin{:});
        end
        FLDS={'FillColor';
              'FaceAlpha';
              'EdgeColor';
              ...
              'LineStyle';
              'LineWidth';
              'LineColor';
              ...
              'Marker';
              'MarkerSize';
              'MarkerFaceColor';
              'MarkerEdgeColor';
              ...
              'FontSize';
              'CI';
        };
        flds=fieldnames(opts);
        for i = 1:length(FLDS)
            f=FLDS{i};
            bInOpts=any(strcmp(FLDS{i},flds));
            bInP=bP && isfield(P,f);
            if bInP
                if ~bInOpts || isempty(opts.(f));
                    opts.(f)=P.(f);
                end
            elseif ~bInOpts
                opts.(f)=[];
            end
        end
    end
    function vargs=parseArgs(type,varargin)
        if nargin < 1 || isempty(type)
            type='line';
        end
        switch type
            case 'line'
                vargs=Plot.parseOutShortHandLine(varargin{:});
                vargs=Plot.setDefaultsLine(vargs{:});
        end
    end
    function [out,vargs,sh]=parseOutShortHandLine(varargin);
        vargs=varargin;
        sh=[];
        if isempty(vargs)
            out=vargs;
            return
        end
        ind=cellfun(@Plot.isMatchShortHandLine,varargin);
        ind=cumprod(ind);
        if ~any(ind)
            out=vargs;
            return
        end
        sh=[vargs{ind}];
        vargs=vargs(~ind);
        if isempty(sh)
            out=vargs;
        else
            out=[Plot.parseShortHandLine(sh) vargs];
        end
    end
end
methods(Static,Access=private)
    function [xq,yq,zq] = scatter_interp_fun(xyzRC,xMinMax,yMinMax,dnk);
    %function handle = Plot.surfRC3(xyzRC,xMinMax,yMinMax,varargin)
        if isempty(dnk)
            dnk=1;
        end

        indEmpty=isnan(xyzRC(:,3)) & ~any(isnan(xyzRC(:,1:2)),2);
        xyEmpty=xyzRC(indEmpty,1:2);

        indNan=any(isnan(xyzRC),2);
        xyzRC=xyzRC(~indNan,:);


        x=xyzRC(:,1);
        y=xyzRC(:,2);
        z=xyzRC(:,3);

        %%
        % remove outliers before MM
        if isempty(xMinMax)
            xr=rmoutliers(x,'grubbs');
            xMinMax=Num.minMax(xr);
        end
        if isempty(yMinMax)
            yr=rmoutliers(y,'grubbs');
            yMinMax=Num.minMax(yr);
        end

        xRng=abs(diff(xMinMax));
        yRng=abs(diff(yMinMax));

        d=xRng/yRng;

        xRes=xRng/(sqrt(numel(x))*d)*dnk;
        yRes=yRng/(sqrt(numel(y))/d)*dnk;

        xlin=xMinMax(1):xRes:xMinMax(2);
        ylin=yMinMax(1):yRes:yMinMax(2);
        [xq,yq] = meshgrid(xlin, ylin);
        F = scatteredInterpolant(x,y,z,'natural','none');
        zq=F(xq,yq);

        indx=knnsearch(transpose(xlin),xyEmpty(:,1));
        indy=knnsearch(transpose(ylin),xyEmpty(:,2));

        zq(indx,indy)=nan;
    end
end
methods(Static,Hidden)
    function D=lineDefaults()
        D=struct();
        D.Color=[0 0 0];
        D.Marker='.';
        D.MarkerEdgeColor='auto';
        D.MarkerFaceColor='auto';
        D.MarkerSize =10;
        D.LineStyle='-';
        D.LineWidth=2;
    end
    function out=setDefaultsLine(varargin)
        out=varargin;
        S=struct(varargin{:});
        D=Plot.lineDefaults();
        flds=fieldnames(D);
        flds(ismember(flds,{'Marker','LineStyle'}))=[];
        for i = 1:length(flds)
            fld=flds{i};
            if ~isfield(S,fld)
                out=[fld,D.(fld),out];
            end
        end
        if ~isfield(S,'LineStyle')
            if isfield(S,'Marker')
                out=['LineStyle', 'none', out];
                out=['Marker', D.Marker, out];
            else
                out=['LineStyle', D.LineStyle, out];
                out=['Marker', 'none', out];
            end
        end
    end
    function out=parseShortHandLine(str)
        [vals,type]=Plot.readShortHandLine(str);
        mind=ismember(type,'m');
        lind=ismember(type,'l');
        cind=ismember(type,'c');
        xind=ismember(type,'x');

        mErr=sum(mind) > 1;
        lErr=sum(lind) > 1;
        cErr=sum(cind) > 1;
        xErr=sum(xind) > 0;


        mCount=0;
        if mErr
            mErrStr=sprintf('  Too many markers definitions: %s \n',[vals{mind}]);
        else
            mErrStr='';
        end
        if lErr
            lErrStr=sprintf('  Too many lines definitions: %s \n',[vals{lind}])
        else
            lErrStr='';
        end
        if cErr
            cErrStr=sprintf('  Too many color definitions: %s \n',[vals{cind}])
        else
            cErrStr='';
        end
        if xErr
            xErrStr=sprintf('  Unknown characters: %s \n',[vals{xind}]);
        else
            xErrStr='';
        end

        ErrStr=[mErrStr lErrStr cErrStr xErrStr];
        if ~isempty(ErrStr)
            error(ErrStr)
        end

        out={};
        if any(mind)
            out{end+1}='Marker';
            out{end+1}=vals{mind};
        end
        if any(lind)
            out{end+1}='LineStyle';
            out{end+1}=vals{lind};
        end
        if any(cind)
            out{end+1}='Color';
            out{end+1}=vals{cind};
        end

    end
    function [out,type]=readShortHandLine(str)
        n=length(str);
        out={};
        type='';
        bSkip=false;
        for i = 1:n
            if bSkip
                bSkip=false;
                continue
            end

            if i+1 <= n
                match=Str.RE.match(str(i:i+1),[Plot.linesRE],0);
            else
                match=Str.RE.match(str(i),['^' Plot.linesRE],0);
            end
            if ~isempty(match)
                out{end+1}=match;
                type(end+1)='l';
                if length(match) > 1
                    bSkip=true;
                end
                continue
            end
            match=Str.RE.match(str(i),Plot.colorsRE);
            if ~isempty(match)
                out{end+1}=match;
                type(end+1)='c';
                continue
            end
            match=Str.RE.match(str(i),Plot.markersRE);
            if ~isempty(match)
                out{end+1}=match;
                type(end+1)='m';
                continue
            end
            out{end+1}=str(i);
            type(end+1)='x';
        end
    end
    function out=shortHandLineRE()
        out=sprintf('^(%s|%s|%s){1,3}$', Plot.colorsRE, Plot.markersRE, Plot.linesRE);
    end
    function out=isMatchShortHandLine(str)
        if ~ischar(str)
            out=false;
            return
        end
        out=Str.RE.ismatch(str,Plot.shortHandLineRE);
    end
    function edges=ctrs2edges(ctrs,extMethod)
        %if ~isvec(ctrs)
        %    error('write code')
        %end
        if ~exist('extMethod','var')
            extMethod='linear';
        end

        edges = conv(ctrs, [0.5, 0.5], 'valid');
        widths=diff(edges);

        if numel(unique(widths))==1
            edges=[edges(1)-widths(1) edges edges(2)+widths(1)];
        end
    end
    function [ctrs,widths]=edges2ctrs(edges)
    % function [ctrs,widths]=edges2ctr(edges)
    % also see
    %histogram('BinEdges',edges,'BinCounts',counts)
        if ~Vec.is(edges)
            error('write code')
        end
        lgwidths=diff(log10(edges(edges ~= 0)));
        if abs(diff(lgwidths)) <= 1e-12
            width=mean(lgwidths);
            ctrs=log10(edges(1:end-1))+width/2;
            ctrs=10.^(ctrs);
        else
            widths=diff(edges);
            ctrs=edges(1:end-1)+widths/2;
        end
    end
    function [ctrs,w,N]=binWidthsFD(data,bRmOutliers)
        % Freedman-Diaconis rule
        % binwidths = 2 * IQR(x)/cuberoot(n)
        % TODO
        % logscale
        %
        % online -> 1. online var
        %           2. var -> iqr

        if ~exist('bRmOutliers','var') || isempty(bRmOutliers)
            bRmOutliers=0;
        end
        if bRmOutliers
            data=rmoutliers(data(:),'quartiles');
        else
            data=data(:);
        end
        data(isinf(data) | isnan(data))=[];
        w=2.*iqr(data)./nthroot(numel(data),3);
        N=(max(data)-min(data))./w;
        ctrs=linspace(min(data),max(data),N);
    end
    function out=scale2hist(xp,yp,xh,yh)
        % XXX not quite
        % dT=totT/N  bin width
        % area=sum(yh)
        dh=numel(xh)./Num.range(xh);
        dp=numel(xp)./Num.range(xp);
        out=yp.*dp./dh.*sum(yh);
    end


end
end
