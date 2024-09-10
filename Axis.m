classdef Axis < handle
methods(Static)
    function limSquare(xl,yl,p);
        if nargin < 3
            p=[];
        end
        if nargin < 1
            xl=xlim();
        end
        if nargin < 2
            yl=ylim();
        end
        xl=Num.minMax([xl yl]);
        if nargout < 1
            Axis.ylim(xl,p);
            Axis.xlim(xl,p);
        end
    end
    function xl=xlim(xl,p)
        if nargin < 1 || isempty(xl)
            xl=xlim;
        end
        if nargin < 2
            p=0;
        end
        xl=Axis.ylim(xl,p);
        if nargout < 1
            xlim(xl);
        end
    end
    function yl=ylim(yl,p)
        if nargin < 1 || isempty(yl)
            yl=ylim;
        end
        if nargin < 2 || isempty(p)
            p=0.05;
        end
        if numel(yl)==1
            yl=sort([-yl yl]);
        elseif numel(yl)>2
            yl=[min(yl(:)) max(yl(:))];
        end
        r=abs(diff(yl));
        m=r*p;
        ind=yl >= 0;
        ind(ind==0)=-1;
        yl=yl+[ind].*[-1 1]*m;
        if nargout < 1
            ylim(yl);
        end
    end
    function fit(a,f)
        if nargin < 1 || isempty(a)
            a=gca;
        end
        if nargin < 2 || isempty(f)
            f=gcf;
        end
        f.Position;
        set(a,'Position',[0 0 1.0 1]);

    end
    function ax=set_full_text_margin(l,b,r,t)
        %XLabel.Postition
        %
        %ax.FontSizeMode='manual';
        %
        %f=figure(1)
        %f.DockControls='off'
        %f.HandleVisibility='off'
        %f.ToolBar='none'

        f=figure(1)
        plot(1,1)
        %xlabel('abc')
        %ylabel('abc')
        ax=gca;

        ax.DataAspectRatioMode='manual';
        ax.PlotBoxAspectRatioMode='manual';
        ax.DataAspectRatio=[1 1 1];
        ax.PlotBoxAspectRatio=[1,1,1];
        ax.CameraViewAngleMode='manual';

        %axis tight;
        %ax.ActivePositionProperty='position';
        %ax.ActivePositionProperty='outerposition'
        %set(gca, 'LooseInset', [0,0,0,0]);
        %set(0,'DefaultAxesLooseInset',[0,0,0,0])
        %
        % innerposition and outerposition are linked

        %set(ax,'WarpToFillMode','manual')
        set(ax,'OuterPositionMode','manual')
        set(ax,'ActivePositionPropertyMode','manual')
        set(ax,'ActivePositionProperty','position')
        %Extent + Margin = text box size
        %ax.XLabel.Margin=30 XXX ??


        oNorm=[0,0,1,1];
        iNorm=[[l,b] [1,1]-[r,t]];

        %set(ax,'Units','points');
        %oPnt=get(ax,'OuterPosition');

        %WPnt=[oPnt(3)-oPnt(1)];
        %HPnt=[oPnt(4)-oPnt(2)];

        %LBRT=[WPnt*l HPnt*t WPnt*r HPnt*b];
        %I=[oPnt(1:2)+LBRT(1:2) oPnt(3:4)-LBRT(3:4)];

        %ax.InnerPosition
        %ax.OuterPosition




        pos=get(ax,'Position');
        w=pos(4)-pos(2);
        h=pos(3)-pos(1);

        set(ax,'LooseInset',get(ax,'TightInset'));
        %[ax.InnerPosition(1:2)-ax.OuterPosition(1:2) ax.OuterPosition(3:4)-ax.InnerPosition(3:4)]
        %set(ax,'Units','normalized');
        %set(gca, 'Position', get(gca, 'OuterPosition') - ...
        %    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        %set(ax,'Position',iNorm,'OuterPosition_I',oNorm,'OuterPosition',oNorm)
        %set(ax,'OuterPosition',oNorm)
        %set(ax,'Position',iNorm)
        %set(ax,'Position',iNorm);
        %set(ax,'LooseInset',get(ax,'TightInset'))
        %ax.PlotBoxAspectRatio=[1 1 1];
        %ax.InnerPosition
        %ax.OuterPosition
        %ax.OuterPosition(4)-ax.OuterPosition(2)
        %ax.OuterPosition(3)-ax.OuterPosition(1)
        %ax.TightInset
        %ax.LooseInset
        ax.OuterPosition
        ax.InnerPosition

        % Outerposition
        % Innerposition
        % looseinset
        % tight inset

        %FontUnits
        %FontSize
        %TitleFontSizeMultiplier
        %LabelFontSizeMultiplier
    end
    function match_res(mult)
        if ~exist('mult','var') || isempty(mult)
            mult=1;
        end
        set(gcf,'Units','Pixels');
        set(gca,'Units','Pixels');
        a=get(gca,'Children');

        width=size(a.CData,2);
        height=size(a.CData,1);

        posF=get(gcf,'Position');
        posI=get(gca,'InnerPosition');
        posO=get(gca,'OuterPosition');


        %ratioIO=posI(3:4)./posO(3:4);
        %ratioFO=posF(3:4)./posO(3:4);

        set(gcf,'Position',[posF(1) posF(2) width*mult height*mult]);
        set(gca,'OuterPosition',[0,0,width*mult, height*mult]);
    end
    function set_border_prcnt(prcnt)
        if ~exist('prcnt','var') || isempty(prcnt)
            prcnt=[0.050 0.050];
        end
        if numel(prcnt)==1
            prcnt=[prcnt prcnt];
        end
        LT=prcnt;
        RB=1-((prcnt-0.005).*2);
        units=get(gca,'Units');
        set(gca,'Units','Normalized');
        %posa=get(gca,'Position');
        set(gca,'Position',[LT(1), LT(2), RB(1), RB(2)]);
        set(gca,'Units',units);
    end
    function set_border_pix(X,Y)
        % TODO
        if ~exist('Y','var') || isempty(Y)
            Y=X;
        end
        units=get(gca,'Units');
        Axis.set_border_prcnt(0);
        set(gca,'Units','Pixels');
        posApix=get(gca,'Position');
        prcnt=[X Y]./(posApix(3:4)+[X Y]);
        Axis.set_border_prcnt(prcnt);
    end
    function lims=padLims(limsORdata,prcnt)
        if ~exist('prcnt','var') || isempty(prcnt)
            prcnt=0.15;
        end

        if numel(limsORdata) > 1
            data=limsORdata;
            lims=Num.minMax(data);
        else
            lims=limsORdata;
        end

        rang=lims(2)-lims(1);
        t=rang*prcnt;
        lims=lims+[-t t];
    end
%% GET LIMS
    function lims=getLims(dataORlims,padPrcnt)
        if nargin < 2
            padPrcnt=[];
        end
        if ~exist('dataORlims') || isempty(dataORlims)
            a=gca;
            dataORlims=a.Children.YData;
        end
        lims=Axis.padLims(dataORlims,padPrcnt);
    end
    function [lims,limsDataMinMax]=getLimsPrcntile(prcnt)
        if ~exist('prcnt','var') || isempty(prcnt)
            prcnt=5;
        end
        a=gca;
        p=prcnt/2;
        Y=a.Children.YData;
        lims=prctile(Y,[0+p, 100-p]);
        if nargout < 2
            return
        end
        ind=Y >= lims(1) & Y <= lims(2);
        limsDataMinMax=Num.minMax(Y(ind));
    end
    function [lims,count,ind]=getLimsOutlier(method)
    % median, gesd
        if ~exist('method','var') || isempty(method)
            method='median';
        end
        a=gca;
        Y=a.Children.YData;
        %ind=~isoutlier(Y,'gesd');
        ind=isoutlier(Y,method);
        count=sum(ind);
        lims=Num.minMax(Y(~ind));
    end
%% SET LIMS
    function varargout=setLims(dataORlims,padPrcnt)
        if ~exist('padPrcnt','var')
            padPrcnt=[];
        end
        if ~exist('dataORlims')
            dataORlims=[];
        end
        lims=Axis.getLims(dataORlims,padPrcnt);

    end
    function varargout=setLimsDiag()
        lims=Num.minMax([xlim(); ylim()]);
        xlim(lims);
        ylim(lims);
        if nargout > 0
            vargout{1}=lims;
        end
    end
    function varargout=setLimsPrcntile(prcnt,padPrcnt);
        if ~exist('prcnt','var')
            prcnt=[];
        end
        if ~exist('padPrcnt','var')
            padPrcnt=[];
        end
        lims=Axis.getLimsPrcntile(prcnt);
        lims=Axis.padLims(lims,padPrcnt);
        if nargout > 0
            varargout{1}=lims;
        end
        ylim(lims);
    end
    function varargout=setLimsDataPrcntile(prcnt,padPrcnt);
        if ~exist('prcnt','var')
            prcnt=[];
        end
        if ~exist('padPrcnt','var')
            padPrcnt=[];
        end
        [~,lims]=Axis.getLimsPrcntile(prcnt);
        lims=Axis.padLims(lims,padPrcnt);
        if nargout > 0
            varargout{1}=lims;
        end
        ylim(lims);
    end
    function varargout=setLimsOutlier(method,padPrcnt)
        if ~exist('method','var')
            method=[];
        end
        if ~exist('padPrcnt','var')
            padPrcnt=[];
        end
        [lims,count,ind]=Axis.getLimsOutlier(method);
        lims=Axis.padLims(lims,padPrcnt);
        if nargout > 0
            varargout{1}=lims;
        elseif nargout > 1
            varargout{2}=count;
        elseif nargout > 2
            varargout{3}=ind;
        end
        ylim(lims);
    end
%% OUTLIERS
    function [x,y]=getHiddenData()
        ylims=ylim();
        a=gca;
        Y=a.Children.YData;
        negInd=Y < ylims(1);
        posInd=Y > ylims(2);

        X=a.Children.XData;
        x=[X(negInd) X(posInd)];
        y=[repmat(ylims(1),1,sum(negInd)) repmat(ylims(2),1,sum(posInd))];

    end
%% RESOLUTION
    function []=setRes(mult)
        if ~exist('mult','var') || isempty(mult)
            mult=1;
        end
        set(gcf,'Units','Pixels');
        set(gca,'Units','Pixels');
        a=get(gca,'Children');
        width=size(a.CData,2);
        height=size(a.CData,1);

        posf=get(gcf,'Position');
        set(gcf,'Position',[posf(1) posf(2) width*mult height*mult]);
        posI=get(gca,'InnerPosition');
        set(gca,'OuterPosition',[0,0,width*mult, height*mult]);
        set(gca,'Units','Normalized');
        posa=get(gca,'Position');
        set(gca,'Position',[0.055, posa(2),0.9, posa(4)]);
    end
    function []=clabel(h,text)
        h.Label.String = text;
        pos = get(h,'Position');
        h.Label.Position = [pos(1)/2 pos(2)-1]; % to change its position
        h.Label.Rotation = 0; % to rotate the text
    end
    function out=isInvalid(h)
        out=isempty(h) || ~isvalid(h) || ~isgraphics(h);
    end
    function writeText(xPositions,yPositions,txtStrings,ratioORabs,fs,hPos,rotDeg,vPos,textColor,zPositions)

        % function Axis.writeText(xPositions,yPositions,txtStrings,ratioORabs,fs,hPos,rotDeg,vPos,textColor,zPositions)
        %
        %   example call: Axis.writeText(.1,.9,{['R^{2} = .5']})
        %
        % writes txtStrings at specified location in figure window
        %
        % xPositions: scalar between 0 and 1 that determines position text will
        %             appear in x (1xn)
        % yPositions: scalar between 0 and 1 that determines position text will
        %             appear in y (1xn)
        % txtStrings: cell array of strings that get written at location
        %             (xPositions,yPositions,zPositions)
        % ratioORabs: indicates whether xPositions & yPositions indicate text
        %             position in the current window as percentage of window size or in
        %             absolute x or y positions
        % fs:         fontsize
        % hPos:      'left','center', 'right' alignment
        % rotDeg:     orientation of text in degrees
        % vPos:       'top', 'bottom', 'middle'
        % textColor:  default: black, [0 0 0]
        % zPositions:

        if (length(xPositions) ~= length(yPositions) | length(xPositions) ~= length(txtStrings))
            error('writeText: all three variables [xPositions,yPositions,txtStrings] must have same number of elements');
        end
        if (~iscell(txtStrings))
            error('writeText: txtStrings must be of type cell');
        end
        if (~exist('ratioORabs','var') | isempty(ratioORabs))
            ratioORabs = 'ratio';
        end
        if (~exist('fs','var') || isempty(fs))
            fs = 18;
        end
        if (~exist('hPos','var'))
            hPos = 'left';
        end
        if (~exist('vPos','var'))
            vPos = 'middle';
        end
        if (~exist('rotDeg','var') || isempty(rotDeg))
            rotDeg = 0;
        end
        if ~exist('textColor','var') || isempty(textColor)
            textColor = [ 0 0 0];
        end
        if ~exist('zPositions','var') || isempty(zPositions)
            if strcmp(ratioORabs,'ratio')
            zPositions = .5*ones(size(xPositions));
            elseif strcmp(ratioORabs,'abs')
            zPositions = 0*ones(size(xPositions));
            end
        end
        % GET AXIS LIMITS
        [xlims ylims] = get_lims;

        if (strcmp(ratioORabs,'abs'))
            xPos = xPositions;
            yPos = yPositions;
            zPos = zPositions;
        elseif (strcmp(ratioORabs,'ratio'))
            xPos = xlims(1) + xPositions.*diff(xlims);
            yPos = ylims(1) + yPositions.*diff(ylims);
            zlims = get(gca,'zlim');
            zPos = zlims(1) + zPositions.*diff(zlims);
        else
            error(['writeText(): invalid ratioORabs value. choose ratio OR abs.']);
        end

        % WRITE TEXT AT INDICATED POSITIONS
        for (s = 1:length(txtStrings))
            text(xPos(s),yPos(s),zPos(s),txtStrings{s},'HorizontalAlignment',hPos,'VerticalAlignment',vPos,'fontsize',fs,'rotation',rotDeg,'color',textColor);
        end
    end


    function [xlims ylims zlims] = get_lims(dim)

        % function [xlims ylims zlims] = Axis.getLims(dim)
        %
        %   x, y, and z axis limits of gca
        %
        % dim: 1 -> returns x lims,
        %      2 -> returns y lims,
        %      3 -> returns z lims
        %      [] or ~exist -> returns xlims, ylims, and zlims

        if nargin < 1
            dim = 1;
        end
        ax=gca;
        x=get(ax,xlim);
        y=get(ax,ylim);
        z=get(ax,zlim);
        if dim == 1
            xlims = [min(x) max(x)];
            ylims = [min(y) max(y)];
            zlims = [min(z) max(z)];
        elseif dim == 2
            xlims = [min(y) max(x)];
        elseif dim == 3
            xlims = [min(z) max(z)];
        else
            error(['Axis.getLims: dim value (' num2str(dim) ') invalid']);
        end
    end
    function [ax,h]=topLeftTitle(txt,fontsize)
        pos=[0, 1.00, 1, 0.05];
        ha='Left';
        va='top';
        tpos=[0.01 0];
        rot=0;

        if nargin < 2
            fontsize=22;
        end
        axes( 'Position', pos);
        set(gca,'Color','None','XColor','None','YColor','None');
        h=text(tpos(1),tpos(2),txt,'FontSize',fontsize,'HorizontalAlignment',ha,'VerticalAlignment',va);
        set(h,'Rotation',rot);
    end
    function [ax,h]=topTitle(txt,fontsize)
        pos=[0, 1.00, 1, 0.05];
        ha='Center';
        va='top';
        tpos=[0.5 0];
        rot=0;

        if nargin < 2
            fontsize=22;
        end
        axes( 'Position', pos);
        set(gca,'Color','None','XColor','None','YColor','None');
        h=text(tpos(1),tpos(2),txt,'FontSize',fontsize,'HorizontalAlignment',ha,'VerticalAlignment',va);
        set(h,'Rotation',rot);
    end
    function [ax,h]=bottomLeftTitle(txt,fontsize)
        ha='Center';
        va='top';
        pos=[0, 0.00, 1, 0.95];
        tpos=[0.5,0];
        rot=0;

        if nargin < 2
            fontsize=22;
        end
        axes( 'Position', pos);
        set(gca,'Color','None','XColor','None','YColor','None');
        h=text(tpos(1),tpos(2),txt,'FontSize',fontsize,'HorizontalAlignment',ha,'VerticalAlignment',va);
        set(h,'Rotation',rot);
    end
    function [ax,h]=bottomTitle(txt,fontsize)
        ha='Center';
        va='top';
        pos=[0, 0.00, 1, 0.95];
        tpos=[0.5,0];
        rot=0;

        if nargin < 2
            fontsize=22;
        end
        axes( 'Position', pos);
        set(gca,'Color','None','XColor','None','YColor','None');
        h=text(tpos(1),tpos(2),txt,'FontSize',fontsize,'HorizontalAlignment',ha,'VerticalAlignment',va);
        set(h,'Rotation',rot);
    end
    function [ax,h]=leftTitle(txt,fontsize)
        ha='center';
        va='top';
        pos=[0.05, 0.05, 1, 0.95];
        tpos=[0,0.5];
        rot=90;

        if nargin < 2
            fontsize=22;
        end
        ax=axes( 'Position',pos);
        set(gca,'Color','None','XColor','None','YColor','None');
        h=text(tpos(1),tpos(2),txt,'FontSize',fontsize,'HorizontalAlignment',ha,'VerticalAlignment',va);
        set(h,'Rotation',rot);
    end
    function test()
        f=gcf;
        subplot()
        a=gca;
        a.Position
        a.XLabel.Position
        a.YLabel.Position

        pos1 = [0.1 0.3 0.3 0.3];
        subplot('Position',pos1)

    end
    function cax=formatIm(g,mp)
        %Defaults for photgraphic images
        if nargin < 2
            mp=[];
        end
        if nargin < 1 || isempty(g)
            g=gca;
        end
        set(g,'xticklabel',{[]});
        set(g,'yticklabel',{[]});
        grid off; axis image;
        set(g,'xtick',[]);
        set(g,'ytick',[]);
        cax=caxis;
    end
    function yticksLog(ax)
        if nargin < 1
            ax=gca;
        end
        YTickLabel = get(ax,'YTick');
        set(ax,'YTickLabel',num2str(YTickLabel'));
    end
    function rmxticklabels(ax)
        if nargin < 1
            ax=gca;
        end
        xt=repmat({''},1,numel(ax.XTick));
        set(ax,'XTickLabel',xt);
    end
    function rmyticklabels(ax)
        if nargin < 1
            ax=gca;
        end
        yt=repmat({''},1,numel(ax.YTick));
        set(ax,'YTickLabel',yt);
    end
    function xticksLog(ax)
        if nargin < 1
            ax=gca;
        end
        XTickLabel = get(ax,'XTick');
        set(ax,'XTickLabel',num2str(XTickLabel'));
    end
    %function formatYLog()
    %end
    function format(varargin)
        %figure(101)
        %clf
        %plot(1,1)
        Axis.format_fun(varargin{:});
        %yticks([0:.2:1])
        %y=ylabel('Proportion Agreement','FontName','Gil Sans')
        %set(y,'FontName','Gil Sans')
        %y
    end
    function format3(varargin)
        if nargin >= 3
            zl=varargin{3};
            varargin(3)=[];
        else
            zl=[];
        end
        Axis.format_fun(varargin{:});
        zlabel(zl);
    end
    function h=format_fun(xlbl,ylbl,ttl,bLogX,bLogY,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD)

        % function Fig.format(xlbl,ylbl,ttl,bLogX,bLogY,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD)
        %
        %   example call: Fig.format('X','Y','Data Plot',1,1,22,18,xtck,ytck,xtckSD,ytckSD))
        %
        % automatically format figure with axis labels, title, log/linear scaling, and specified fontsize
        %
        % xlbl:     x-axis label
        % ylbl:     y-axis label
        % ttl:      title
        % bLogX:    log scale x-axis
        %           0 -> linear scale (default)
        %           1 -> log    scale
        % bLogX:    log scale y-axis
        %           0 -> linear scale (default)
        %           1 -> log    scale
        % fsLbl:   fontsize for labels
        % fsTck:   fontsize for ticks
        % xtck:    x tick values
        % ytck:    y tick values
        % xtckSD:  number of significant digits for xtck values
        % ytckSD:  number of significant digits for ytck values

        if ~exist('ttl','var')     || isempty(ttl),    ttl = [];    end
        if ~exist('h','var')       || isempty(h)       h = gca;     end
        if ~exist('bLogX','var')   || isempty(bLogX)   bLogX = 0;   end
        if ~exist('bLogY','var')   || isempty(bLogY)   bLogY = 0;   end
        if ~exist('fsLbl','var')   || isempty(fsLbl)   fsLbl = 22;  end
        if ~exist('fsTck','var')   || isempty(fsTck)   fsTck = 18;  end

        %set(h,'FontName','Helvetica','FontSize',fsTck,'FontWeight','normal');
        set(h,'FontSize',fsTck,'FontWeight','normal',...
              'LineWidth',2,'TickLength',[0.03, 0.02]);

        try xlabel(xlbl,'fontsize',fsLbl); catch, end
        try ylabel(ylbl,'fontsize',fsLbl); catch, end
        title(ttl,'fontsize',fsLbl);
        if bLogX == 1
            set(h,'xscale','log');
        end
        if bLogY == 1
            set(h,'yscale','log');
        end
        set(h,'fontsize',fsTck,'fontWeight','normal');
        try set(h,'XColor','k'); catch, end
        try set(h,'YColor','k'); catch, end
        %set(gcf,'color','w');
        box on;

        % SET TICKS AND SINGIFICANT DIGITS
        if exist('xtck','var')   & ~isempty(xtck)   set(h,'xtick',[xtck]); end
        if exist('ytck','var')   & ~isempty(ytck)   set(h,'ytick',[ytck]); end
        if exist('xtckSD','var') & ~isempty(xtckSD) set(h,'xticklabel',num2str(get(h,'xtick')',['%.' num2str(xtckSD) 'f'] )); end
        if exist('ytckSD','var') & ~isempty(ytckSD) set(h,'yticklabel',num2str(get(h,'ytick')',['%.' num2str(ytckSD) 'f'] )); end

        % UNBOLD TITLE DEFAULTS FROM MATLABv2015 or later
        v = version('-release');
        if str2num(v(1:4)) >= 2015
            set(h,'TitleFontWeight','normal');
        end
        %axis square;

    end
%% UTIL

end
end
