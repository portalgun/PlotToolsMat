classdef Fig < handle
methods(Static)
%- NEW FIG
    function align(varargin)
        if isnumeric(varargin{1}) && numel(varargin{1}) > 1
            figs=varargin{1};
        end

        set(groot,'CurrentFigure',figs(1)); % No focus XXX
        set(gcf,'Units','pixels');
        f=gcf;
        oPos=get(f,'OuterPosition');
        iPos=get(f,'InnerPosition');

        for i = 1:length(figs)
            set(groot,'CurrentFigure',figs(i));

            f=gcf;
            set(f,'Units','pixels','InnerPosition',iPos,'OuterPosition',oPos);
            p=get(f,'OuterPosition');
            set(f,'Units','normalized');
        end
    end
    function hideMenu(f)
        if nargin < 1 || isempty(f)
            f=gcf;
        end
        set(f, 'MenuBar', 'none');
        set(f, 'ToolBar', 'none');
    end
    function showMenu(f)
        if nargin < 1 || isempty(f)
            f=gcf;
        end
        set(f, 'MenuBar', 'figure');
        set(f, 'ToolBar', 'auto');
    end
    function fit(f,a,dnk)
        if nargin < 1 || isempty(f)
            f=gcf;
        end
        if nargin < 2 || isempty(a)
            a=gca;
        end
        if nargin < 3 || isempty(dnk)
            dnk=1;
        end

        Axis.formatIm();
        [n,m,l]=size(a.Children.CData);
        Fig.setWH(f,m/dnk,n/dnk);
        Axis.fit(a,f);
    end
    function fitClean(f,dnk)
        if nargin < 1 || isempty(f)
            f=gcf;
        end
        if nargin < 2 || isempty(dnk)
            dnk=1;
        end

        Fig.hideMenu(f);
        Fig.fit(f,[],dnk);
    end
    function setXY(varargin)
        if nargin > 0 && isgraphics(varargin{1});
            f=varargin{1};
            varargin=varargin(2:end);
        else
            f=gcf;
        end
        if nargin > 0 numel(varargin{1}) == 2
            X=varargin{1}(1);
            Y=varargin{2}(1);
        else
            if length(varargin) > 0
                X=varargin{1};
            else
                X=4;
            end

            if length(varargin) > 1
                X=varargin{2};
            else
                Y=0;
            end
        end
        set(f,'Position',[X,Y,f.Position(3:4)]);
    end
    function setWH(varargin)
        if isgraphics(varargin{1});
            f=varargin{1};
            varargin=varargin(2:end);
        else
            f=gcf;
        end
        if numel(varargin{1}) == 2
            W=varargin{1}(1);
            H=varargin{1}(2);
        else
            if length(varargin) > 0
                W=varargin{1};
            end

            if length(varargin) > 1
                H=varargin{2};
            else
                H=f.Position(4);
            end
        end
        cl=Fig.tmpSet(f,'Units','pixels');

        Pos=f.Position(1:2);
        set(f,'Position',[Pos,W,H]);
    end
    function names=ls()
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        n=length(FigList);
        names=cell(n,1);
        for i = 1:n
            names{i}=FigList(i).Name;
        end
    end
    function saveAll(typ,bDateDir)
        if nargin < 1
            typ={'fig','epsc','png'};
        end
        if ~iscell(typ)
            typ={typ};
        end
        if nargin < 2 || isempty(bDateDir)
            bDateDir=true;
        end
        if bDateDir
            dire=[getenv('PX_CUR_MEDIA') 'raw' filesep Date.dateStr filesep];
        else
            dire=[getenv('PX_CUR_MEDIA') 'raw' filesep];
        end
        if ~Dir.exist(dire)
            Dir.mk(dire);
        end
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        set(0, 'DefaultFigureRenderer', 'painters');
        for iFig = 1:length(FigList)
              FigHandle = FigList(iFig);
              FigName   = strrep(get(FigHandle, 'Name'),' ','_');
              if isempty(FigName)
                  FigName=num2str(FigHandle.Number);
              end
              for i = 1:length(typ)
                  fname=[dire FigName '.' typ{i}];
                  if strcmp(typ,'fig')
                      savefig(FigHandle, fname);
                  else
                      if ismember(typ{i},{'svg','eps','epsc'})
                          set(FigHandle, 'RendererMode', 'manual','Renderer','painters');
                      else
                          set(FigHandle, 'RendererMode', 'auto');
                      end
                      saveas(FigHandle,fname,typ{i});
                  end
              end
        end
    end
    function saveCur(typ,bDateDir)
        if nargin < 1
            typ={'fig'};
        end
        if ~iscell(typ)
            typ={typ};
        end
        if nargin < 2 || isempty(bDateDir)
            bDateDir=true;
        end
        if bDateDir
            dire=[getenv('PX_CUR_MEDIA') Date.dateStr filesep];
        else
            dire=getenv('PX_CUR_MEDIA');
        end
        if ~Dir.exist(dire)
            Dir.mk(dire);
        end
        FigList = gcf;
        set(0, 'DefaultFigureRenderer', 'painters');
        for iFig = 1:length(FigList)
            FigHandle = FigList(iFig);
            FigName   = strrep(get(FigHandle, 'Name'),' ','_');
            for i = 1:length(typ)
                fname=[dire FigName '.' typ{i}];
                if strcmp(typ,'fig')
                    savefig(FigHandle, fname);
                else
                    if ismember(typ{i},{'svg','eps'})
                        set(FigHandle, 'RendererMode', 'manual');
                    else
                        set(FigHandle, 'RendererMode', 'auto');
                    end
                    saveas(FigHandle,fname,typ{i});
                end
            end
        end
    end
    function getName(nameORopts)
        if nargin < 2 || isempty(YMD)
            bDateDir=false;
            depth=2;
        else
            dateStr=dat;
            depth=1;
        end
        if ischar(nameORopts)
            name=nameORopts;
        elseif iscell(nameORopts)
            name=nameFromOpts(nameORopts{:});
        else
            name=nameFromOpts(nameORopts);
        end
    end
    function [out,num]=exist(nameORopts,YMD)
        if ischar(nameORopts)
            name=nameORopts;
        elseif iscell(nameORopts)
            name=nameFromOpts(nameORopts{:});
        else
            name=nameFromOpts(nameORopts);
        end
        if nargin < 2 || isempty(YMD)
            bDateDir=false;
            depth=2;
        else
            dateStr=dat;
            depth=1;
        end
        if bDateDir
            dire=[getenv('PX_CUR_MEDIA') Date.dateStr(YMD) filesep];
        else
            dire=getenv('PX_CUR_MEDIA');
        end
        [~,files]=Fil.find(dire,[name '\..*'],depth);
        num=numel(files);
        out =num >= 1;
    end
    function rm(nameORopts,YMD)
        % mat, fig, other
        if nargin < 1
            typ={'fig'};
        end
        if nargin < 2 || isempty(YMD)
            bDateDir=false;
            depth=2;
        else
            dateStr=dat;
            depth=1;
        end
        if ischar(nameORopts)
            name=nameORopts;
        elseif iscell(nameORopts)
            name=nameFromOpts(nameORopts{:});
        else
            name=nameFromOpts(nameORopts);
        end
        if bDateDir
            dire=[getenv('PX_CUR_MEDIA') Date.dateStr(YMD) filesep];
        else
            dire=getenv('PX_CUR_MEDIA');
        end
        [~,files]=Fil.find(dire,[name '\..*'],depth);
    end
    function F=sel(name)
        f=findobj('Name',name);
        if isempty(f) && nargout > 0
            F=[];
        elseif isempty(f)
            error('Cannot find figure %s',name);
        end
        h=figure(f.Number);
        if nargout > 0
            F=h;
        end
    end
    function name=nameFromOpts(varargin)
        if nargin == 1
            Opts=varargin{1};
            args=Struct.toCell(Opts);
        else
            args=reshape(varargin,2,numel(varargin)/2)';
        end
        bInd=cellfun(@isnumeric,args(:,2));
        args(bInd,2)=cellfun(@Num.toStr,args(bInd,2),'UniformOutput',false);
        a=args';
        name=strjoin(a(:),' ');
    end
    function [F,bNew]=renew(varargin)
        [bClear,varargin,bSuccess]=Args.getPair('bClear',varargin{:});
        varargin=[varargin 'bClear' false ];
        [f,bNew]=Fig.new(varargin{:});
        if nargout > 0
            F=f;
        end
    end
    function [F,bNew]=newIfBase(varargin)
        s=dbstack;
        if length(s) <= 2
            [F,bNew]=Fig.new(varargin{:});
        else
            F=gcf;
            bNew=false;
        end

    end
    function [F,bNew]=new(varargin)
    % Fig.new()
    %   next unopened
    % f=Fig.new()
    %   return handle
    % Fig.new('myffigname')
    %   select if exists by name or create
    % Fig.new(opts)
    % Fig.new(opts{:})
    %   select if exists by name (generated by params) or create

        if nargin > 0
            [bClear,varargin,bClearG]=Args.getPair('bClear',varargin{:});
            [bFocus,varargin,bFocusG]=Args.getPair('bFocus',varargin{:});
            [Start,varargin,bStartG]=Args.getPair('Start',varargin{:});
        else
            bClearG=false;
            bFocusG=false;
            bStartG=false;
        end
        narg=nargin;
        if ~bClearG
            bClear=true;
        else
            narg=narg-2;
        end

        if ~bFocusG
            bFocus=true;
        else
            narg=narg-2;
        end

        if ~bStartG
            Start=true;
        else
            narg=narg-2;
        end

        bNew=true;
        if narg < 1
        % NO PARAMS
            f=figure(Fig.next);
            return
        elseif narg == 1 && ischar(varargin{1})
        % BY NAME
            name=varargin{1};
            varargin=varargin(2:end);
        else
        % BY OPTIONS
            name=Fig.nameFromOpts(varargin{:});
        end
        f=findobj('Name',name);

        % NEW
        if ~isempty(f)
            bNew=false;
            %close(obj.f);
            fn=get(f,'Number');
            if bClear
                clf(f,'reset');
            end
            f.Name=name;
            if ~bFocus
                set(groot,'CurrentFigure',fn,varargin{:}); % No focus XXX
            else
                f=figure(f,varargin{:});
            end
        else
            fn=get(f,'Number');
            if ~bFocus
                set(groot,'CurrentFigure',fn,varargin{:}); % No focus XXX
            else
                f=figure('Name',name,varargin{:});
            end

        end
        if nargout > 0
            F=f;
        end
    end
    function [num] = next()
        %Finds the next non-existsing image number
        % example Call
        %   Figure(nFn)

        figs=findobj('type','figure');
        if length(figs)==0
            num=1;
            return
        end

        [nums,ind]=sort([figs(:).Number]');

        i=1;
        while true
            if ~ismember(i,nums)
                num=i;
                return
            end
            i=i+1;
        end
    end
    function [varargout]=newU(nameORnum)
        % XXX OLD
        if nargin < 1
            name=Fig.gen_name(0,1);
        elseif isnumeric(nameORnum)
            num=nameORnum;
            name=Fig.gen_name(num,1);
        else
            name=nameORnum;
        end
        num=Fig.nextU(name);
        f=figure(num);
        f.Name=name;
        if nargout >= 1
            varargout{1}=f;
            if nargout >= 2
                varargout{2}=f.Name;
            end
        else
            figure(num);
        end

    end
    function name=genName()
        % XXX OLD
        name=Fig.gen_name(1);
    end
    function num=sameORnext(name)
        % XXX OLD
        num=Fig.find(name);
        if isempty(num)
            num=Fig.next();
        end
    end
    function [num] = nextU(name)
        % XXX OLD
        %Finds the next non-existsing image number
        % example Call
        %   Figure(nFn)

        names=Fig.get_names();

        if isempty(names)
            num=1;
        elseif ismember(name,names)
            num=find(ismember(names,name),1,'first');
        else
            num=Fig.next();
        end

    end
%%- CLOSE
    function reClose(re)

    end
%%- FORMAT
    function ax=formatIm(varargin)
        Error.warnSoft('Call formatIm from axis!');
        Axis.formatIm(varargin{:});
    end
    function format(varargin)
        Error.warnSoft('Call format from axis!');
        Axis.format(varargin{:});
    end
    function [] = move(figNum,bMoveElseWhere,bMoveElseWhereAfter,Xdestination,Ydestination,bPrimary)
        if ~exist('bPrimary','var') || isempty(bPrimary)
            bPrimary=1;
        end
        if ~exist('figNum','var') || isempty(figNum)
            fig=get(gcf);
            figNum=fig.Number;
        end
        if ~exist('bMoveElseWhere','var') || isempty(bMoveElseWhere)
            bMoveElseWhere=0;
        end
        if ~exist('bMoveElseWhereAfter','var') || isempty(bMoveElseWhereAfter)
            bMoveElseWhere=0;
        end

        dire=curFunctionDir;
        cmd=[dire filesep 'moveFigure.sh ' num2str(figNum) ' ' num2str(bMoveElseWhere) ' ' num2str(bMoveElseWhereAfter) ' ' Xdestination ' ' Ydestination ' ' bPrimary];
        Sys.run(cmd);
    end
    function []=setRes(figNumOrObj,mult)
        %MAXIMIZE SUBPLOTS
        if isnumeric(figNumOrObj)
            fig=figure(figNumOrObj);
        else
            fig=figNumOrObj;
        end
        if ~exist('mult','var') || isempty(mult)
            mult=1;
        end
        figure(fig.Number);
        [dimsRC,subsRC,f]=Plot.getDims(fig.Number);
        maxResXY=[1 1];
        dimsXY=fliplr(dimsRC);
        width=maxResXY(1)/dimsXY(1);
        height=maxResXY(2)/dimsXY(2);
        if numel(f.Children)==1
            Axis.setRes(mult);
        else
            for i = 1:length(f.Children)
                subPlot(dimsRC,subsRC(i,1),subsRC(i,2));
                hold on;
                X=(subsRC(i,2)-1)*width;
                Y=(dimsRC(1)-subsRC(i,1))*height;
                %subsRC(i,:):
                %[X Y]
                set(gca,'OuterPosition',[X,Y,width,height]);
            end
        end
    end
    function [dimsXY,subsRC,f] = getDims(figNumOrObj)
        % GET FIGURE PARAMS
        if isnumeric(figNumOrObj)
            fig=figure(figNumOrObj);
        else
            figNumber=figNumOrObj;
            fig=figure(figNumber);
        end

        set(gcf,'Units','Normalized');
        f=get(gcf);
        X=zeros(length(f.Children),1);
        Y=zeros(length(f.Children),1);
        h=zeros(length(f.Children),1);
        %W=zeros(length(f.Children),1);
        %H=zeros(length(f.Children),1);
        bSGTitle=true;
        sgInd=[];
        for i = 1:length(f.Children)
            child=f.Children(i);
            if isa(child,'matlab.graphics.illustration.subplot.Text')
                bSGTitle=true;
                sgInd=[sgInd; i];
                continue
            end
            X(i)=child.Position(1);
            Y(i)=child.Position(2);
            h(i)=child;
            %H(i)=f.Children(i).Position(3);
            %W(i)=f.Children(i).Position(4);
        end
        X(sgInd)=[];
        Y(sgInd)=[];
        [~,~,rX]=unique(X);
        [~,~,rY]=unique(1-Y);
        dimsXY=[max(rX) max(rY)];
        subsRC=[rY rX];
    end
    function save(filename,addenum,bPrompt)
        % fig
        if ~exist('addenum','var') || isempty(addenum)
            addenum='';
        else
            addenum='_';
        end
        if ~exist('bPrompt','var')
            bPrompt=[];
        end

        f=gcf;
        out=Fig.getInfo();

        % defaultName
        if ~exist('filename','var') || isempty(filename)
            defName=[Fig.get_default_name() addenum];
            if isempty(bPrompt)
                bPrompt=true;
            end
        end

        % change directory
        currentFolder = pwd;
        c = onCleanup(@()cd(currentFolder));
        dire=Env.var('PX_CUR_MEDIA');
        cd(dire);


        if bPrompt
            [filename,extOut,exitflag] = Fig.imPut('*.fig',[defName '.fig']);
            filename=filename(1:end-4);
            if exitflag
                return
            end
        elseif ~bPrompt && ~exist('filename','var') || isempty(filename)
            filename=defName;
        end

        %saveas(fig,filename,formattype)
        savefig(f,[filename '.fig']);
        Fig.saveTypes(f,filename);

    end
    function info=getInfo(f)
        if ~exist('f','var') || isempty(f)
            f=gcf;
        end
        if ~exist('bRec','var') || isempty(bRec)
            bRec=false;
        end

        info=struct();
        [info.dimsXY,info.subsRC]=Fig.getDims(f);
        [info.titles,info.titleInd,titleSubNorm]=Fig.get_titles(f);
        [info.sgTitles,info.sgInd]=Fig.get_sgTitles(f);
        info.root=Fig.get_root_caller_name();
        info.last=Fig.get_last_caller_name();
    end
%% SG TITLES
    function SGTitle(varargin)
        f=gcf;
        if nargin > 0 && ischar(varargin{1})
            text=varargin{1};
            varargin=varargin(2:end);
        elseif nargin > 1 && ischar(varargin{2})
            text=varargin{2};
            varargin=[varargin(1) varargin(3:end)];
        end
        sgTitles=Fig.get_sgTitles(f);

        if numel(sgTitles)>0
            sgTitles(1).String=text;
        else
            sgtitle(text,varargin{:});
        end
    end
    function supTitle(varargin)
        Fig.SGTitle(varargin{:});
    end
    function suptitle(varargin)
        Fig.SGTitle(varargin{:});
    end
    function sgtitle(varargin)
        Fig.SGTitle(varargin{:});
    end
    function [sgTitles,sgInd]=get_sgTitles(f)
        sgInd=[];
        j=0;
        for i = 1:length(f.Children)
            j=j+1;
            child=f.Children(i);
            if isa(child,'matlab.graphics.illustration.subplot.Text')
                sgInd=[sgInd; i];
            end

        end
        sgTitles={};
        if numel(sgInd) > 0
            sgTitles=f.Children(sgInd);
        end
    end
%% TITLES
    function [titles,titleInd,titleSubNorm]=get_titles(f)
        titleInd=[];
        titleIndNorm=[];
        j=0;
        for i = 1:length(f.Children)
            j=j+1;
            child=f.Children(i);
            if isa(child,'matlab.graphics.axis.Axes') && ~isempty(child.Title.String)
                titleInd=[titleInd; i];
                titleIndNorm=[titleIndNorm; j];
            elseif isa(child,'matlab.graphics.illustration.subplot.Text')
                j=j-1;
            elseif  ~isa(child,'matlab.graphics.axis.Axes')
                error(['unhandled type ' class(child)]);
            end

        end
        if numel(titleInd) > 0
            [dimsXY,~]=Fig.getDims(f);
            titles=transpose({f.Children(titleInd).Title});
            titles=cellfun(@(x) x.String, titles, 'UniformOutput',false);
            [n,m]=ind2sub(fliplr(dimsXY),titleIndNorm);
            titleSubNorm=[n,m];
        else
            titles={};
            titleSubNorm=[];
        end
    end
%% MAIN TITLE
    function out=get_default_name(info)
        if ~exist('f','var') || isempty(f)
            f=gcf;
        end
        if ~exist('info','var') || isempty(info)
            info=Fig.getInfo(f);
        end

        if numel(info.sgTitles)==1
            out=info.sgTitles{1}.String;
            return
        elseif prod(info.dimsXY)==1
            out=info.titles{1}.String;
        elseif strcmp(info.root,info.last)
            out=info.root;
        else
            out=[info.root '__' info.last];
        end
    end
%% CALLER
    function out=get_root_caller_name()
        out=dbstack;
        out=out(end).name;
    end
    function src=get_last_caller_name()
        dbs=dbstack;
        thisFile=[mfilename '.m'];

        src='BASE';
        for i = 1:length(dbs)
            file=dbs(i).file;
            if ~strcmp(file,thisFile);
                src=dbs(i).name;
                break
            end
        end
    end
%% SAVE
    function saveTypes(fig, fulldir, varargin)

        if length(varargin) < 1
            try
                types=strsplit(Env.var('FIG_TYPES'));
            catch
                types={'eps','png'};
            end
        else
            types=varargin;

        end
        for i = 1:length(types)
            saveas(fig,[fulldir '.' types{i}],types{i});
        end
    end
    function autoSave(fig,opts,varargin)
    %function save(fig,opts,varargin)
    % varargin=types

        % FIG
        if ~exist('fig','var') || isempty(fig)
            fig=gcf;
        end

        % OPTS
        if ~exist('opts','var') || isempty(opts)
            opts=struct();
        elseif ~isfield(opts,'bSave')
            opts.bSave=1;
        elseif opts.bSave==0
            return
        end

        % FNAME
        name=get_last_caller();

        % DIRE
        fulldir=[opts.saveLocation fname filesep];
        if ~exist(fulldir,'dir')
            mkdir(fulldir);
        end

        Fig.saveTypes(fig, varargin{:});
    end


%end
%methods(Static, Access=private)
    function [filename,extOut,exitflag] = imPut(filterSpec,varargin)
        %IMPUTFILE Save Image dialog box.
        %   [FILENAME, EXT, USER_CANCELED] = IMPUTFILE displays the Save Image
        %   dialog box for the user to fill in and returns the full path to the
        %   file selected in FILENAME.  Additionally the file extension is returned
        %   in EXT.  If the user presses the Cancel button, USER_CANCELED will
        %   be TRUE. Otherwise, USER_CANCELED will be FALSE.
        %
        %   The Save Image dialog box is modal; it blocks the MATLAB command line
        %   until the user responds.
        %
        %   See also IMFORMATS, IMTOOL, IMGETFILE, IMSAVE.

        %   Copyright 2007-2014 The MathWorks, Inc.

        % Get filter spec for image formats

        if ~exist('filterSpec','var') || isempty(filterSpec)
            filterSpec = Fig.create_fig_filter_spec();
        elseif ischar(filterSpec)
            filterSpec={filterSpec};
        end

        msg=getString(message('images:fileGUIUIString:putFileWindowTitle'));

        [filename, pathname, filterindex] = uiputfile(filterSpec,msg,varargin{:});

        exitflag = (filterindex == 0);

        if ~exitflag

            filename = fullfile(pathname,filename);
            % If there are variants of extension name, return first extenstion in list
            selectedExt = filterSpec{filterindex,1};
            ind = strfind(selectedExt,';');
            if isempty(ind)
                extOut = selectedExt(3:end);
            else
                extOut = selectedExt(3:ind-1);
            end

        else
            filename = '';
            extOut = '';
        end
    end


    function filterSpec = create_fig_filter_spec()
        %   Copyright 2007-2014 The MathWorks, Inc.

        [desc, ext , ~, write_fcns] = iptui.parseImageFormats();
        nformats = length(desc);

        % Grow filter spec dynamically to avoid need to hardcode knowledge of
        % number of supported file formats or formats with write_fcns.
        filterSpec = {};

        % Formats that we want to disable in imputfile
        excluded_exts = {'gif','hdf','pcx','pnm','xwd','dcm','rset'};

        for i = 1:nformats

            format_is_writable = ~isempty(write_fcns{i});
            format_is_excluded = any(ismember(ext{i},excluded_exts));

            if format_is_writable && ~format_is_excluded
                thisExtension = ext{i};
                numExtensionVariants = length(thisExtension);
                thisExtensionString = strcat('*.',thisExtension{1});
                for j = 2:numExtensionVariants
                    thisExtensionString = strcat(thisExtensionString,';*.',thisExtension{j});
                end

                % Populate individual file extension and descriptions
                filterSpec{end+1,1} = thisExtensionString; %#ok<AGROW>
                filterSpec{end,2} = desc{i};

            end

        end
    end
%%-  UTIL
    function names=get_names()
        figs=findobj('type','figure');
        if length(figs)==0
            names={};
            return
        end

        [nums,ind]=sort([figs(:).Number]');
        names={figs(:).Name}';
        names=names(ind);
    end
end
methods(Static, Access=private)
    function prj=get_caller_name(layer)
        if nargout > 1
            layer=0;
        end
        prj=Prj.name([],2+layer);
    end
    function name=gen_name(num,layer)
        [names,funs,lines]=Prj.names(layer+1);
        base=strjoin(strcat( names' ,':', funs' ));
        if num ==0
            num=lines(end);
        end
        num=[num2str(num) '  '];
        name=[num base];
    end
    %function num=get_num(names,base)
    %    allnames=Fig.get_names();
    %    inds=Str.RE.ismatch(allnames,['^[0-9]+__' base]);
    %    nums=str2double(strrep(allnames(inds),['__' base],''));
    %    i=0;
    %    while true
    %        i=i+1;
    %        if ~ismember(i,nums)
    %            num=i;
    %            break
    %        end
    %    end
    %end
    function cl=tmpSet(obj,fld,val);
        oval=obj.(fld);
        obj.(fld)=val;
        cl=onCleanup(@() set_fun(obj,fld,oval));
        function set_fun(obj,fld,oval)
            obj.(fld)=oval;
        end
    end
end
end
