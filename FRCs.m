classdef FRCs < handle
% NOTE VALS ARE IN PARENT
properties
    SP
    figs

    FRC=[0 0 0] % N
    frc=[0 0 0] % sel
    fig         % sel
    sp

    sFRCFlds    % flds ordered (static)
    FRCFlds     % flds in spec plot order (variable)
end
properties(Access=protected)
    Parent
end
methods
    function obj=FRCs(parent,figNames,RC,spArgs, sFlds,flds)
        if nargin >= 6
            obj.setFlds(sFlds,flds);
        end
        obj.Parent=parent;

        if numel(RC)==3
            obj.FRC=RC;
        elseif numel(RC)==2
            obj.FRC=[numel(figNames) RC];
        end

        % SUBPLOT
        obj.SP=SubPlot.empty;
        for i = 1:obj.FRC(1)
            obj.figs(i)=Fig.new(figNames{i});
            obj.SP(i)=SubPlot(obj.FRC(2:end),spArgs{:});
            hold off;
            axis off;
        end
    end
    function sel(obj,ir,ii,ij)
        irLast=obj.frc(1);
        obj.frc=[ir,ii,ij];
        if ir ~= irLast
            figure( obj.figs(ir));
            obj.fig=obj.figs(ir);
            obj.sp =obj.SP(ir);
        end
        obj.sp.sel([ii,ij]);
    end
%- FLDS
    function setFlds(obj,sFlds,flds)
        obj.sFRCFlds=sFlds;
        obj.FRCFlds=flds;

        % XXX ?
        %for i = 1:3
        %    obj.FRC(i)=numel(obj.Parent.(fld));
        %end
    end
    function varargout=unpackFRC(obj)
        bTest=false;
        varargout=cell(1,3);
        sflds=obj.sFRCFlds;
        for i = 1:3
            sfld=sflds{i};
            ind=find(strcmp(sfld,obj.FRCFlds));
            fld=obj.FRCFlds{ind};
            vals=obj.Parent.(fld);
            val=vals(obj.frc(ind));

            %% TEST
            if bTest && all(obj.frc==[3 2 1]);
                disp([num2str(i) ' ' fld ' ' num2str(val)])
            end

            varargout{i}=val;
        end
    end
end
end
