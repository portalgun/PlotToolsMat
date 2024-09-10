classdef Hold < handle
properties
    g
    f
    bHold
end
methods
    function obj=Hold(g,f)
        if nargin < 1 || isempty(g)
            g=gca;
        end
        if nargin < 2 || isempty(f)
            f=gcf;
        end
        obj.g=g;
        obj.f=f;
        obj.bHold=ishold();
    end
    function delete(obj)
        set(0,'CurrentFigure',obj.f);
        set(obj.f,'CurrentAxes',obj.g);
        if obj.bHold
            hold on;
        else
            hold off;
        end
    end
end
end
