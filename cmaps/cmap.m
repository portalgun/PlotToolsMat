classdef cmap < handle
% cmap(name,n)
properties
    name
    fullmap
    map
    n % size of whatever plotting
    N % size of cmap

    color
    i=0;
    bInvert=0
end
methods
    function obj=cmap(name,n,bInvert)
        obj.name=name;
        obj.fullmap=feval(name);
        obj.N=size(obj.fullmap,1);
        if exist('bInvert','var')
            obj.bInvert=bInvert;
        end
        if obj.bInvert
            obj.invert_cmap();
        end
        if exist('n','var')
            obj.n=n;
            obj.shrink_map();
        end
    end
    function obj=shrink_map(obj,n)
        if ~exist('n','var') || isempty(n)
            n=obj.n;
        else
            obj.n=n;
        end
        inds=round(linspace(1,obj.N,n));
        obj.map=obj.fullmap(inds,:);
    end
    function obj=c(obj)
        obj.reset();
    end
    function obj=reset(obj)
        obj.i=0;
        obj.color=[];
    end
    function color=u(obj)
        color=obj.update();
    end
    function color=update(obj,ind)
        if ~exist('ind','var') || isempty(ind)
            obj.update_helper();
        else
            obj.i=ind;
            obj.color=obj.map(obj.i,:);
        end
        color=obj.color;
    end
    function colors=output(obj,inds)
        colors = zeros(length(inds),3);
        for i = 1:length(inds) 
            colors(i,:)=obj.update(inds(i));
        end
    end
    function obj=update_helper(obj)
        obj.i=obj.i+1;
        obj.color=obj.map(obj.i,:);
    end
    function obj=invert_cmap(obj)
        obj.fullmap=repmat([1 1 1],obj.N,1)-obj.fullmap;
    end
end
end
