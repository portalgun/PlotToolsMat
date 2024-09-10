function out=nextColor(h)
    if nargin < 1 || isempty(h)
        h=gca;
    end
    colors = get(h,'ColorOrder');
    index  = get(h,'ColorOrderIndex');
    if index > size(colors,1);
        index = 1;
    end
    out = colors(index,:);
