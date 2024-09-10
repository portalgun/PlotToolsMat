function plotCmap(cmap,bInvert)
    if ~exist('bInvert','var') || isempty(bInvert)
        bInvert=0;
    end
    if bInvert
        cmap=invertRGB(cmap);
    end
    ind=repelem(1:size(cmap,1),3,1);
    cmap=permute(cmap(ind,:),[1 3 2]);
    cmap=repmat(cmap,1,size(cmap,1));
    imshow(cmap)
end
