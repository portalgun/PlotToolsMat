function Cmap=invertRGB(Cmap)
if isa(Cmap,'cmap')
    Cmap=Cmap.invert_cmap();
    return
end
bFlag=0;
if all(isint(Cmap))
    Cmap=Cmap./255;
    bFlag=1;
end
Cmap=repmat([1 1 1],size(Cmap,1),1)-Cmap;
if bFlag
    Cmap=Cmap.*255;
end
