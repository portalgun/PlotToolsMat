function varargout = subPlot(RC,r,c,bSimplify, xticklbls,yticklbls)
%function varargout = subPlot(RC,r,c,bSimplify,xlabel,ylabel)
%Sane subplot.
i=sub2ind(fliplr(RC),c,r);
ax1=subplot(RC(1),RC(2),i);

nOutputs=nargout;
varargout = cell(1,nOutputs);
for k = 1:nOutputs;
    varargout{k}=ax1;
end

if ~exist('bSimplify','var') || bSimplify==0
    return
end

R=RC(1);
C=RC(2);
if r~=R
    xticklabels([]);
elseif ~exist('xticklbls','var')
    return
else
    xticks(xticklbls);
    xticklabels(xticklbls);
end

if c~=1
    yticklabels([]);
elseif ~exist('yticklbls','var')
    return
else
    yticks(yticklbls);
    yticklabels(yticklbls);
end

