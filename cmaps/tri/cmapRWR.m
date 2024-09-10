function newmap = cmapRWR(m)
% function cmapRWR   
% 
%   Red, white, and red color map.
%   
%   cmapRWR(M) returns an M-by-3 matrix containing a red to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for plotting images surface plots 
%   of errors. Positive and negative values with the same absolute value have 
%   the same color
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(cmapRWR(256)), colorbar
% 
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(cmapRWR), colorbar
% 
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(cmapRWR), colorbar
% 
%   figure
%   surf(peaks)
%   colormap(cmapRWR)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

newmap = cmapBWR(m);

newmap(1:(size(newmap,1)./2),:) = flipud(newmap( (size(newmap,1)./2+1):end,:));