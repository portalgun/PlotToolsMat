function newmap = cmapRBR(X)
% function cmapRBR 
% 
%   Red, blue, red
%   
%   cmapBWR(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  cmapBWR, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(cmapBWR(256)), colorbar
% 
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(cmapBWR), colorbar
% 
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(cmapBWR), colorbar
% 
%   figure
%   surf(peaks)
%   colormap(cmapBWR)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.


if nargin < 1
   X = 1:64;
end

bottom    = [1.0 0.0 0.0];
botmiddle = [0.0 0.8 0.0];
middle    = [0.0 0.0 1.0];
topmiddle = [0.0 0.8 0.0];
top       = [1.0 0.0 0.0];

M = [bottom; botmiddle; middle; topmiddle; top];

newmap = interp1([1:size(M,1)]',M,linspace(1,size(M,1),length(X)));
