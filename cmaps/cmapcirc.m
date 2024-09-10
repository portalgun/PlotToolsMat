function Y = cmapcirc(varargin)
% function Y = cmapcirc(varargin)
%
%    Y = cmapcirc(cmapType) creates a circular colormap type cmapType. Circular colormaps
%    are built so that the first and last colors are the same. This kind of
%    colormap is useful for ploting angle images which have a circular behavior
%    (0 = 2*pi).
%    
%    Y = COLORMAPC(cmapType,N) creates the circular colormap cmapType with size N.
%
%    cmapType = 1 -> blk, yel, red, gre, blk        cmapType = 2 -> cya, blu, mag, cya (cool)
%    cmapType = 3 -> blu, cya, gre, yel, blu        cmapType = 4 -> blk, red, yel, blk (hot)
%

switch nargin
case 0
   cmapType = 1;
   N = 256;
case 1
   cmapType = varargin{1};
   N = 256;
case 2
   cmapType = varargin{1};
   N = varargin{2};
otherwise      
      error('Too many parameters.')
end

Y = zeros(N,3);

switch cmapType
case 1 
   Y = cmapc(N,[0 0 0],[1 0 0],[1 1 0],[0 1 0]);
case 2 
   Y = cmapc(N,[0 1 1],[0 0 1],[1 0 1]);
case 3 
   Y = cmapc(N,[0 0 1],[0 1 1],[0 1 0],[1 1 0]);
case 4 
   Y = cmapc(N,[0 0 0],[1 0 0],[1 1 0]);
otherwise
   error('Unknown colormap type.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = cmapc(varargin)
%if  nargin < 3
%   error('Too few inputs.')
%end
N = varargin{1};
C = zeros(nargin,3);

C = cat(1,varargin{2:end},varargin{2});

Y = zeros(N,3);
ls = 1;
r = N/(size(C,1)-1); % Length of each color transition.
for i = 1 : size(C,1)-1
   li = ls;
   ls = round(i*r);
   d = ls - li + 1;
   Y(li:ls,:) = [linspace(C(i,1),C(i+1,1),d); linspace(C(i,2),C(i+1,2),d); linspace(C(i,3),C(i+1,3),d)]';
end