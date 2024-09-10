function cmap = cmapColorCircle(saturation)

% function cmap = cmapColorCircle
%
% color circle color map: BLUE -> GREEN -> YELLOW -> RED -> BLUE
%
% saturation:  saturation of color wheel. value must lie on (0 1]
% %%%%%%%%%%%%%%%
% cmap:   256x3 colormap
    
if ~exist('saturation','var') || isempty(saturation)
   saturation = 1; 
end
if saturation >  1, saturation = 1.0; end
if saturation <= 0, saturation = 0.2; end
    

% cmap = interp1(1:5, [ 0 0 1; 0 1 0; 1 1 0; 1 0 0; 0 0 1], linspace(1,5,256),'pchip');
cmap = interp1(1:5, [ 0 0 1; 0 1 0; 1 1 0; 1 0 0; 0 0 1], linspace(1,5,256),'pchip');

cmap = saturation.*cmap + (1-saturation)./2;