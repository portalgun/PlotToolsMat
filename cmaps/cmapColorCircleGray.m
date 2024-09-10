function cmap = cmapColorCircleGray(saturation)

% function cmap = cmapColorCircleGray
%
% color circle color map: BLACK -> GRAY -> WHITE -> GRAY -> BLACK
%
% saturation:  saturation of color wheel. value must lie on (0 1]
% %%%%%%%%%%%%%%%
% cmap:   256x3 colormap
    
if ~exist('saturation','var') || isempty(saturation)
   saturation = 1; 
end
if saturation >  1, saturation = 1.0; end
if saturation <= 0, saturation = 0.2; end

vals = [0 0 0; 0.5 0.5 0.5; 1 1 1; 0.5 0.5 0.5; 0 0 0];
cmap = interp1(1:5, vals, linspace(1,5,256),'pchip');

cmap = saturation.*cmap + (1-saturation)./2;

