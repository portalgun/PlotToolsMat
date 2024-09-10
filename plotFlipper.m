function t =plotFlipper(t,n)
% Example Call:
%    i=1;
%    while i <= length(Inds);
%           imagesc(Image(:,:,Inds(t)));
%           i=plotFlipper(i);
%    end
keys.space=32; %Space
keys.L=28; %LEFT
keys.R=29; %RIGHT
keys.U=30; %UP
keys.G=103;
keys.D=31; %DOWN
keys.r=114; %DOWN
keys.esc=27;
keys.zero=48;

drawnow
while true %LOOP UNTIL VALID KEY IS PRESSED
    kk=waitforbuttonpress;
    k = double(get(gcf,'CurrentCharacter')); % NOTE uncomment to see keycodes on keypress

    if isempty(k)
        continue
    end
    if  k==keys.R | k==keys.space
        t=t+1;
        break
    elseif k==keys.L && t > 1
        t=t-1;
        break
    elseif k==keys.G && exist('n','var') && ~isempty(n)
        t=Input.range(n);
        break
    end
end
