function cmap=cmap5(minmaxd)

mind = min(minmaxd);
maxd = max(minmaxd);
% return a colormap the has mind-maxd entries

cmap=[...
        zeros(1,-mind),          0, ones(1,maxd)          ;...%red
        [1:-1/(-mind-1):0].^0.3, 0, [0:1/(maxd-1):1].^0.7 ;...%green
        [0:1/(-mind-1):1].^0.3,  0, zeros(1,maxd)         ;...%blue
        ]';
