function color = lastColor();
    a = get(gca,'Children');
    color=a(1).Color;
