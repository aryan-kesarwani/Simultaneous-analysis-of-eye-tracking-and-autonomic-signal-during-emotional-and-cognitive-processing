function plotLed(cond,col,pos)

xl = get(gca,'xlim');
yl = get(gca,'ylim');
b = (yl(2)-yl(1))/12;

if pos == 'b'
    for i = 1:length(cond)
        text(xl(2),yl(1)+(length(cond)+1-i)*b,cond{i},'HorizontalAlignment','right','Color',col{i},'FontWeight','bold')
    end
elseif pos == 't'
    for i = 1:length(cond)
        text(xl(2),yl(2)-(i)*b,cond{i},'HorizontalAlignment','right','Color',col{i},'FontWeight','bold')
    end
elseif pos == "lt"
    for i = 1:length(cond)
        vl = xl(2)-xl(1);
        text(xl(1)+0.1*vl,yl(2)-(i)*(yl(2)-yl(1))/10,cond{i},'HorizontalAlignment','left','Color',col{i},'FontWeight','bold')
    end
elseif pos == "lb"
    for i = 1:length(cond)
        vl = xl(2)-xl(1);
        text(xl(1)+0.1*vl,yl(1)+(length(cond)+1-i)*b,cond{i},'HorizontalAlignment','left','Color',col{i},'FontWeight','bold')
    end
end
end