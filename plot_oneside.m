
function p = plot_oneside(data12,areav,vside,fh,labelst,col,shiftup,nylabel,nxlabel)
[nb,xout] = hist(data12, 100);

figure(fh);
nnb=nb/sum(nb);
p=plot(xout,nb/sum(nb),col,'linewidth',1.5);
[v1,ind1]=find(vside*xout>vside*0.5);
hold on
if numel(ind1)
    if vside>0
        Dnnb = nnb(ind1(1))-nnb(max(1,ind1(1)-1));
        Dxout = xout(ind1(1))-(xout(max(1,ind1(1)-1)));
        nnb0 = nnb(ind1(1))-(xout(ind1(1))-0.5)*Dnnb/Dxout;
        harea=area([0.5 xout(ind1)],[nnb0 nnb(ind1)]);
    else
        Dnnb = nnb(min(ind1(end)+1,length(nb)))-nnb(ind1(end));
        Dxout = xout(min(ind1(end)+1,length(nb)))-(xout(ind1(end)));
        nnb0 = nnb(min(ind1(end)+1,length(nb)))-(xout(ind1(end))-0.5)*Dnnb/Dxout;
        harea=area([0 xout(ind1)],[0 nnb(ind1)]);
    end
else
    if vside>0
        harea = area([0.5 vside*xout(1,1)+1],[0 0])
    else
        harea = area([vside*xout(1,1) 0.5 ],[0 0])
    end
end
% harea=area(xout(ind1),nnb(ind1));
set( harea, 'FaceColor', col, 'EdgeColor',col);
alpha(.3);

[v,ind]=max(nb/sum(nb));
ylabel(nylabel);
xlabel(nxlabel);
a = axis;
width = a(2)-a(1);
ht = a(4)-a(3);
text(shiftup*width,ht*0.5,[labelst, sprintf('%.3f', round(areav*1000)/1000)],'Color', col);
% textbp([labelst, sprintf('%.3f', round(areav*1000)/1000)])

end