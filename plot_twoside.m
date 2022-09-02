
% --------------------------------------------------------------
function p = plot_twoside(data12,xl,xr,fh,col,shiftup,nylabel,nxlabel)
[nb,xout] = hist(data12, 80);

figure(fh)
nnb=nb/sum(nb);
p=plot(xout,nb/sum(nb),col,'linewidth',1.5);
hold on
ylabel(nylabel);
xlabel(nxlabel);
%             line([xl xr],[0 0],'Color',col,'LineWidth',10);
a = axis;
width = a(2)-a(1);
ht = a(4)-a(3);
text(xl+0.0*width,ht*0.03*shiftup,sprintf('%.3f', round(xl*1000)/1000))
text(xr+0.0*width,ht*0.03*shiftup,sprintf('%.3f', round(xr*1000)/1000))

[v1,ind1]=find(xout>xl);
[v1,ind2]=find(xout<xr);
ind1=intersect(ind1,ind2);
hold on
if isempty(ind1)==0
harea=area(xout(ind1),nnb(ind1));
set( harea, 'FaceColor', col)
alpha(.3)
end

end
