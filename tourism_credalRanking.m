clc
clear


load tourismData
data = Weights;
cri_no = length(nameOfCriteria);

nameOfCriteria = {'Tangibles','Reliability','Responsiveness','Assurance','Empathy'};

probs = zeros(cri_no);
for i=1:cri_no
    for j=i+1:cri_no
        
        %{
        [~,p] = ttest_bf(avgDev_array(i,j)/avgDev_array(j,i),dm_no);
        
        if avgDev_array(i,j) > 1
            probs(i,j) = p;
        else
            probs(j,i) = p;
        end
        %}
        
        [p, h] = isignrank(log(data(:,i)),log(data(:,j)));
        probs(i,j) = mean(p);
        probs(j,i) = 1 - probs(i,j);
    end
end


probs(probs < .5) = 0; 

nonZeroProb = length(find(probs > 0.5));
edgeTable = zeros(nonZeroProb,3);
counter = 1;
for i=1:size(probs,1)
    for j=1:size(probs,1)
        if probs(i,j) > 0.5
           edgeTable(counter,:) = [i j probs(i,j)];
           counter = counter+1;
        end
    end
end

figure()
directedGraph = probs > 0.5;

set(0,'defaultAxesFontSize', 14,'defaultAxesFontWeight','b')
G = digraph(directedGraph);
p = plot(G,'k','NodeLabel',nameOfCriteria,'EdgeLabel',round(edgeTable(:,3),2));
p.Marker = 's';
p.NodeColor = 'r';
p.MarkerSize = 7;
p.EdgeFontSize = 14;
p.EdgeFontWeight = 'Bold';
p.ArrowSize = 13;
%p.LineStyle = '--';
axis off
NameArray = {'linewidth','markersize'};
ValueArray ={3,10};
set(p,NameArray,ValueArray);
set(gcf,'color','w');

p.NodeLabel = '';
xd = get(p, 'XData');
yd = get(p, 'YData');
text(xd, yd, nameOfCriteria, 'FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','bottom')

C = {[0	255	255]/255,'r','g','b',[238 154 0]/255,[139 58 58]/255,[.2 .5 .7],[.6 .2 .7],[255 0 255]/255,[.5 .7 .2]};
C = C(1:cri_no);
for j=1:cri_no
    T = [];
    for i=1:size(G.Edges.EndNodes,1)
        if G.Edges.EndNodes(i,1) == j
           %T = [T;G.Edges.EndNodes(j,1:2)];
            highlight(p,G.Edges.EndNodes(i,1:2),'EdgeColor',C{j},'LineWidth',3)
        end
    end
    
end

    
    

