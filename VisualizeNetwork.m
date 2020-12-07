function VisualizeNetwork(X,p,U,pos)
    
    % X is array of all the nodal quantities
    totalnum = 2*p.row*p.col;
    gap = 0.25;
    biasgap = 0.3;
    
    % currents
    sc = 1:1:totalnum;
    tc = totalnum+1:1:2*totalnum;
    currents = zeros(1,totalnum);
    for i = 1:totalnum
        if mod(i,2)==1
            currents(1,i) = abs((X(i)-p.v1)/p.REmitter); 
        else
            currents(1,i) = abs((X(i)-p.v2)/p.RCollector);
        end
    end
    newcurrents = reshape(currents,7,[]);
    sumc = sum(newcurrents,1); % sum up all currents in each bias line
    
    % Drawing positions of all the nodes
    x =[];
    y = [];
    xc = [];
    for j = 1:p.col
        x1 = repmat([(j-1), (j-1)+gap],1,p.row);
        xc1 = repmat([(j-1)-biasgap, (j-1)+gap+0.3],1,p.row);
        x = [x x1];
        xc = [xc xc1];
        y11 = repmat(0:-1:-p.row+1,2,1);
        y1 = reshape(y11,1,[]);
        y = [y y1];
    end
    
    G = digraph(sc,tc,currents); % arrow for current dir
    Gc = flipedge(G,1:2:totalnum); 

    yc = [y y];
    xc = [x xc];
    
    % Line widths are prop to currents, node size prop to nodal voltage
    LWidths = 3*abs(Gc.Edges.Weight/max(Gc.Edges.Weight))+0.5;
    MSize = [X; zeros(length(X),1)];
    MSize(1:2:length(X))=abs(MSize(1:2:length(X))-5);
    MSize = abs(MSize) * 20 + 0.5;
    LWidths(LWidths>5)=5;
    % reduce node growing rate after size 10
    MSize(MSize>10) = 10+2*log(MSize(MSize>10)-9); 
%     MSize(MSize>10) = 10
    
    imgfile = strcat('bowtie49',pos,'.jpg'); % load image
    array = imread(imgfile);
    imagesc([-1.2 7.05], [-7.15 1.3], array); % scale them to overlay
    hold on
%     openfig('firstfig.fig','reuse')
    pg = plot(Gc,'XData',xc,'YData',yc,'NodeLabel',[],'LineWidth',LWidths, 'EdgeColor','r','MarkerSize',MSize);
    x = -1.2+0.9;
    y = 1.1;
    sumc=round(sumc,2);
    gap = gap+0.04;
    text(-1.2,1.1-0.2,'Current from','Color','r','FontSize',15)
    text(-1.2,1.1,'bias line','Color','r','FontSize',15)
    
    % print currents from bias line
    for jj = 1:p.col  
        text(x,y,[num2str(sumc(2*jj-1))],'Color','r','FontSize',17)
        x = x+0.72-gap/2;
        text(x,y,[num2str(sumc(2*jj))],'Color','r','FontSize',17)
        x = x+3*gap/2;
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'visible','off')
    
   
   



