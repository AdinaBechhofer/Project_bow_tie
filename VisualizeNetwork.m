function VisualizeNetwork(X,p,U)
   
    
    % X is array of all the nodal quantities
    totalnum = 2*p.row*p.col;
    gap = 0.25;
    biasgap = 0.3;
    % currents
    sc = 1:1:totalnum;
    tc = totalnum+1:1:2*totalnum;
    currents = zeros(1,totalnum);
    for i = 1:totalnum
        if mod(i,1)==1
            currents(1,i) = (X(i)-p.v1)/p.REmitter; 
        else
            currents(1,i) = (X(i)-p.v2)/p.RCollector;
        end
    end
    
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
%     disp(currents)
    G = digraph(sc,tc,currents);
    Gc = flipedge(G,1:2:totalnum);
%     Gc = addnode(Gc,totalnum)
    yc = [y y];
    xc = [x xc];
    
    LWidths = 5*abs(Gc.Edges.Weight/max(Gc.Edges.Weight))+0.5;
    
%     nodelabels = string(1:1:totalnum);
%     emptylabels = repmat("",1,totalnum);
    pg = plot(Gc,'XData',xc,'YData',yc,'NodeLabel',[],'LineWidth',LWidths, 'EdgeColor','r');
    % Subtract steady state for better visualization
%     Xnew = zeros(length(X),1);
%     Xnew(1:2:end) = X(1:2:end) - U(1);
%     Xnew(2:2:end) = X(2:2:end) - U(2);
    Xnew = X;
%     pg = plot(Gc,'XData',xc,'YData',yc,'NodeLabel',[]);
    for i = 1:totalnum     
        highlight(pg,i,'MarkerSize',(abs(Xnew(i))+1)*1);
        highlight(pg,i+totalnum,'MarkerSize',1,'NodeColor','r');    
    end
    trix1 = [-sqrt(3)/2 sqrt(3)/2 -sqrt(3)/2 -sqrt(3)/2]*0.15*0.6;
    triy1 = [1 0 -1 1]*0.15;
    trix2 = gap+[sqrt(3)/2 -sqrt(3)/2 sqrt(3)/2 sqrt(3)/2]*0.15*0.6;
    triy2 = [1 0 -1 1]*0.15;
    hold on
    for j = 1:p.col
        plot([(j-1)-biasgap (j-1)-biasgap],[0 -p.row], 'r', 'LineWidth',1)
        plot([(j-1)+gap+biasgap (j-1)+gap+biasgap],[0 -p.row], 'r', 'LineWidth',1)
        for i = 1:p.row
           plot(trix1+(j-1),triy1 - (i-1), 'k', 'LineWidth', 2)
           plot(trix2+(j-1),triy2 - (i-1), 'k', 'LineWidth', 2)
        end
    end
 
    drawnow
    hold off 



