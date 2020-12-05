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
        if mod(i,2)==1
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

    yc = [y y];
    xc = [x xc];
    
    LWidths = abs(Gc.Edges.Weight/max(Gc.Edges.Weight))+0.5;
    MSize = [X; zeros(length(X),1)];
    MSize(1:2:length(X))=abs(MSize(1:2:length(X))-5);
    MSize = abs(MSize) * 20 + 0.5;
    LWidths(LWidths>5)=5;
    MSize(MSize>10) = 10+2*log(MSize(MSize>10)-9); % reduce growing rate after 10
    
%     MSize(MSize>10) = 10
    
%     nodelabels = string(1:1:totalnum);
%     emptylabels = repmat("",1,totalnum);
    
    array = imread('bowtie49yellow.jpg');
    imagesc([-1.2 7.05], [-7.15 1.3], array);
    hold on
%     openfig('firstfig.fig','reuse')
    pg = plot(Gc,'XData',xc,'YData',yc,'NodeLabel',[],'LineWidth',LWidths, 'EdgeColor','r','MarkerSize',MSize);
    
    % Subtract steady state for better visualization
%     Xnew = zeros(length(X),1);
%     Xnew(1:2:end) = X(1:2:end) - U(1);
%     Xnew(2:2:end) = X(2:2:end) - U(2);

%     pg = plot(Gc,'XData',xc,'YData',yc,'NodeLabel',[]);
%     for i = 1:totalnum     
%         highlight(pg,i,'MarkerSize',(abs(Xnew(i))+1)*1);
%         highlight(pg,i+totalnum,'MarkerSize',1,'NodeColor','r');    
%     end
%     trix1 = [-sqrt(3)/2 sqrt(3)/2 -sqrt(3)/2 -sqrt(3)/2]*0.15*0.6;
%     triy1 = [1 0 -1 1]*0.15;
%     trix2 = gap+[sqrt(3)/2 -sqrt(3)/2 sqrt(3)/2 sqrt(3)/2]*0.15*0.6;
%     triy2 = [1 0 -1 1]*0.15;
%     hold on
%     for j = 1:p.col
%         plot([(j-1)-biasgap (j-1)-biasgap],[0 -p.row], 'r', 'LineWidth',1)
%         plot([(j-1)+gap+biasgap (j-1)+gap+biasgap],[0 -p.row], 'r', 'LineWidth',1)
%         for i = 1:p.row
%            plot(trix1+(j-1),triy1 - (i-1), 'k', 'LineWidth', 2)
%            plot(trix2+(j-1),triy2 - (i-1), 'k', 'LineWidth', 2)
%         end
%     end
    
    
    hold off
    x0=10;
    y0=10;
    width=500;
    height=400;
%     set(gcf,'position',[x0,y0,width,height])
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'visible','off')
    
   
   



