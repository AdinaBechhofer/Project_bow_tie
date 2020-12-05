totalnum = 2*p.row*p.col;
gap = 0.25;
biasgap = 0.3;
trix1 = [-sqrt(3)/2 sqrt(3)/2 -sqrt(3)/2 -sqrt(3)/2]*0.15*0.6;
triy1 = [1 0 -1 1]*0.15;
trix2 = gap+[sqrt(3)/2 -sqrt(3)/2 sqrt(3)/2 sqrt(3)/2]*0.15*0.6;
triy2 = [1 0 -1 1]*0.15;
figure;
hold on

for j = 1:p.col
    plot([(j-1)-biasgap (j-1)-biasgap],[0 -p.row], 'r', 'LineWidth',1)
    plot([(j-1)+gap+biasgap (j-1)+gap+biasgap],[0 -p.row], 'r', 'LineWidth',1)
    for i = 1:p.row
%        plot(trix1+(j-1),triy1 - (i-1), 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
%        plot(trix2+(j-1),triy2 - (i-1), 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
       fill(trix1+(j-1),triy1 - (i-1), 'y')
       fill(trix2+(j-1),triy2 - (i-1), 'y')
       
    end
end
set(gca,'visible','off')
x0=10;
y0=10;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'bowtie49yellow.jpg') ;
% close(gcf) ;    
hold off;
% figure;
% array = imread('bowtie7x7.jpg');
% imshow(array);
% hold on;

