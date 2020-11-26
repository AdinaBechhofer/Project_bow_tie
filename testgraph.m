p.row = 5;
p.col = 5;
u.vEmitter = 5;
u.vCollector = 0;
p.RCollector = 1;
p.REmitter = 1;
X = rand(p.row*p.col*2,10);
figure;
% hold on
for t = 1:1:10
    VisualizeNetwork(X(:,t),p,u)
end
