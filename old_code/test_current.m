p.NumBowties = 100;
p.REmitter = 1; % 1 ohm 
p.RCollector = 1; % 1 ohm 
p.Area = 5e-18; % 5 nm^2
p.Beta = 25; % enhancement factor 
p.Distance = 10; % 10 nm
p.workFunction = 5.1; % work function of gold 
p.CemitterCollector = 0.1; % nano farad
p.Cparasitic = 0.08; % 0.1 nano farad
p.Radius = 10; % 1 nm?
p.taby = csvread('rspa20140811supp3.csv');    
p.Col = 10;
p.Row = 10;
p.Cwire = 0.08*1e5;
p.Rwire = 1*1e-4;
u.Wire1Bias = 10/p.Rwire;
u.Wire2Bias = 0;
u.jnano = 1;

p.Area*Jnano(p.workFunction,p.Beta*5/p.Distance,p.Radius, p.taby)