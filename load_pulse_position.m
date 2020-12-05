function [b,demopt] = load_pulse_position(p,optionflag)
    
    unitb = [1, 0, 0, 0; 0, 1, 0, 0];
    b = repmat(unitb,p.NumBowties,1); % bias always on

if strcmp(optionflag, 'centre-only')
    if mod(p.NumBowties,2) == 0
        error("use odd no. of bowties")
    else
        % turn on centre bowtie
        b(p.NumBowties, 3) = 1;
        b(p.NumBowties+1, 4) = 1;
    end
    demopt = 150:1:250;

elseif strcmp(optionflag,'centre-neighbours')
    if mod(p.NumBowties,2) == 0
        error("use odd no. of bowties")
    else
        b(p.NumBowties, 3) = 1;
        b(p.NumBowties+1, 4) = 1;
        % neighbours get half
        b(p.NumBowties-2, 3) = 0.5;
        b(p.NumBowties+1-2, 4) = 0.5;
        b(p.NumBowties+2, 3) = 0.5;
        b(p.NumBowties+1+2, 4) = 0.5;
        b(p.NumBowties-p.row*2, 3) = 0.5;
        b(p.NumBowties+1-p.row*2, 4) = 0.5;
        b(p.NumBowties-p.row*2-2, 3) = 0.25;
        b(p.NumBowties+1-p.row*2-2, 4) = 0.25;
        b(p.NumBowties-p.row*2+2, 3) = 0.25;
        b(p.NumBowties+1-p.row*2+2, 4) = 0.25;
        b(p.NumBowties+p.row*2, 3) = 0.5;
        b(p.NumBowties+1+p.row*2, 4) = 0.5;
        b(p.NumBowties+p.row*2-2, 3) = 0.25;
        b(p.NumBowties+1+p.row*2-2, 4) = 0.25;
        b(p.NumBowties+p.row*2+2, 3) = 0.25;
        b(p.NumBowties+1+p.row*2+2, 4) = 0.25;
    end
    demopt = 150:1:250;
elseif strcmp(optionflag, 'left-side')
    for ii = 1:2*p.row
        if mod(ii,2) == 1
            b(ii,3) = 1;
        else
            b(ii,4) = 1;
        end
    end
    demopt = 100:1:200;
    
elseif strcmp(optionflag,'left-corner')
    b(1,3) = 1;
    b(2,4) = 1;
    demopt = 150:1:250;
end
