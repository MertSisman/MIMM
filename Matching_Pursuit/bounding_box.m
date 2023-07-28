function [lims] = bounding_box(Mask)
s = size(Mask);

AP1 = 1;
AP2 = s(1);
LR1 = 1;
LR2 = s(2);
SI1 = 1;
SI2 = s(3);

while sum(sum(Mask(AP1,:,:))) == 0
    AP1 = AP1 + 1;
end

if AP1 > 4 
    AP1  = AP1 - 3;
end

while sum(sum(Mask(AP2,:,:))) == 0
    AP2 = AP2 - 1;
end

if AP2 < s(1) - 3 
    AP2  = AP2 + 3;
end

while sum(sum(Mask(:,LR1,:))) == 0
    LR1 = LR1 + 1;
end

if LR1 > 4 
    LR1  = LR1 - 3;
end

while sum(sum(Mask(:,LR2,:))) == 0
    LR2 = LR2 - 1;
end

if LR2 < s(2) - 3 
    LR2  = LR2 + 3;
end


while sum(sum(Mask(:,:,SI1))) == 0
    SI1 = SI1 + 1;
end

if SI1 > 4 
    SI1  = SI1 - 3;
end

while sum(sum(Mask(:,:,SI2))) == 0
    SI2 = SI2 - 1;
end

if SI2 < s(3) - 3 
    SI2  = SI2 + 3;
end


lims = [{AP1:AP2}, {LR1:LR2}, {SI1:SI2}];

end