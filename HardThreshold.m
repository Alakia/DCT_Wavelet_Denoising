function HT = HardThreshold(input,H);
HT = input;
for i = 1:numel(input)
    if abs(input(i)) <= H
        HT(i) = 0;
    end
end
end