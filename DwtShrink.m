function shrinkage = DwtShrink(input,blocksize,rate);
[num ind]          = sort(input(:),1,'ascend');
% mask               = zeros(blocksize,blocksize);
for i = 1:numel(input)*2/3
    input(ind(i))  = 0;
end
shrinkage          = input;
% shrinkage          = input.*mask;
end