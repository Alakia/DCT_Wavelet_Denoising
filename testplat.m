clc;clear all;
image_index       = [101];
image_index       = [num2str(image_index),'.png']; 
[O_image]         = double((imread(num2str(image_index))));
L                 = 1;
i                 = sqrt(-1);
s                 = zeros(size(O_image));
for k = 1:L
    s             = s + abs(randn(size(O_image)) + 1i * randn(size(O_image))).^2 / 2;
end
N_image           = O_image .* sqrt(s / L);%%%%%%%ÔëÉùÍ¼
mean_noise        = gamma(L+0.5)*(1/L)^(1/2)/gamma(L);%%%%%%%ÔëÉùµÄ¾ùÖµ
SN_image          = N_image/mean_noise;
[m,n]             = size(N_image);
p_image           = padarray(N_image,[3 3],'symmetric');
e_image           = zeros(m,n);
for  k = 1:m
    for j = 1:n
        e_image(k,j) = mean(mean( p_image(k:k+6,j:j+6)));
    end
end
matrix_sigma      = ((e_image/mean_noise).*(4/pi-1)^0.5)./L.^0.5;
matsig_col        = im2col(matrix_sigma,[8 8],'distinct');
mask_col          = zeros(size(matsig_col));
var_col           = var(matsig_col,0,1);
var_colall        = sum(var_col,2)/size(var_col,2);
for i = 1:size(var_col,2)
    if var_col(:,i) <= var_colall
        mask_col(:,i) = 1;
    end
end
maskfinal         = col2im(mask_col,[64 1],[256 256],'distinct');
imshow(maskfinal);