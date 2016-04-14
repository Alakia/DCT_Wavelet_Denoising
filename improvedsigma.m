for image_index = [  201 202 1301 1302 ]
    image_index = [num2str(image_index),'.png']; 
        for L = [  1 2 4 8  ]; %%%%视数
            for iteration = 1:1
            [O_image]      = double((imread(num2str(image_index))));%%%%读一幅图

randn('seed',0)
s = zeros(size(O_image));   
for k = 1:L
    s = s + abs(randn(size(O_image)) + 1i * randn(size(O_image))).^2 / 2;
end
N_image = O_image .* sqrt(s / L);%%%%%%%噪声图
% N_image = O_image;
nim = N_image;

if L == 1
tableA=[0.653997,1.40002,0.746019,0.208349;   %sigma =  0.5
            0.578998,1.50601,0.927012,0.255358;      %         0.6
            0.496999,1.63201,1.13501,0.305303;     %         0.7 
            0.403999,1.79501,1.39101,0.361078;     %         0.8
            0.286000,2.04301,1.75701,0.426375;     %         0.9
            0.203000,2.25998,2.05698,0.466398] ;    %         0.95    
                I_1 = 0.286;I_2 = 2.043;newstdnoise = 0.426375;
        %%4_look%%% I1    I2   I2-I1  yita      
    elseif L == 2
         tableA=[0.760,1.263,0.503,0.139021;   %sigma =  0.5
            0.705,1.332,0.627,0.169777;      %         0.6
            0.643,1.412,0.769,0.206675;     %         0.7 
            0.568,1.515,0.947,0.244576;     %         0.8
            0.467,1.673,1.206,0.291070;     %         0.9
            0.387,1.812,1.425,0.319955] ;    %         0.95
                I_1 = 0.467;I_2 = 1.673;newstdnoise = 0.291070;
    elseif L == 4
        tableA=[0.832,1.179,0.347,0.0894192;   %sigma =  0.5
            0.793,1.226,0.433,0.112018;      %         0.6
            0.747,1.279,0.532,0.139243;     %         0.7 
            0.691,1.347,0.656,0.167771;     %         0.8
            0.613,1.452,0.839,0.201036;     %         0.9
            0.548,1.543,0.995,0.222048] ;    %         0.95
        I_1 = 0.613;I_2 = 1.452;newstdnoise = 0.201036;
    elseif L == 8
        I_1 = 0.7300;I_2 = 1.3220;newstdnoise = 0.182048;
end
        

stdnoise       =    sqrt((4/pi-1)/L);
% stdnoise       =    sqrt(1/L);
tic
output = GRSL_Improvesigma(nim,1,3,stdnoise,newstdnoise,I_1,I_2,7);
time = toc;

PSNROut           = 20*log10(255/sqrt(mean((output(:)-O_image(:)).^2)))
MSE = sqrt(mean((output(:)-O_image(:)).^2))
    K3 = [0.01 0.03];
                window3 = ones(8);
                L3 =max(output(:));
                [mssim, ssim_map] = ssim(O_image,output, K3, window3, L3);
                mssim
            end
save(['F:\MatlabWork\new\code\2013.7.14-OMP部分代码效果调试 -2\data\',num2str(image_index),'\12-13-sigma-0-9-trueimprovedsigma-amplitude-randseed-iter1','--L=',num2str(L),'-Psnr',num2str(PSNROut),'-mse',num2str(MSE),'-ssim',num2str(mssim),'-Time=',num2str(time),'.mat']','output');
        end
end