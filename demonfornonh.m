% addpath('ppbNakagami')
clear all
        bb  = 8; 
        K   = 256;
        C   = 1.1;
       
for image_index = [   1502 1503 1504 1505  ]
    image_index = [num2str(image_index),'.png']; 
        for L = [  1  ]; %%%%视数
            psnrtemp           = zeros(2,1);IFinal = [];
            dictemp = zeros(64*256,15);coeftemp = zeros(256*(256-8+1)*(256-8+1),1);
            psnrper = zeros(1,1);mseper = zeros(1,1);ssimper = zeros(1,1);imageper25 = [];mseblk = zeros(1,15625);
            recount = zeros(1,1);
            for iteration = 1:1
                disp(iteration);
                save('index.mat','image_index','L','iteration','bb','K','C','psnrtemp','dictemp','coeftemp','recount','psnrper','mseper','ssimper','imageper25','mseblk');
                clear all;
                load  index.mat;
                [O_image]      = double((imread(num2str(image_index))));%%%%读一幅图
                slidingDis     = 2;   
%                 randn('seed', 0); 
%                 i              = sqrt(-1);
%                 s              = zeros(size(O_image));
%                     for k = 1:L
%                         s      = s + abs(randn(size(O_image)) + 1i * randn(size(O_image))).^2 / 2;
%                     end
%                 N_image        = O_image .* sqrt(s / L);%%%%%%%噪声图
%                 N_image        = O_image .* (s / L);
N_image = O_image;
                mean_noise     = gamma(L+0.5)*(1/L)^(1/2)/gamma(L);%%%%%%%噪声的均值
%                 mean_noise = 1;
                SN_image       = N_image/mean_noise;%%%为什么只对这里进行噪声均值归一化
                PSNRIn         = 20*log10(255/sqrt(mean((SN_image(:)-O_image(:)).^2)));
% PSNRIn = 0;
                %%%%%%%%%%求权值矩阵
                [m,n]          = size(N_image);%%%这里为什么不用SN_image
                p_image        = padarray(N_image,[3 3],'symmetric');
                e_image        = zeros(m,n);
                    for  k = 1:m
                        for j = 1:n
                            e_image(k,j) = mean(mean( p_image(k:k+6,j:j+6)));
                        end
                    end
                matrix_sigma   = ((e_image/mean_noise).*(4/pi-1)^0.5)./L.^0.5;
%                 matrix_sigma   = (e_image/mean_noise)./sqrt(L);
%                 tic
                [mseblk0 timeout ksvd_step_output psnr_step D_step IOut Dictionary Coefs] = denoiseImageKSVDNONH(O_image,SN_image,matrix_sigma,bb,K,C,slidingDis,L);
                time              = timeout;
                PSNROut           = 20*log10(255/sqrt(mean((IOut(:)-O_image(:)).^2)));MSE = sqrt(mean((IOut(:)-O_image(:)).^2));
                 K3 = [0.01 0.03];
                window3 = ones(8);
                L3 =max(IOut(:));
                [mssim, ssim_map] = ssim(O_image,IOut, K3, window3, L3);
                ratioimg = IOut./SN_image;
%                 mseblk(iteration,:) = mseblk0;
                
% for number = 1:15
%     psnrper(iteration,number) = 20*log10(255/sqrt(mean((ksvd_step_output(:,number)-O_image(:)).^2)));
%     mseper(iteration,number) = sqrt(mean((ksvd_step_output(:,number)-O_image(:)).^2));
%     [mssim, ssim_map] = ssim(O_image,IOut);
%     ssimper(iteration,number) = mssim;
% end
% imageper25(:,:,iteration) = ksvd_step_output;

                
                
%                 for itn = 1:16
%                     ksvdtemp = ksvd_step_output(:,itn);
%                 recount(itn,iteration) = 20*log10(255/sqrt(mean((ksvdtemp(:)-O_image(:)).^2)));
%                 end
% PSNROut = 0;
                psnrtemp(1,iteration) = PSNROut;
                psnrtemp(2,iteration) = time;
                psnrtemp(3,iteration)=mssim;
%                 outputtemp(:,iteration) = IOut(:);
                dictemp(:,iteration) = Dictionary(:);
%                 coeftemp(:,iteration) = Coefs(:);
                [tp1 tp2] = size(Coefs);
%                 if PSNROut >= max(psnrtemp(1,:))
%                     save('7-25maxwaveletzero.mat','IOut');
%                     IFinal = IOut;
%                 end
%                 if PSNROut <= max(psnrtemp(1,:))
%                     save('7-25minwaveletzero.mat','IOut');
%                 end
                disp(PSNRIn);disp(PSNROut);disp(time);
            end
            PSNRFinal = mean(psnrtemp(1,:));
            TFinal    = mean(psnrtemp(2,:));
%             tempp1 = psnrtemp(1,:);
%             [num bes] = max(tempp1);
%             bestoutput = reshape(outputtemp(:,bes),[256 256]);
%             bestdic = reshape(dictemp(:,bes),[64 256]);
%             bestcoef = reshape(coeftemp(:,bes),[tp1 tp2]);
%             [coefplots] = coefplot(bestcoef);
            save(['F:\MatlabWork\new\code\2013.7.19-二阶段缩减效果实验-小波3\data\',num2str(image_index),'\+mse1-18-newim2col-collectDmse-t8d32-tssc-trainper2-sqrt2-iter1','--L=',num2str(L),'-Psnr',num2str(PSNRFinal),'-mse',num2str(MSE),'-ssim',num2str(mssim),'-trainingTime=',num2str(TFinal),'.mat']','ksvd_step_output','psnr_step','D_step','psnrtemp','IOut','PSNRFinal','TFinal','Dictionary','Coefs','dictemp','ksvd_step_output','psnr_step','D_step','recount','psnrper','imageper25','ssimper','mseblk','ratioimg');
        end
end

% end
