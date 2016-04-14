function  [ksvd_step_output psnr_step D_step Dictionary] = KSVDNONH(Data,DCT,BETA,errorGoal,L,Image,BETA_sig,sigma,matrix_sigma)
% =========================================================================
sigma_ap       =    sqrt((4/pi-1)/L);
% sigma_ap = 1/sqrt(L);
mean_noise     = gamma(L+0.5)*(1/L)^(1/2)/gamma(L);
% mean_noise = 1;
% sigma_ap       = (sigma_ap/mean_noise);
% sigma_ap = sigma_ap*0.5
Dictionary                              = DCT;
%1 对初始的字典进行归一化处理
Dictionary                              = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
Dictionary                              = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.
iter                                    = 15;
ksvd_step_output                        = zeros(256*256,iter);
psnr_step                               = zeros(1,iter);
D_step                                  = zeros(64*256,iter);
% IOutstep                                = DenoisingByStep(Image,BETA_sig,Dictionary,sigma);
% ksvd_step_output(:,1)                   = IOutstep(:);
% psnr_step(1)                            = 20*log10(255/sqrt(mean((IOutstep(:)-Image(:)).^2)));
% D_step(:,1)                             = Dictionary(:);
[Tfor1,Tfor2,Tfor3,Tinv1,Tinv2,Tinv3]   = nonsumple_wavelet_matrix(8);
mean_noise     = gamma(L+0.5)*(1/L)^(1/2)/gamma(L);
C = 1.1;
% sigmau = ((4/pi)-1)/L; sigmau = sigmau/(1+sigmau); r = 0.6;
%2 对字典进行20次的更新
for iterNum = 1:iter
    disp(iterNum)
    % 2.1 利用字典，进行一次OMP求解，得到系数矩阵coefmatrix 每列为一个像素点patch对应于字典每一列的系数
    CoefMatrix                          = OMPerrNONHtest(Dictionary,Data,BETA,errorGoal,8);
   
%     CoefMatrix                     = OMPnew(Dictionary,Data,8,BETA);
% TSSC二次追踪     
% CoefMatrix_firstdenoiseblk = OMP(DCT,Data-Dictionary* CoefMatrix,10);
%     CoefMatrix_firstdenoiseblk          = OMPerrNONHtest(DCT,Data-Dictionary* CoefMatrix,BETA,0.85*errorGoal,53);
%     DA2                                 = DCT*CoefMatrix_firstdenoiseblk;

    DA1 = Dictionary*CoefMatrix;
%     OriRes                              = Data-DA1;
    Residual                            = (Data-Dictionary*CoefMatrix).*BETA;
%     %% column directly
% %     STC = 3*sigma_ap;
% %     for i3 = 1:size(Residual,2)
% %         dwttemp = Residual(:,i3);
% %         [DC DL] = wavedec(dwttemp,3,'db8');
% %         DC(DL(1)+1:end,1) = SoftThreshold(DC(DL(1)+1:end,1),STC);
% %         Residual(:,i3) = waverec(DC,DL,'db8');
% %     end
% %% 
    for i = 1:size(Data,2)
%         %%
%     
%   %%
        Dwttemp                         = reshape(Residual(:,i),[8 8]);
        cof1                            = 2*(Tfor1*Dwttemp*Tinv1); %第一层小波分解系数
        cof2                            = 2*(Tfor2*cof1(1:8,1:8)*Tinv2);%第二层小波分解系数
        cof3                            = 2*(Tfor3*cof2(1:8,1:8)*Tinv3); %第三层小波分解系数
% %   LL3=cof2(1:8,1:8);
  LL3 = cof3(1:8,1:8);
  LH1=cof1(9:16,1:8);HL1=cof1(1:8,9:16);HH1=cof1(9:16,9:16);
  LH2=cof2(9:16,1:8);HL2=cof2(1:8,9:16);HH2=cof2(9:16,9:16);
  LH3=cof3(9:16,1:8);HL3=cof3(1:8,9:16);HH3=cof3(9:16,9:16);
  %% SoftThresholding
%   Med1 = median(HH1(:))/0.6745;
%   Med2 = median(HH2(:))/0.6745;
%   Med3 = median(HH3(:))/0.6745;
% %   HT1  = 3*Med1;
thd = sqrt(2*log10(64))*sigma_ap;
  HT1 = thd;
  HT2 = thd;
  HT3 = thd;
% %   HT1 = thselect(Dwttemp,'minimaxi');
% %   HT2  = 3*Med1;
% %   HT3  = 3*Med1;
  LH1 = SoftThreshold(LH1,HT1);HL1 = SoftThreshold(HL1,HT1);HH1 = SoftThreshold(HH1,HT1);
  LH2 = SoftThreshold(LH2,HT2);HL2 = SoftThreshold(HL2,HT2);HH2 = SoftThreshold(HH2,HT2);
  LH3 = SoftThreshold(LH3,HT3);HL3 = SoftThreshold(HL3,HT3);HH3 = SoftThreshold(HH3,HT3);
  LL3 = SoftThreshold(LL3,HT1);
% %   cof2(1:8,1:8)=LL3;
  cof3(1:8,1:8)=LL3;
  cof1(9:16,1:8)=LH1;cof1(1:8,9:16)=HL1;cof1(9:16,9:16)=HH1;
  cof2(9:16,1:8)=LH2;cof2(1:8,9:16)=HL2;cof2(9:16,9:16)=HH2;
  cof3(9:16,1:8)=LH3;cof3(1:8,9:16)=HL3;cof3(9:16,9:16)=HH3;

  %%  Shrinking
%   DZ=(LL3.^2+LH1.^2+LH2.^2+LH3.^2+HL1.^2+HL2.^2+HL3.^2+HH1.^2+HH2.^2+HH3.^2)/10;
%   DV = ((Residual(:,i)).^2)*sigmau;DV = reshape(DV,[8 8]);
%   shrinkcoef=r*max(((DZ-DV)./DZ),0);
% %   shrinkcoef = 1/3;
%   cof1(9:16,1:8)=cof1(9:16,1:8).*shrinkcoef;cof2(9:16,1:8)=cof2(9:16,1:8).*shrinkcoef;cof3(9:16,1:8)=cof3(9:16,1:8).*shrinkcoef;
%   cof1(1:8,9:16)=cof1(1:8,9:16).*shrinkcoef;cof2(1:8,9:16)=cof2(1:8,9:16).*shrinkcoef;cof3(1:8,9:16)=cof3(1:8,9:16).*shrinkcoef;
%   cof1(9:16,9:16)=cof1(9:16,9:16).*shrinkcoef;cof2(9:16,9:16)=cof2(9:16,9:16).*shrinkcoef;cof3(9:16,9:16)=cof3(9:16,9:16).*shrinkcoef;
%   cof3(1:8,1:8)=cof3(1:8,1:8).*shrinkcoef;
%             LL3                         = cof3(1:8,1:8);
%             cof1                        = DwtShrink(cof1,8,3);
%             cof2                        = DwtShrink(cof2,8,3);
%             cof3                        = DwtShrink(cof3,8,3);
%             cof3(1:8,1:8)               = LL3;
%%
        cof2(1:8,1:8)                   =    (Tinv3*cof3*Tfor3)/2;
        cof1(1:8,1:8)                   =    (Tinv2*cof2*Tfor2)/2;
        inv_Block                       =    (Tinv1*cof1*Tfor1)/2;
        Residual(:,i)                   = reshape(inv_Block,[64 1]);
    end
    DA2                                 = Residual./BETA;
%             DA2 = [];
     %2.2 对字典进行逐列更新
%      disp('start update :')
    rPerm                               = randperm(size(Dictionary,2));
%     DA2=zeros(3,3);
    for j = rPerm
        [betterDictionaryElement,CoefMatrix] = I_findBetterDictionaryElement(Data,Dictionary,BETA,j,CoefMatrix, DA2 );
        Dictionary(:,j)                      = betterDictionaryElement;
    end  
%     IOutstep                                 = DenoisingByStep(Image,BETA_sig,Dictionary,sigma);
%     ksvd_step_output(:,iterNum)               = IOutstep(:);
%     psnr_step(iterNum)                        = 20*log10(255/sqrt(mean((IOutstep(:)-Image(:)).^2)));
%     D_step(:,iterNum)                         = Dictionary(:);
    
     %% update beta
%     if mod(iter,5) == 0
%     sigma                      = min(matrix_sigma(:));
%     errorGoal                  = sigma*C;
%     IOuttemp                   = DenoisingByStep(Image,BETA_sig,Dictionary,sigma);
%     [m,n]          = size(IOuttemp);%%%这里为什么不用SN_image
%     p_image        = padarray(IOuttemp,[3 3],'symmetric');
%     e_image        = zeros(m,n);
%     for  k = 1:m
%         for j = 1:n
%             e_image(k,j) = mean(mean( p_image(k:k+6,j:j+6)));
%         end
%     end
%     matrix_sigma   = ((e_image/mean_noise).*(4/pi-1)^0.5)./L.^0.5;
%     BETA_sig                   = sigma*(ones(size(matrix_sigma))./matrix_sigma);
%     [BETA]                    = my_im2col(BETA_sig,[8,8],1);
%     end
    %%
    
end
end



