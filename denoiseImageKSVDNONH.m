function [mseblk timeout ksvd_step_output psnr_step D_step IOut Dictionary Coefs] = denoiseImageKSVDNONH(O_image,Image,matrix_sigma,bb,K,C,slidingDis,L)
%==========================================================================
[NN1,NN2]                    = size(Image);
%1 �õ�beta����ÿһ��Ϊһ�����ص�path��beta
  %1.1�ҳ���С�ı�׼��sigma
  sigma                      = min(matrix_sigma(:));
  errorGoal                  = sigma*C;
  %1.2 ��sigma��ÿ�����ص�ķ�����й�һ��
  BETA_sig                   = sigma*(ones(size(matrix_sigma))./matrix_sigma);
  %1.3 �����ص�ķ���ȡ������л�
  %2 �����ص�ȡ������л����õ����ݾ���blkmatrix,ÿһ��Ϊһ�����ص��path
  [BETA]                    = my_im2col(BETA_sig,[bb,bb],slidingDis);
  [blkMatrix]               = my_im2col(Image,[bb,bb],slidingDis);
  
%3 �õ���ʼ��DCT�����ֵ�
  %3.1 �õ���ʼDCT
  Pn                        = ceil(sqrt(K));
  DCT                       = zeros(bb,Pn);
  for k = 0:1:Pn-1,
      V                     = cos([0:1:bb-1]'*k*pi/Pn);
      if k>0
          V                 = V-mean(V); 
      end;
      DCT(:,k+1)            = V/norm(V);
  end;
  DCT                       = kron(DCT,DCT);
%3.2 ��DCT���о�ֵ����
  reduceDC = 0;
  if (reduceDC)
      vecOfMeans            = mean(blkMatrix);
      blkMatrix             = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
  end
% 4 ���ֵ�DCT����ѵ��
tic
  [ksvd_step_output psnr_step D_step Dictionary]              = KSVDNONH(blkMatrix,DCT,BETA,errorGoal,L,Image,BETA_sig,sigma,matrix_sigma);
  timeout = toc;
  
  DCT                       = Dictionary;
  disp('finished Trainning dictionary');
% ksvd_step_output = [];psnr_step = [];D_step = [];
  Dictionary                = DCT;
%5 ��ѵ���õ��ֵ�DCT����ȥ��

[OblkMatrix]               = my_im2col(O_image,[bb,bb],slidingDis);
OCoefs                     = OMPerrNONHtest(DCT,blkMatrix,BETA,errorGoal,32);%%%%%%%%%���Ӧ��ϵ��
mseblk0 = sqrt((DCT*OCoefs-OblkMatrix).^2);mseblk = mean(mseblk0);

  slidingDis                = 1;
  errT                      = sigma*C;
     [blocks,idx]              = my_im2col(Image,[bb,bb],slidingDis);%%%%%%%%%%%%%%�����ص�ȡ������л�
%     [O_blocks,idx]              = my_im2col(O_image,[bb,bb],slidingDis);%%%%%%%%%%%%%%�����ص�ȡ������л�
  [blocksbeta,idxbeta]      = my_im2col(BETA_sig,[bb,bb],slidingDis);%%%%%%�����ص�ķ���ȡ������л�
  if (reduceDC)
      vecOfMeans            = mean(blocks);
      blocks                = blocks- repmat(vecOfMeans,size(blocks,1),1);
  end
  Coefs                     = OMPerrNONHtest(DCT,blocks,blocksbeta,errT,32);%%%%%%%%%���Ӧ��ϵ��
%   Coefs                     = OMPnew(DCT,blocks,8,blocksbeta);
  if (reduceDC)
      blocks                = DCT*Coefs + ones(size(blocks,1),1) * vecOfMeans;
  else
      blocks                = DCT*Coefs ;%%%%%%%%%%%%%�õ�ȥ���Ŀ�
  end
%   mseblk = sqrt((blocks-O_blocks).^2);
%   mseblk = mean(mseblk);
  count                     = 1;
  Weight                    = zeros(NN1,NN2);
  IMout                     = zeros(NN1,NN2);
  [rows,cols]               = ind2sub(size(Image)-bb+1,idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ÿ�����ؽ���ƽ��
  for i = 1:length(cols)
      col                   = cols(i); row = rows(i);        
      block                 = reshape(blocks(:,count),[bb,bb]);
      IMout(row:row+bb-1,col:col+bb-1)  = IMout(row:row+bb-1,col:col+bb-1)+block;
      Weight(row:row+bb-1,col:col+bb-1) = Weight(row:row+bb-1,col:col+bb-1)+ones(bb);  
      count                 = count+1;
  end;
  IOut                      = IMout./Weight;%%%%%%%%%ȥ����ͼ




