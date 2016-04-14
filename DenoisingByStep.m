function IOut = DenoisingByStep(Data,BETA_sig,Dictionary,sigma)
  slidingDis                = 1;
  C                         = 1.1;
  bb                        = 8;
  reduceDC                  = 0;
  [NN1,NN2]                 = size(Data);
  errT                      = sigma*C;
  [blocks,idx]              = my_im2col(Data,[bb,bb],slidingDis);%%%%%%%%%%%%%%�����ص�ȡ������л�
  [blocksbeta,idxbeta]      = my_im2col(BETA_sig,[bb,bb],slidingDis);%%%%%%�����ص�ķ���ȡ������л�
  if (reduceDC)
      vecOfMeans            = mean(blocks);
      blocks                = blocks- repmat(vecOfMeans,size(blocks,1),1);
  end
  Coefs                     = OMPerrNONHtest(Dictionary,blocks,blocksbeta,errT,256);%%%%%%%%%���Ӧ��ϵ��
%   Coefs                     = OMPnew(Dictionary,blocks,8,blocksbeta);
  if (reduceDC)
      blocks                = Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
  else
      blocks                = Dictionary*Coefs ;%%%%%%%%%%%%%�õ�ȥ���Ŀ�
  end
  count                     = 1;
  Weight                    = zeros(NN1,NN2);
  IMout                     = zeros(NN1,NN2);
  [rows,cols]               = ind2sub(size(Data)-bb+1,idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ÿ�����ؽ���ƽ��
  for i = 1:length(cols)
      col                   = cols(i); row = rows(i);        
      block                 = reshape(blocks(:,count),[bb,bb]);
      IMout(row:row+bb-1,col:col+bb-1)  = IMout(row:row+bb-1,col:col+bb-1)+block;
      Weight(row:row+bb-1,col:col+bb-1) = Weight(row:row+bb-1,col:col+bb-1)+ones(bb);  
      count                 = count+1;
  end;
  IOut                      = IMout./Weight;%%%%%%%%%ȥ����ͼ
end
