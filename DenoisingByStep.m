function IOut = DenoisingByStep(Data,BETA_sig,Dictionary,sigma)
  slidingDis                = 1;
  C                         = 1.1;
  bb                        = 8;
  reduceDC                  = 0;
  [NN1,NN2]                 = size(Data);
  errT                      = sigma*C;
  [blocks,idx]              = my_im2col(Data,[bb,bb],slidingDis);%%%%%%%%%%%%%%对像素点取块进行列化
  [blocksbeta,idxbeta]      = my_im2col(BETA_sig,[bb,bb],slidingDis);%%%%%%对像素点的方差取块进行列化
  if (reduceDC)
      vecOfMeans            = mean(blocks);
      blocks                = blocks- repmat(vecOfMeans,size(blocks,1),1);
  end
  Coefs                     = OMPerrNONHtest(Dictionary,blocks,blocksbeta,errT,256);%%%%%%%%%求对应的系数
%   Coefs                     = OMPnew(Dictionary,blocks,8,blocksbeta);
  if (reduceDC)
      blocks                = Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
  else
      blocks                = Dictionary*Coefs ;%%%%%%%%%%%%%得到去噪后的块
  end
  count                     = 1;
  Weight                    = zeros(NN1,NN2);
  IMout                     = zeros(NN1,NN2);
  [rows,cols]               = ind2sub(size(Data)-bb+1,idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%对每个像素进行平均
  for i = 1:length(cols)
      col                   = cols(i); row = rows(i);        
      block                 = reshape(blocks(:,count),[bb,bb]);
      IMout(row:row+bb-1,col:col+bb-1)  = IMout(row:row+bb-1,col:col+bb-1)+block;
      Weight(row:row+bb-1,col:col+bb-1) = Weight(row:row+bb-1,col:col+bb-1)+ones(bb);  
      count                 = count+1;
  end;
  IOut                      = IMout./Weight;%%%%%%%%%去噪后的图
end
