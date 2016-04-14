function [output] = GRSL_Improvesigma(input0,f1,f2,stdnoise,newstdnoise,IorA1,IorA2,Tk)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  input: image to be filtered
 %  f1: 第一阶段估计均值用的窗口半径
 %  f2: 第二阶段用的窗口半径
 %  stdnoise: 修改前的噪声标准差
 %  newstdnoise: 修改后的噪声标准差
 %  IorA1: 表示估计的sigma范围的下限(幅度或强度)
 %  IorA2： 表示估计的sigma范围的上限（幅度或强度）
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [x y] = size(input0);
 
 % Replicate the boundaries of the input image
 input = padarray(input0,[f2 f2],'symmetric');
 
 output = zeros(x,y);
 n = (2*f2+1)^2;
 m = (2*f1+1)^2;
 mer = zeros(n,1);
 
 totalpixel = x*y;
 totalpixel02 = totalpixel*0.02;
 [N,posi] = hist(input0,100);
 numpixel = 0;
 for k1 = 100:-1:1
     for k2 = 1:y
         numpixel = numpixel+N(k1,k2);
     end
     if (numpixel > totalpixel02)
         break;
     end
 end
 
 Z98 = posi(k1-1);   
 

 
 for i=1:x
     for j=1:y
         
         i1 = i+f2;
         j1 = j+f2;
                                  
         block = input(i1-f1:i1+f1,j1-f1:j1+f1);
         W1 = block(:);
         
         K=0;
         if (input(i1,j1) >= Z98)
             for d=1:m
                 if (W1(d) >= Z98)
                     K = K+1;
                 end
             end
             if (K >= Tk)
                 output(i,j) = input(i1,j1);
                 continue;
             end
         end
         
             
         IA1 = IorA1;
         IA2 = IorA2;
         
         %MMSE
         meanW1 = mean(W1);                                        %均值
         varW1 = var(W1);                                          %方差
         varx = max(0,(varW1-(meanW1*stdnoise)^2)/(1+stdnoise^2));
         if(varW1~=0)
              b = varx/varW1;
              estimatorx1 = (1-b)*meanW1+b*input(i1,j1);
         else
              estimatorx1 = meanW1;
         end
         
         
         
%          直接对f1*f1窗求均值
%          estimatorx1 = mean(W1);
         

         IA1 = IA1*estimatorx1;
         IA2 = IA2*estimatorx1;
         
         
         block = input(i1-f2:i1+f2,j1-f2:j1+f2);
         W2 = block(:);
         
         meanW2 = 0;
         varW2 = 0;
         
         selectnumber = 0;
         for d = 1:n
             if (W2(d)>=IA1 && W2(d)<=IA2)
                 meanW2 = meanW2+W2(d);
                 selectnumber = selectnumber+1;
             end
         end
         
         if (selectnumber>0)
             
              meanW2 = meanW2/selectnumber;
         else
              output(i,j)=input(i1,j1);
              continue;
         end
         
         for d = 1:n
             if (W2(d)>=IA1&W2(d)<=IA2)
                 varW2 = varW2+(W2(d)-meanW2)^2;                
             end
         end
        if(selectnumber ~= 1)
             varW2 = varW2/(selectnumber-1);
        else
            varW2 = 0;
        end
                
        varx2 = max(0,(varW2-(meanW2*newstdnoise)^2)/(1+newstdnoise^2));
        if (varW2==0)
            output(i,j)=meanW2;
            continue;
        end
             
         b2 = varx2/varW2;
         output(i,j) = (1-b2)*meanW2+b2*input(i1,j1);
     end
 end
 
 