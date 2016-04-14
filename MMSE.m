function   [x,flag] = MMSE( z,orseed,yita)%%%%%%%%N*1����
%%%%%%%%%%%%%%%%%%%%zͼ������ת����n*1���У�seed����Ҫ��������ĵ㣬yita����
%for one-look amplitude data, ��v=0.5227 
%for one-look intensity, ��v=1
% [M,N]=size(z);

% varx=max(0,(var(z,0,1)-mean(z,1)^2*yita^2)/(1+yita^2));
% if var(z,0,1)==0
%     x=mean(z);
% else
%     b=varx/var(z,0,1);
%     if b>=1
%         b=1;
%     elseif b<=0
%         b=0;
%     end
%     x=(1-b)*mean(z)+b*orseed;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���˵�MMSE    
flag=0;         
 selectnumber=size(z);           
if (selectnumber(1)>0)
             
              
              meanW2 = mean(z,1);


                if(selectnumber(1) ~= 1)
                     varW2 = var(z,0,1);
                else
                    varW2 = 0;         
                end

                varx2 = max(0,(varW2-(meanW2*yita)^2)/(1+yita^2));
                 b2 = varx2/varW2;
                 x = (1-b2)*meanW2+b2*orseed;      
                 if (varW2==0)
                    x=meanW2; 
                end
 else
              x=orseed;
              flag=1;
             
              
end

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

