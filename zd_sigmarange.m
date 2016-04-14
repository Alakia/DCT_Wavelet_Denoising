% function sigma_range=sigmarange(input)

clear all;
    L=8;
    syms v ;    
    sig = int((2*L^L* v^(2*L-1)/gamma(L))*exp(-(L*(v^2))),v);
    sig2 = int((2*L^L* v^(2*L)/gamma(L))*exp(-(L*(v^2))),v);
%     sig = int(L^L* v^(L-1)/gamma(L)*exp(-(L*v)),v);
%     sig2 = int(L^L* v^(L)/(gamma(L))*exp(-(L*(v))),v);

sigma0=0.9;
%    
        lamda1 = 0.4;
        lamda2 = 2.0;
%         lamda1 = 0.8;
%         lamda2 = 1.4;

   while 1
        tempA= 0;
        tempB=0;
        flag1 = 1;

        while(flag1)
            sigma =double(subs(sig,v,lamda2))-double(subs(sig,v,lamda1));
            if sigma > sigma0
                lamda2 = lamda2-0.001;
                tempA=1;
            elseif sigma < sigma0
                lamda2 = lamda2+0.001;
                tempB=1;
            else 
                flag1 = 0;
            end
            if tempA==1 && tempB==1
                flag1 = 0;
            end
        end
        tempA= 0;
        tempB=0;

        flag2 = 1;
        flag1= 1;
%         lamda1= 0.654 ;     
%         lamda2 = 1.40;
%         sigma0 = 0.5;
        a=0;
        while(flag2)
             ave =(double(subs(sig2,v,lamda2))-double(subs(sig2,v,lamda1))) /sigma0;
%                 ave = double(int((2*L^L* v^(2*L)/gamma(L)/ave0^L)*exp(-(L*((v/ave0)^2))),v,lamda1,lamda2))/sigma0(k);
             if ave> 1.01
                lamda1 = lamda1-0.01;
                tempA=1;a=a+1;
                break
            elseif ave < 0.99
                lamda1 = lamda1+0.01;
                tempB=1;a=a+1;
                break
            else
                 flag2 = 0;
                if a==0
                    flag1=0;
%                     sigma_range(1,n)=lamda1;
%                     sigma_range(2,n)=lamda2;
                    break
                end
            end

        end
       if flag1==0 && flag2==0
           break
       end
   end
% end
% sigma_range
lamda1
lamda2
