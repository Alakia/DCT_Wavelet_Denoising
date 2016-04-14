function [A]=OMPerrNONH(D,X,BETA,errorGoal); 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);
[n,K]=size(D);

D_beta=zeros(n,K);

maxNumCoef = n/2;
%%%%%%%%%
%%%%%%%%%%%
A = sparse(size(D,2),size(X,2));
for k=1:1:P,
    
    x=X(:,k);
    beta=BETA(:,k);
    b_n=n-length(find(beta==0));
    x_beta=beta.*x;
   
     D_beta=zeros(size(D));
    for i=1:K
        D_beta(:,i)=D(:,i).*beta;
    end
   
    E2 = errorGoal^2*b_n;
    residual= x_beta;
	indx = [];
	a = [];
	currResNorm2 = sum(residual.^2);
	j = 0;
    while currResNorm2>E2 && j < maxNumCoef,
		j = j+1;
        proj=D_beta'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D_beta(:,indx(1:j)))*x_beta;
        residual=x_beta_factor-D_beta(:,indx(1:j))*a;
		currResNorm2 = sum(residual.^2);
   end;
   if (length(indx)>0)
       A(indx,k)=a;
   end
end;
return;
