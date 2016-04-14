function [A]=OMPnew(D,X,L,BETA); 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  L - the maximal number of coefficient for representation
%                  of each signal.
% output arguments: A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);
[n,K]=size(D);
for k=1:1:P
    beta = BETA(:,k);
    a=[];
    x=X(:,k).*beta;
    D = D.*repmat(beta,1,K);
    residual=x;
    indx=zeros(L,1);
    for j=1:1:L,
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
    end;
    temp=zeros(K,1);
    temp(indx)=a;
    A(:,k)=sparse(temp);
end;
return;
