function [A]=OMPerrNONHtest(D,X,BETA,errorGoal,LL)
%=============================================

P=size(X,2);
K=size(D,2);
A=sparse(K,P); 

%��P�����ݵ�ÿһ�н���OMP���㣬�õ����Ӧ��ϵ��
for i=1:P
   
    x=X(:,i);
    beta=BETA(:,i);
    
    %1 �õ����beta�ĸ�������
    x_beta=beta.*x;
    
    %2 �õ����beta���ֵ�
    D_beta=zeros(size(D));
    for j=1:K
        D_beta(:,j)=D(:,j).*beta;
    end
    %3 �Ը������ݺ��ֵ����һ��OMP �õ��������ݶ�Ӧ�ֵ��ϵ��
    [a,indx]=omp_single(D_beta,x_beta,errorGoal,LL);
    
    %4 ����ʾ�������ݵ�ϵ������ϵ������ 
    if (~isempty(indx))
           A(indx,i)=a;
    end
end
end
