function [betterDictionaryElement,CoefMatrix] = I_findBetterDictionaryElement(Data,Dictionary,BETA,j,CoefMatrix,DA2)

%�ҳ������ֵ��еĸ��б����Ǽ�������
relevantDataIndices = find(CoefMatrix(j,:));
if (length(relevantDataIndices)<1)
    ERR=BETA.*(Data-Dictionary*CoefMatrix);
    sumerr=sum(ERR.^2);
    [indx]=find(sumerr==max(sumerr));
    pos=indx(1);
    betterDictionaryElement=(Data(:,pos)./BETA(:,pos))/norm((Data(:,pos)./BETA(:,pos)));
    
     
     CoefMatrix(j,:) = 0;

    return;
end
%% improved weighed lowrank
% [n,m]                                                          = size(BETA(:,relevantDataIndices));
% betterDictionaryElement                                        = zeros(n,1);
% betaVector                                                     = zeros(m,1);
% iter                                                           = 20;
% for i = 1:iter
% tmpCoefMatrix                                                  = CoefMatrix(:,relevantDataIndices); 
% tmpBETA                                                        = BETA(:,relevantDataIndices);
% tlabel                                                         = zeros(1,length(relevantDataIndices));
% dictemp                                                        = Dictionary(:,j);
% error                                                          = zeros(size(Data(:,relevantDataIndices)));
% for i3 = 1:length(relevantDataIndices)
%     tlabel(1,i3)                                               = dictemp'*diag(tmpBETA(:,i3))*dictemp;
%     error(:,i3)                                                = tmpBETA(:,i3).*(Data(:,i3)-Dictionary*tmpCoefMatrix(:,i3))+tlabel(1,i3)*dictemp*tmpCoefMatrix(j,i3);
% end
% ttemp                                                          = sum(BETA,2);
% t1                                                             = diag(ttemp)/length(relevantDataIndices);
% dictemp                                                        = Dictionary(:,j);
% t                                                              = dictemp'*t1*dictemp;
% save('test.mat','t');
% error                                                          = tmpBETA.*(Data(:,relevantDataIndices)-Dictionary*tmpCoefMatrix)+dictemp*tmpCoefMatrix(j,:);
% error                                                          = tmpBETA.*(Data(:,relevantDataIndices));
% for ii=1:size(t,2)
%     error(:,ii)=error(:,ii)+t(ii)*dictemp(:,ii)*tmpCoefMatrix(j,ii);
% end
% [betterDictionaryElement,singularValue,betaVector_t]           = svds(error,1);
% betaVector                                                     = betaVector_t*singularValue;
% CoefMatrix(j,relevantDataIndices)                              = betaVector';
% end


%%   weighted low rank and t method
%            tmpCoefMatrix = CoefMatrix(:,relevantDataIndices);
%            tmpCoefMatrix(j,:) = 0;
%            tmpBETA            = BETA(:,relevantDataIndices);
%            Error = (Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix); 
%            betaVector = CoefMatrix(j,relevantDataIndices);
%            DCTVector = Dictionary(:,j);
%            betterDictionaryElement_new = zeros(size(Dictionary(:,1)));
%            betaVector_new = zeros(size(tmpCoefMatrix(j,:)));
%            iter=10;
%            error = zeros(size(tmpBETA));
%            for i2 = 1:iter
%            error = zeros(size(errors));%%�����errorӦ���ǲв�����ݣ����òв����ԡ�
%            tao = zeros(1,size(betaVector,2));
%            tao = DCTVector'*diag(mean(tmpBETA,2))*DCTVector;
%            for i1 = 1:size(betaVector,2)
%                tao(:,i1) = DCTVector'*diag(tmpBETA(:,i1)/size(tmpBETA,2))*DCTVector;
%                error(:,i1) = tmpBETA(:,i1).*(Error(:,i1)-DCTVector*betaVector(:,i1))+tao(:,i1)*DCTVector*betaVector(1,i1);
%            end
%            error = tmpBETA.*Error+(ones(size(tmpBETA))-tmpBETA).*(DCTVector*betaVector);   %% weighted low rank          
%            error = tmpBETA.*(Error-DCTVector*betaVector)+tao*DCTVector*betaVector;   %% t method
           
%            [betterDictionaryElement_new,singularValue,betaVector_t]=svds(error,1);
%            betaVector_new=betaVector_t*singularValue;betaVector_new=betaVector_new';
%            DCTVector = betterDictionaryElement_new;betaVector = betaVector_new;
%            end
%            betterDictionaryElement = betterDictionaryElement_new;
%            CoefMatrix(j,relevantDataIndices) =betaVector_new;






% %�õ��⼸�����ݶ�Ӧ����Ӧ�ֵ��и��е�ϵ��
% tmpCoefMatrix = CoefMatrix(:,relevantDataIndices); 
% 
% %�����ȥ���⼸�����ݵı���еĸ��еĳɷֺ��⼸�����������ǵı��֮������
% tmpCoefMatrix(j,:) = 0;
Error = Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix; 
% Error= DA2(:,relevantDataIndices)+Dictionary(:,j)* CoefMatrix(j,relevantDataIndices); %%%Ӧ��û������
% %���µĸ��У�ʹ����ڼ�ȥ���еı�����С���������Ӧ��ϵ����
% % tmpBETA=BETA(:,relevantDataIndices);
[u,s,v]=svds(Error,1);
betterDictionaryElement=u;
          betaVector=s*v;
% % a=CoefMatrix(j,relevantDataIndices) ;
% % [betterDictionaryElement,betaVector]=weirank(tmpBETA,errors,errorGoal,DCT);
% 
% %����ϵ������
CoefMatrix(j,relevantDataIndices) =betaVector';% *signOfFirstElem
% disp(j)
end    
