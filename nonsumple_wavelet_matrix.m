function [Tfor1,Tfor2,Tfor3,Tinv1,Tinv2,Tinv3] = nonsumple_wavelet_matrix(blocksize)
f1     =   2*blocksize;
[Tforward1, Tinverse1] = getTransfMatrix (blocksize, 'db8',1,0);%�����²�����һ��ֽ���󣬴�СΪ8*8
[Tforward2, Tinverse2] = getTransfMatrix (blocksize/2, 'db8',1,0);%�����²����ڶ���ֽ���󣬴�СΪ4*4
[Tforward3, Tinverse3] = getTransfMatrix (blocksize/4, 'db8',1,0);%�����²���������ֽ���󣬴�СΪ2*2

Tfor1  = zeros(f1,blocksize);%��һ��ֽ������ͼ�����ˣ����²�����һ��ֽ����ѭ��ƽ�Ƶõ�
Tfor2  = zeros(f1,blocksize);%�ڶ���ֽ������LL1��ˣ����²����ڶ���ֽ�����㡢ѭ��ƽ�Ƶõ�
Tfor3  = zeros(f1,blocksize);%������ֽ������LL2��ˣ����²���������ֽ�����㡢ѭ��ƽ�Ƶõ�
s1     = [1:2:blocksize];
s2     = [1:4:blocksize];
%���������²�������ͨ�������ѭ��ƽ�Ƶõ����²�������
Tfor1(1,:)          =           Tforward1(1,:);
for k = 2:blocksize
    Tfor1(k,:)      =            circshift(Tfor1(k-1,:),[0 1]);
end
Tfor1(blocksize+1,:)          =            Tforward1(blocksize/2+1,:);
for p = 10:f1
    Tfor1(p,:)      =            circshift(Tfor1(p-1,:),[0 1]);
end
Tinv1               =             pinv(Tfor1);

Tfor2(1,s1)         =            Tforward2(1,:);
for k = 2:blocksize
    Tfor2(k,:)      =            circshift(Tfor2(k-1,:),[0 1]);
end
Tfor2(blocksize+1,s1)         =            Tforward2(blocksize/4+1,:);
for p = 10:f1
    Tfor2(p,:)      =            circshift(Tfor2(p-1,:),[0 1]);
end
Tinv2               =            pinv(Tfor2);

Tfor3(1,s2)         =            Tforward3(1,:);
for k = 2:blocksize
    Tfor3(k,:)      =            circshift(Tfor3(k-1,:),[0 1]);
end
Tfor3(blocksize+1,s2)         =            Tforward3(blocksize/8+1,:);
for p = 10:f1
    Tfor3(p,:)      =            circshift(Tfor3(p-1,:),[0 1]);
end
Tinv3               =            pinv(Tfor3);