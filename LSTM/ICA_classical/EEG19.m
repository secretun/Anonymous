clear;
close all;

%***********�˶���������2��������������16����***************
s=loadcnt('D1.cnt','data','int16');%����cnt�ļ�
N=2500;
% N=length(s.data);
x=[];
fs=250;
time=N/fs;
b=0:N-1; 
f=fs*b/N; %Ƶ������
t=26;
for i=1:t
    x(i,:)=s.data(i,53*2500:53*2500+N);
end

figure(10); %ԭʼ�ź�ͼ��
for i=1:t
    subplot(t/2,2,i);
    plot(x(i,:));
end

figure(2);
for i=1:t
    subplot(6,5,i);
  topoplot(abs(x(i,:)),'eeglab_chan26.locs', 'maplimits','maxmin');
  title(i);
end
% % 
% figure(2);%ԭʼ�źŸ���Ҷת�����ͼ��
% for i=1:t
%     F1=abs(fft(x(i,:)));
% end
% for i=1:t
%     subplot(t/2,2,i);
%     plot(f,F1);
% end
% figure(3);%ԭʼ�źŹ�����
% for i=1:t
%     subplot(t/2,2,i);
%     plot(f,F1.^2);
% end
% % %��ͨ�˲�
% Wn=[0.5 *2 30*2]/fs;%����ͨ��
% [k,l]=butter(1,Wn);
% x=filtfilt(k,l,x);
% f1_design=8; f2_design=30;
% window_k=kaiser(17,0.02);
% b_design=fir1(16,[f1_design,f2_design]/(fs/2),window_k);
% 
% x=double(x);
% x=filtfilt(b_design,1,x')';
% 
% figure(3);
% for i=1:t
%     subplot(15,2,i);
%     plot(x(i,:));
% end
% figure(4);%����Ҷת�����ͼ��
% for i=1:t
%     subplot(15,2,i);
%     plot(f,abs(fft(x(i,:))));
% end
% 
%%*******ICA�㷨��********
% figure(2);  
%fsatICA
 [S,A,W] = fastica(x);  

%  [A,W,S]=jade_wu(x,26);     %�˶�����1
%  [A,W,S]=jade_wu(x,16);   %��������2
%  [A,W,S]=jade_wu(x,17);   %�Ե��ź� 3
% [A,S]=Infomax_ICA_teacher(x,5000,0.03)%ԭʼ�źţ�����������ѧϰ����
 
figure(3);   %�����ͼ��
for i=1:t
    subplot(t/2,2,i);
    plot(S(i,:));
end
%
figure(4);
for i=1:t
    subplot(6,5,i);
  topoplot(abs(S(i,:)),'eeglab_chan26.locs', 'maplimits','maxmin');
  title(i);
end

% 
F=[];
for i=1:t
    F(i,:)=abs(fft(S(i,:)));
end
figure(5);
for i=1:t
    subplot(t/2,2,i);
    plot(f,F(i,:));
end

figure(6)
for i=1:t
    subplot(t/2,2,i);
    plot(f,F(i,:).^2);
end
% % 
% %****************��ͨ�˲���*******************
 filtercoefficient2=fir1(8,[0.1/(fs/2),30/(fs/2)],'bandpass');  %��ͨ�˲�
 for i=1:t
 S(i,:)=filtfilt(filtercoefficient2,1,S(i,:));%x1�˲�������
 end
figure(7); 
for i=1:t 
    subplot(t/2,2,i);
    plot(S(i,:));
end
figure(8); 
for i=1:t 
    subplot(t/2,2,i);
    plot(f,abs(fft(S(i,:)).^2));
end

%  filtercoefficient2=fir1(8,[0.1/(fs/2),30/(fs/2)],'bandpass');  %��ͨ�˲�
%  for i=1:t
%  x(i,:)=filtfilt(filtercoefficient2,1,x(i,:));%x1�˲�������
%  end
% figure(7); 
% for i=1:t 
%     subplot(t/2,2,i);
%     plot(x(i,:));
% end
% figure(8); 
% for i=1:t 
%     subplot(t/2,2,i);
%     plot(f,abs(fft(x(i,:)).^2));
% end

% 
% %%******ȥ���۵硪�������Ͷ�ֵ,գ��ͨ���Ͷ���󣬹ʽ��Ͷ�����ͨ������*******
% d=[];
% for i=1:t
%    d(i,:)=kurt(S(i,:));
% end
% max=d(1,:);
%   for i=1:t
%       if d(i,:)>max
%           max=d(i,:);
%       end
%   end

% figure(9); %ȥ���۵���ź�ͼ��
% for i=1:t 
%     subplot(t/2,2,i);
%     plot(S(i,:));
% end


% %%**********ȥ���ĵ�**************
% %����QRS����Ȳ�����0.12s(30��������)
% %PR��Χ0.12~0.20s��QT��Χ0.32~0.44s
% %�õ�����ͨ���ĵĲ���Ͳ��ȣ������ĵ��������
% %1.���岨�ȳɶԳ���
% %2.���岨�ȵķ�ֵ��һ������
% %3.���岨�ȵ���������������й�
% %***************�󲨷岨����*********** 
%    for a=1:t
%      ecg=S(a,:);
% %     ecg=S(22,:);
%     %IndMax��IndMinΪ�����жϵĲ���Ͳ��Ⱥ�����
%     IndMin=find(diff(sign(diff(ecg)))>0)+1;
%     IndMax=find(diff(sign(diff(ecg)))<0)+1;
% 
%     k1=0;  %�жϽϴ��ֵ�Ĳ���������ecg��Ϊ�ж������Ĵ�����
%     k2=0;
%     a1=[]; %�õ����崦��Ӧ�ķ�ֵ
%     a2=[]; %�õ����ȴ���Ӧ�ķ�ֵ
%     b1=[]; %�õ����ȴ���Ӧ�ĺ�����
%     b2=[];
%     for i=1:N
%        a1=ecg((IndMax)); 
%        a2=ecg(IndMin); 
%     end
%      m=length(IndMax); %mΪ�������Ͳ������еĽ�Сֵ
%     if length(IndMin)<m
%         m=length(IndMin);
%     else
%         m=length(IndMax);
%     end
%     for i=1:m
%         if a1(i)>2.5
%           k1=k1+1; 
%           b1(k1)=IndMax(i);
%         else if a2(i)<-0.45
%           k2=k2+1;
%           b2(k2)=IndMin(i);
%             end
%        end
%     end
%     %�õ�һ����ֵ��Χ�ڵĲ��岨���� ���⼸�����ݽ��д���
%     %������12������㣬�ж���12����ĺ������Ƿ������250���ڣ�������˵������ͬһ����������ȥ���ظ��ĵ�
%     k3=0;%������
%     s1=[];%�������
%     for i=1:k1-1
%         for j=i+1:k1
%             p=b1(j)-b1(i);
% %           if p>=252
%              if p>=150
%               k3=k3+1;
%               s1(k3)=b1(j);
%               i=j;
%           end 
%         end
%         break;
%     end
%     
%     k4=0;%������
%     s2=[];%�������
%     for i=1:k2-1
%         for j=i+1:k2
%             p=b2(j)-b2(i);
% %           if p>=252
%               if p>=140
%               k4=k4+1;
%               s2(k4)=b2(j);
%               i=j;
%           end 
%         end
%         break;
%     end
%     
%     if (abs(k3-k4)<=3)&&k3>=time&&k3<=(time*2)
%         disp('this is ecg');   %�жϳ��ĵ�ͨ��
%         S(a,:)=0;
%     else
%         disp('this is not ecg');
%     end 
%    end  
 
%   S(26,:)=0;
% figure(10);%ȥ�����ź���ź�
% for i=1:t 
%     subplot(15,2,i);
%     plot(S(i,:));
% end

 S(1,:)=0;
% S(6,:)=0;
% % S(4,:)=0;
% S(26,:)=0;
% S(19,:)=0;eegla

%******�Դ����Ĳ��ν�����ICA*******
figure(10);
Y=A*S;
for i=1:t  
    subplot(t/2,2,i);
    plot(Y(i,:));
end
figure(11);
for i=1:t
    subplot(6,5,i);
  topoplot(abs(Y(i,:)),'eeglab_chan26.locs', 'maplimits','maxmin');
  title(i);
end

% 
%%***********�Ե��ź����� ��17������������·�����źţ�***************
% s=load('EEG19.mat');
% N=length(s.EEG19);
%  fs=512;
% time=N/fs;
% b=0:N-1; 
% f=fs*b/N; %Ƶ������
% t=17;
% for i=1:19
%     x(i,:)=s.EEG19(i,1:2048);
% end
% x(4,:)=[];
% x(9,:)=[];