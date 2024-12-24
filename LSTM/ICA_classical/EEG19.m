clear;
close all;

%***********运动想象数据2（不含脉搏，共16导）***************
s=loadcnt('D1.cnt','data','int16');%导入cnt文件
N=2500;
% N=length(s.data);
x=[];
fs=250;
time=N/fs;
b=0:N-1; 
f=fs*b/N; %频率序列
t=26;
for i=1:t
    x(i,:)=s.data(i,53*2500:53*2500+N);
end

figure(10); %原始信号图像
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
% figure(2);%原始信号傅里叶转换后的图像
% for i=1:t
%     F1=abs(fft(x(i,:)));
% end
% for i=1:t
%     subplot(t/2,2,i);
%     plot(f,F1);
% end
% figure(3);%原始信号功率谱
% for i=1:t
%     subplot(t/2,2,i);
%     plot(f,F1.^2);
% end
% % %带通滤波
% Wn=[0.5 *2 30*2]/fs;%设置通带
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
% figure(4);%傅里叶转换后的图像
% for i=1:t
%     subplot(15,2,i);
%     plot(f,abs(fft(x(i,:))));
% end
% 
%%*******ICA算法们********
% figure(2);  
%fsatICA
 [S,A,W] = fastica(x);  

%  [A,W,S]=jade_wu(x,26);     %运动想象1
%  [A,W,S]=jade_wu(x,16);   %运行想象2
%  [A,W,S]=jade_wu(x,17);   %脑电信号 3
% [A,S]=Infomax_ICA_teacher(x,5000,0.03)%原始信号，迭代次数，学习步长
 
figure(3);   %分离后图像
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
% %****************带通滤波器*******************
 filtercoefficient2=fir1(8,[0.1/(fs/2),30/(fs/2)],'bandpass');  %带通滤波
 for i=1:t
 S(i,:)=filtfilt(filtercoefficient2,1,S(i,:));%x1滤波后数据
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

%  filtercoefficient2=fir1(8,[0.1/(fs/2),30/(fs/2)],'bandpass');  %带通滤波
%  for i=1:t
%  x(i,:)=filtfilt(filtercoefficient2,1,x(i,:));%x1滤波后数据
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
% %%******去除眼电——计算峭度值,眨眼通道峭度最大，故将峭度最大的通道置零*******
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

% figure(9); %去除眼电干扰后图像
% for i=1:t 
%     subplot(t/2,2,i);
%     plot(S(i,:));
% end


% %%**********去除心电**************
% %正常QRS波宽度不超过0.12s(30个采样点)
% %PR范围0.12~0.20s，QT范围0.32~0.44s
% %得到所有通道的的波峰和波谷，满足心电的条件：
% %1.波峰波谷成对出现
% %2.波峰波谷的幅值有一定限制
% %3.波峰波谷的数量与采样点数有关
% %***************求波峰波谷数*********** 
%    for a=1:t
%      ecg=S(a,:);
% %     ecg=S(22,:);
%     %IndMax和IndMin为初步判断的波峰和波谷横坐标
%     IndMin=find(diff(sign(diff(ecg)))>0)+1;
%     IndMax=find(diff(sign(diff(ecg)))<0)+1;
% 
%     k1=0;  %判断较大幅值的波数（对于ecg即为判断心跳的次数）
%     k2=0;
%     a1=[]; %得到波峰处对应的幅值
%     a2=[]; %得到波谷处对应的幅值
%     b1=[]; %得到波谷处对应的横坐标
%     b2=[];
%     for i=1:N
%        a1=ecg((IndMax)); 
%        a2=ecg(IndMin); 
%     end
%      m=length(IndMax); %m为波峰数和波谷数中的较小值
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
%     %得到一定阈值范围内的波峰波谷数 对这几个数据进行处理
%     %例如有12个波峰点，判断这12个点的横坐标是否相差在250以内，若在则说明属于同一个波，则需去除重复的点
%     k3=0;%计数器
%     s1=[];%起点坐标
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
%     k4=0;%计数器
%     s2=[];%起点坐标
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
%         disp('this is ecg');   %判断出心电通道
%         S(a,:)=0;
%     else
%         disp('this is not ecg');
%     end 
%    end  
 
%   S(26,:)=0;
% figure(10);%去除干扰后的信号
% for i=1:t 
%     subplot(15,2,i);
%     plot(S(i,:));
% end

 S(1,:)=0;
% S(6,:)=0;
% % S(4,:)=0;
% S(26,:)=0;
% S(19,:)=0;eegla

%******对处理后的波形进行逆ICA*******
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
%%***********脑电信号数据 共17导（不包括两路无用信号）***************
% s=load('EEG19.mat');
% N=length(s.EEG19);
%  fs=512;
% time=N/fs;
% b=0:N-1; 
% f=fs*b/N; %频率序列
% t=17;
% for i=1:19
%     x(i,:)=s.EEG19(i,1:2048);
% end
% x(4,:)=[];
% x(9,:)=[];