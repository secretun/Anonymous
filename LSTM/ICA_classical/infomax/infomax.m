%load d:\mat\laughter;
%s7=y(1:20000)';
%load d:\mat\handel;
%s8=y(1:20000)';
%s3=rand(1,20000);
%sr=rand(1,20000);
%rs=reshape(sr,400,50);
%s4=reshape(sort(rs),1,20000);
f1=fopen('d:\eeg\dat1','r');
date11=fread(f1,7500,'int16');
f2=fopen('d:\eeg\dat2','r');
date12=fread(f2,7500,'int16');
f3=fopen('d:\eeg\dat3','r');
date13=fread(f3,7500,'int16');
b=fir1(25,120/125);
for i=1:1500,
   date1(i)=date11(i+6000);
end
s1=wkeep(filter(b,1,date1),1500);

for i=1:1500,
   date2(i)=date12(i+6000);
end
s2=wkeep(filter(b,1,date2),1500);

for i=1:1500,
   date3(i)=date13(i+6000);
end
s3=wkeep(filter(b,1,date3),1500);

ss=[s1'/(sqrt(cov(s1))),s2'/(sqrt(cov(s2))),s3'/(sqrt(cov(s3)))]';

S1=ss(:,1:1500);
%S1(1,1)=17;
%S1(1,980)=6;
%S1(1,500)=10;
%S1(1,2800)=12;
%S1(1,4500)=10;
%S1(1,4000)=7;
%S1(2,2300)=10;
%S1(2,1700)=7;
%S1(2,500)=10;
%S1(2,2500)=8;
%S1(2,40)=11;
%S1(2,300)=10;
figure(1)
for i=1:3,
   subplot(3,1,i)
   plot(S1(i,:));
end
%S1(2,3000)=15;
a=[0.43,0.57,0.43
   0.78,0.65,-0.51
   -0.41,0.55,0.49];
wc=remmean(S1);
figure(2)
for i=1:3,
   subplot(3,1,i)
   plot(wc(i,:));
end
[N,P]=size(wc);
w=rand(3,3)-0.5;
u=zeros(N,P);
BI=40*eye(N);
NL=400;
for i=1:NL
   u=w*wc;
   %for j=1:N,
     %u1(j,:)=medfilt1(u(j,:));
      %end
        m2=mean(u'.^2).^2; 
        m4= mean(u'.^4);
        K= sign(diag((m4./m2)-3.0));  
        dw=(BI-K*tanh(u)*u'-u*u')*w;
        w=w+0.0005*dw;
        tw{i}=w;
        tk{i}=diag(K);
     end
     figure(3)
     for i=1:3,
   subplot(3,1,i)
   plot(u(i,:));
end
e=zeros(1,NL);
for i=1:NL,
   e1=0;
e2=0;

   p=tw{i}*a;
   ws=double(p);
   lm=max(abs(ws));
   hm=max(abs(ws'));
   for k=1:N,
      h=hm(k);
      l=lm(k);
      for j=1:N,
         e1=e1+abs(p(k,j))/h;
         e2=e2+abs(p(j,k))/l;
       end
       e1=e1-1;
       e2=e2-1;
   end
   e(i)=e1+e2;
   
end
figure(4)
plot((e),'r');
hold on
p=zeros(N,NL);
figure(5)
for i=1:NL,
   p(:,i)=tk{i};
end

   for j=1:N,
     subplot(N,1,j)
     plot(p(j,:))
   end

