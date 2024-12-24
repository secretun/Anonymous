%load d:\mat\laughter;
%s7=y(1:20000)';
f1=fopen('d:\wu\dat1','r');
s6=fread(f1,6000,'int16');
f2=fopen('d:\wu\dat2','r');
s7=fread(f2,6000,'int16');
s6=s6';
s7=s7';
load d:\mat\handel;
s8=y(1:6000)';
s3=rand(1,6000);
sr=rand(1,6000);
rs=reshape(sr,300,20);
s4=reshape(sort(rs),1,6000);
s2=(sin(2*pi*0.05*(1:6000)/10));
%s3=(cos(2*pi*0.1*(1:6000)/10));

ss=[s2'/(sqrt(cov(s2))),s3'/(sqrt(cov(s3)))]';

S1=ss(:,1:6000);
[N,P]=size(S1);
wc=zeros(N,P);
%S1(1,200)=7;
%S1(1,980)=6;
%S1(1,500)=4;
%S1(1,2800)=6;
%S1(1,3500)=8;
%S1(1,3000)=7;
%S1(2,2300)=6;
%S1(2,1700)=7;
%S1(2,700)=5;
%S1(2,2500)=6;
%S1(2,40)=7;
%S1(2,300)=5;
figure(1)
for i=1:N,
   subplot(N,1,i)
   plot(S1(i,:));
end
a1=[0.13,0.37
   0.15,0.25];
a2=[0.28,0.57
  0.73,0.35];

wc(:,1:3000)=remmean(a2*S1(:,1:3000));
wc(:,3001:6000)=remmean(a1*S1(:,3001:6000));

figure(2)
for i=1:N,
   subplot(N,1,i)
   plot(wc(i,:));
end
w=rand(N);
u=zeros(N,P);
BI=eye(N);
NL=P;
ku=0;
m2=0;
m4=0;
dwo=0;
for i=1:P
   u=w*wc(:,i);
   %for j=1:N,
     %u1(j,:)=medfilt1(u(j,:));
     %end
        u2=u.^2;   
        m2=(u2+m2*(i-1))/i; 
        u4=u.^4;
        m4=(u4+m4*(i-1))/i;
        ku=(m4./(m2.^2))-3.0;
        K= sign(diag(ku));  
        dw=(BI-K*tanh(u)*u'-u*u')*w;
        w=w+0.005*dw;
        dwo=dw;
        tw{i}=w;
        tk{i}=diag(K);
        us(:,i)=w*wc(:,i);
     end
     U=w*wc;
     figure(3)
     for i=1:N,
   subplot(N,1,i)
   plot(U(i,:));
end
e=zeros(1,NL);
figure(5)
for i=1:NL,
   p(:,i)=tk{i};
end

   for j=1:N,
     subplot(N,1,j)
     plot(p(j,:))
     axis([1,P,-2,2]);
   end

figure(6)
   for j=1:N,
     subplot(N,1,j)
     plot(us(j,:))
   end
   figure(7)
   for j=1:P,
      wp=double(tw{j});
      wpp=reshape(wp,1,N*N);
      wp1(j)=wpp(1);
      wp2(j)=wpp(2);
      wp3(j)=wpp(3);
      wp4(j)=wpp(4);
   end
   subplot(4,1,1)
   plot(wp1);
   subplot(4,1,2)
   plot(wp2);
   subplot(4,1,3)
   plot(wp3);
   subplot(4,1,4)
   plot(wp4);
   



      