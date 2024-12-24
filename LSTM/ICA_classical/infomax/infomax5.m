load d:\mat\laughter;
s7=y(1:20000)';
load d:\mat\handel;
s8=y(1:20000)';
s3=rand(1,20000);
sr=rand(1,20000);
rs=reshape(sr,400,50);
s4=reshape(sort(rs),1,20000);
s2=(sin(2*pi*0.05*(1:20000)/10));
ss=[s2'/(sqrt(cov(s2))),s3'/(sqrt(cov(s3)))]';


S1=ss(:,1:4000);
[N,P]=size(S1);
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
for i=1:N,
   subplot(N,1,i)
   plot(S1(i,:));
end
wc=zeros(N,P);
%S1(2,3000)=15;
a1=[0.13,0.37
   0.15,0.25];
wc(:,1:2000)=remmean(a1*S1(:,1:2000));
a2=[0.28,0.57
  0.73,0.35];
wc(:,2001:4000)=remmean(a2*S1(:,2001:4000));
figure(2)
for i=1:N,
   subplot(N,1,i)
   plot(wc(i,:));
end
[N,P]=size(wc);
w=eye(N);
u=zeros(N,P);
BI=40*eye(N);
NL=1000;
dwo=0;
for i=1:NL
   u=w*wc;
   %for j=1:N,
     %u1(j,:)=medfilt1(u(j,:));
      %end
        m2=mean(u'.^2).^2; 
        m4= mean(u'.^4);
        K= sign(diag((m4./m2)-3.0));  
        dw=(BI-K*tanh(u)*u'-u*u')*w;
        w=w+0.0005*dw+0.2*dwo;
        dwo=dw;
        tw{i}=w;
        tk{i}=diag(K);
     end
     figure(3)
     for i=1:N,
   subplot(N,1,i)
   plot(u(i,:));
end
e=zeros(1,NL);
for i=1:NL,
   e1=0;
e2=0;
   if(i<=2000)   
      p=tw{i}*a1;
   else     
      p=tw{i}*a2;
   end
   
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
figure(5)
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

