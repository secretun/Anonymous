load d:\mat\laughter;
s7=y(1:20000)';
load d:\mat\handel;
s8=y(1:20000)';
s3=rand(1,20000);
sr=rand(1,20000);
rs=reshape(sr,400,50);
s4=reshape(sort(rs),1,20000);

ss=[s3'/(sqrt(cov(s3))),s4'/(sqrt(cov(s4))),s7'/(sqrt(cov(s7))),s8'/(sqrt(cov(s8)))]';
figure(1)
for i=1:4,
   subplot(4,1,i)
   plot(ss(i,1:4000));
end
S1=ss(:,1:4000);
S1(1,2000)=30;
S1(2,3000)=30;
a=rand(4,4);
wc=remmean(a*S1);
figure(2)
for i=1:4,
   subplot(4,1,i)
   plot(wc(i,1:4000));
end
[u,W]=ext_icao(wc);
figure(3)
for i=1:4,
   subplot(4,1,i)
   plot(u(i,:));
end
xlabel('Informax');
