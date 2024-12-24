wc2=demosig;
A=rand(4,4)
wc1=A*wc2;
wc=wc1(1:3,:);

figure(1)
for i=1:3,
subplot(3,1,i)
plot(wc(i,1:500));
end
xlabel('observed')
[w,s,u]=runica(wc,'extended',-3);
figure(2)
for j=1:3,
   subplot(3,1,j)
   plot(u(j,1:500));
end
figure(3)
for j=1:4,
   subplot(4,1,j)
   plot(wc1(j,1:500));
end
