wd=double(wdd{3});
for i=1:1999,
   a1=norm(wd(:,i+1)-wd(:,i));
   a2=norm(wd(:,i+1)+wd(:,i));
   am(i)=min(a1,a2);
end
hold on
figure(4)
plot(am);

for i=1:2000,
   y1=wd(:,i)'*nv;
   k(i)=mean((y1.^4)')-3*(mean((y1.^2)'))^2;
end

figure(5)
hold on
plot(k)