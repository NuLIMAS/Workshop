f1=load('postProcessing/Profiles/5/Profile1_p_tauXZ.xy');

d=f1(:,1)/50;
p1=f1(:,2)/18544;
plot(p1,d,'r','linewidth',1.5)
figure(1)
analytical = load('myFile.txt');
hold on
plot(analytical.temp1(:,1),analytical.temp1(:,2),'--b','linewidth',1.5)

legend('BiotFoam','Analytical')

xlabel('p/p0')
ylabel('z/h')
saveas(figure(1),'pressure.png')

f2=load('postProcessing/Profiles/5/Profile2_p_tauXZ.xy');

figure(2)
t1=-f2(:,3)/18544;
plot(t1,d,'r','linewidth',1.5)


hold on
plot(analytical.temp2(:,1),analytical.temp2(:,2),'--b','linewidth',1.5)
legend('BiotFoam','Analytical')
xlabel('tauXZ/p0')
ylabel('z/h')
saveas(figure(2),'stress.png')
