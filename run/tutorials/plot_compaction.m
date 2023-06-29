%close all
pp=false;

dir= 'case3_compaction'
R1=load(strcat(dir,'/postProcessing/Probes/0.4/pMean'));
R=load(strcat(dir,'/postProcessing/Probes/0/p'));

%R=[R;R11];
fac=2;
st=8;
%pkg load io
c='rgbk';
%c='cmyk';

figure(1)
hold on
plot(R1(:,1)/60,R1(:,2)/1000,'r','lineWidth',1.25)
plot(R1(:,1)/60,R1(:,3)/1000,'g','lineWidth',1.25)
plot(R1(:,1)/60,R1(:,4)/1000,'b','lineWidth',1.25)
plot(R1(:,1)/60,R1(:,5)/1000,'k','lineWidth',1.25)

xlabel('Time(m)','fontSize', 22)
ylabel('P (kPa)','fontSize', 22)
h= legend('H=0.1m','H=0.15m','H=0.20m','H=0.25m')
set(h,'fontSize', 22)

figure(2)
for i= 1:size(R,2)-1

sample =R(:,i+1)/1000;
step = st;

sample = sample(1:floor(end/step)*step);

avg = mean(reshape(sample,step,[]));
% Let's take a look!
hold on

plot((step/2:step:length(sample))*0.1*fac,avg,strcat(c(i),'-'),  'linewidth',4)
%text(100,max(Exp1(:,i*2+1)),strcat('R2 = ',num2str(R2)));
xlim([0 30])
%ylim([0 15])
xlabel('Time(s)','fontSize', 22)
ylabel('P (kPa)','fontSize', 22)
h= legend('Z=0.1m','Z=0.15m','Z=0.20m','Z=0.25m');
set(h,'fontSize', 22);
set(gca,'fontSize', 18)
 end
%print -djpg figure2


depth = load(strcat(dir,'/liquefactionDepth.txt'));
figure(3)

plot(depth(:,1),depth(:,2),'linewidth',1.6)
hold on
%plot([0 6],[-0.06 -0.06], 'k--')
xlabel('Time(s)','fontSize', 22)
ylabel('Z (m)','fontSize', 22)
set(gca,'fontSize', 18)
xlim([0 30])
%ylim([-0.3 0])
%print -djpg figure1



