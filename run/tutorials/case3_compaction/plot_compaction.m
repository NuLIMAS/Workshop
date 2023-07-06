%Compaction

R=load('postProcessing/Probes/0/p');%% Read samples of pore pressure

fac=0.2;
step=8;
c='rgbk';


figure(1) %% Pore pressure vs time
for i= 1:size(R,2)-1
sample =R(:,i+1)/1000;
sample = sample(1:floor(end/step)*step);
avg = mean(reshape(sample,step,[]));
hold on
xt= (step/2:step:length(sample))*0.2;
plot(R(:,1),R(:,i+1)/1000,strcat(c(i),'-'),'HandleVisibility','Off')
plot(xt,avg,strcat(c(i),'-'),  'linewidth',4)
xlim([0 30])
ylim([-0.5 2])
xlabel('Time(s)','fontSize', 22)
ylabel('P (kPa)','fontSize', 22)
h= legend('Z=0.1m','Z=0.15m','Z=0.20m','Z=0.25m');
set(h,'fontSize', 22);
set(gca,'fontSize', 18)
end


depth = load('liquefactionDepth.txt');

figure(2) %% Liquefaction depth

plot(depth(:,1),depth(:,2),'linewidth',1.6)
hold on
xlabel('Time(s)','fontSize', 22)
ylabel('Z (m)','fontSize', 22)
set(gca,'fontSize', 18)
xlim([0 30])





