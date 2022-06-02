A=load('~/Projects/NuLIMAS_GIT/run/postLiquefaction/testCaseBE/postProcessing/Probes/0/pE');
B=load('~/Projects/NuLIMAS_GIT/run/postLiquefaction/testCasePL/postProcessing/Probes/0/pE');


for i= 1:6
plot(A(:,1),A(:,i),'r')
hold on
plot(B(:,1),B(:,i),'b')
end
hold off
legend('Biot', 'Navier-Stokes ')
