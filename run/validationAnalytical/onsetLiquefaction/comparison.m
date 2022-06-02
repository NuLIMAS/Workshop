
f1 = figure;
firstax = axes (f1, 'FontSize', 16); 
c1=load('postProcessing/Profiles/120/Profile1_pE_liqueCriteria.xy' );
c2=load('postProcessing/Profiles/240/Profile1_pE_liqueCriteria.xy' );
c3=load('postProcessing/Profiles/480/Profile1_pE_liqueCriteria.xy' );
c4=load('postProcessing/Profiles/900/Profile1_pE_liqueCriteria.xy' );

plot(c1(:,3), c1(:,1),'b','lineWidth',2)
hold on
plot(c2(:,3), c2(:,1),'r','lineWidth',2)
plot(c3(:,3), c3(:,1),'g','lineWidth',2)
plot(c4(:,3), c4(:,1),'m','lineWidth',2)

d1=load('t2.csv');
d2=load('t4.csv');
d3=load('t8.csv');
d4=load('t15.csv');

plot(d1(:,1), d1(:,2),'b--','lineWidth',2)

plot(d2(:,1), d2(:,2),'r--','lineWidth',2)
plot(d3(:,1), d3(:,2),'g--','lineWidth',2)
plot(d4(:,1), d4(:,2),'m--','lineWidth',2)

plot([1,1],[0,-1],'k--')

%xlabel({'$\overline{p}/\sigma^{'}_0$'},'Interpreter','latex')
x1=xlabel({ strcat('P/ \sigma',"'" ,'_{0} ')},'Interpreter','tex')
y1=ylabel('z/h')
lgd=legend({'t=2min','t=4min','t=8min','t=15min'},'location','southeast' );
set(lgd,'fontsize',16)
 
set(x1, 'FontSize', 20)
set(y1, 'FontSize', 20)
set(x1,'FontWeight','bold')
set(y1,'FontWeight','bold')
 
 
 
p = 0:0.5:5;
k = 2*p + 3;
o = 0:0.5:5;
r = 2*o +3.2;
 
 hold on 
%secondax = copyobj(firstax, gcf);
%delete( get(secondax, 'Children'))
%H1 = plot(k, '-', 'LineWidth', 2, 'Color', [0 0 0]/255, 'Parent', secondax, 'Visible', 'off');
%H2 = plot(r, '--', 'LineWidth', 2, 'Color', [0 0 0]/255, 'Parent', secondax, 'Visible', 'off');
%set(secondax, 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off','fontsize',16) 


%legend (secondax, [H1 H2], 'Numerical', 'Analytical', 'location', 'southeast');
t=title({'Comparison of Numerical (-) and Analytical^1 (--) solution'} , 'Interpreter', 'tex')
set(t,'fontsize',16)
saveas(figure(1), 'pbysigma0_z2.png')
