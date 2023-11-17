% Load pressure and stress data from BiotQSFoam simulation
f1 = load('postProcessing/Profiles/5/Profile1_p_tauXZ.xy');
f2 = load('postProcessing/Profiles/5/Profile2_p_tauXZ.xy');

% Extract data for pressure profile
d = f1(:, 1) / 50;   % Depth normalized by h
p1 = f1(:, 2) / 18544;  % Pressure normalized by p0

% Extract data for stress profile
t1 = -f2(:, 3) / 18544;  % Negative of stress normalized by p0

% Load analytical solution
analytical = load('analyticalSolution.txt');

% Plot pressure profile
figure(1)
plot(p1, d, 'r', 'linewidth', 1.5)
hold on
plot(analytical.temp1(:, 1), analytical.temp1(:, 2), '--b', 'linewidth', 1.5)
legend('BiotQSFoam', 'Analytical', 'fontsize', 18)
xlabel('p/p0', 'fontsize', 18)
ylabel('z/h', 'fontsize', 18)
xlim([0 1])
saveas(figure(1), 'pressure.png')

% Plot stress profile
figure(2)
plot(t1, d, 'r', 'linewidth', 1.5)
hold on
plot(analytical.temp2(:, 1), analytical.temp2(:, 2), '--b', 'linewidth', 1.5)
legend('BiotQSFoam', 'Analytical', 'fontsize', 18)
xlabel('tauXZ/p0', 'fontsize', 18)
ylabel('z/h', 'fontsize', 18)
xlim([min(t1) max(t1)])  % Adjusted to the range of t1
set(gca, 'fontsize', 18)
saveas(figure(2), 'stress.png')

