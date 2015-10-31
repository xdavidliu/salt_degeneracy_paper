ds = [0 20 40 60 80 100];

ws = [3.8658991916226051  3.8658993348787880  3.8658994624233678 ...
    3.8658995830920593 3.8658996979613574 3.8658998076518456]; 
% lasing passive pole
w0s = [3.8658991916226051  3.8658993342776946  3.8658994601439440 ...
    3.8658995782505094 3.8658996897952793 3.8658997955197303]; 
% other passive pole

imw0s = [0  5.6694391e-10  1.878e-09 3.915e-09 6.49e-09 9.55e-09]; 
% imaginary part

% at d = 20, c = 578. c^2 = a * d. a = c^2 / d
% a = 578^2 / 20; d = 40; c = sqrt(a * d);
set(0, 'DefaultAxesFontSize',20);

plot(ds, ws - w0s, 'bo-', ds, imw0s, 'bo-', 'linewidth', 2);
ylim([0, 1.5e-8]);
xlabel('$d$', 'interpreter', 'latex');
ylabel('$\omega$ splitting', 'interpreter', 'latex');

annotation('textbox', [0.61  0.54 0.3 0.2],'fontsize', 24, ...
    'linestyle', 'none', 'interpreter', 'latex', ...
    'string', '$\delta \omega - \mathrm{Re} \delta \omega^{\prime}$' );

annotation('textbox', [0.755 0.274 0.3 0.2],'fontsize', 24, ...
    'linestyle', 'none', 'interpreter', 'latex', ...
    'string', '$\mathrm{Im} \delta \omega^{\prime}$' );


print -painters -dpdf ~/Desktop/thesis/figures/even_cylinder_stability.pdf