theta = [1.62 1.61 1.6 1.59 1.58 1.57 1.56 1.55 1.54 1.53];
w0 = [3.83234175947220 3.83234114838034 3.83234068519080 3.83234036219887 3.83234017225512 3.83234010871266 3.83234016538772 3.83234033650689 3.83234061668336 3.83234100087863];
w1 = [3.83627464497351 3.83627525786377 3.83627572237908 3.83627604624852 3.83627623664426 3.83627630023395 3.83627624322236 3.83627607140030 3.83627579017233 3.83627540459329];
set(0, 'DefaultAxesFontSize',20);

plot(theta, w0-w0(6), 'b-', theta, w1-w1(6), 'b-', 'linewidth', 2);

xlim([min(theta), max(theta)]);
xlabel('$\theta$', 'interpreter', 'latex');
ylabel('$\omega_\mu$ (shifted)', 'interpreter', 'latex');

text(1.576, 8.1e-07, '$\omega_2 - \omega_2|_{\theta=1.57}$', ...
    'fontsize', 20, 'interpreter', 'latex');

text(1.576, -7.6e-07, '$\omega_1 - \omega_1|_{\theta=1.57}$', ...
    'fontsize', 20, 'interpreter', 'latex');

print -painters -dpdf ~/Desktop/thesis/figures/even_cylinder_stationary.pdf