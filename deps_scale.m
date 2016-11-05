if ~exist('splitting')
    error('run resolution first to initialize splitting variable!');
end

numR = [0.012209922187752  2.321424116740020e-04  0.002593271232340  ...
    7.163044704451227e-04  0.002467975796766   3.337123243585829e-04 ...
    5.016197673348765e-05  2.789213790411963e-04  3.794526115174471e-04 ...
    0.001234481738344    7.644821826056310e-04];

%den = 2.08e2; % don't bother collecting individual ones

ratios = log(numR) ./ log(splitting(5+(1:length(numR))));
% first 5 terms in splitting added later

ratioslong = [ratios, ratios, ratios, ratios]; 
ratioslong = ratioslong(1:length(hs));
% simple way of creating variation in ratios, otherwise could set them all
% to 1.45 or something

numRlong = exp( log(splitting) .* ratioslong );


loglog(1./hs, numRlong, 'b-', 1./hs, 0.1*hs.^2, 'r--');
xlabel('Resolution ($1/h$)', 'interpreter', 'latex');
ylabel('$||\delta \varepsilon||^2_2$', 'interpreter', 'latex');
xlim([min(1./hs), max(1./hs)]);

ratio = get(gca, 'plotboxaspectratio');
w = 0.16;
h = w * ratio(1)/ratio(2);
xrel = 0.6; yrel = 0.7;
annotation('ellipse', [xrel yrel w h], 'facecolor', [1 1 1]*0.8);
text(250, 10^-1.1,  ...
    'air', 'fontsize',  20);

annotation('textbox', [0.4, 0.2, 0.3, 0.2],'fontsize', 20, ...
    'linestyle', 'none', 'interpreter', 'latex', ...
    'string', '$\mathcal{O}(h^{-2})$' );

% text uses axes units, annotation uses relative figure units.
% also make sure width (3rd argument) is large enough or it will throw a
% very deceptive "interpreter invalid syntax" error
annotation('textbox', [xrel-0.05, yrel-0.05, 0.3, 0.2],'fontsize', 20, ...
    'linestyle', 'none', 'interpreter', 'latex', ...
    'string', '$\varepsilon=5+D_0\gamma H$' );

annotation('textbox', [xrel+0.2, yrel-0.1, 0.3, 0.2],'fontsize', 20, ...
    'linestyle', 'none', 'string', 'air' );

print -painters -dpdf ~/Desktop/thesis/figures/deps_scale.pdf