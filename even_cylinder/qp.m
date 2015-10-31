D0s = [0.0006708591461, 0.0006342428932, 0.000634240667265435];
D1s = [0.0005983088105, 0.0006342384441, 0.000634240667310395];

iterations = [0, 1, 2];

w0s = [3.8680985, 3.8658991828946876, 3.86589919162261];
w1s = [3.8637302 3.8658992003312949 3.8658991916221];

% keep data for deps in case need to plot it
% issue: only two data points while 3 data points for splitting
denom = 288724;
% integral of eps^2
deps01 = 0.393650403992694;
% integral of first iteration's deps
deps12 = 6.247877277207378e-12;
% integral of change in deps from first to second iteration


% notation for dD / D, i.e. fractional change. for D in denom use average
dlnD = abs( D1s - D0s ) .*2 ./ ( D1s + D0s );
dlnw = abs( w1s - w0s ) .*2 ./ ( w1s + w0s );

lw = 2.0; % linewidth
set(0, 'DefaultAxesFontSize',20);

figure;
% hold on does not seem to work with semilogy
semilogy(iterations, dlnD, 'bo-', iterations, dlnw, 'bo-',...
    'linewidth', lw, 'markersize', 10);
set(gca, 'xtick', [0 1 2]);
xlabel('QP iterations');
ylabel('relative splitting');
box on;

text(1.3, 8.9e-12, '$\omega_t$', 'fontsize', 24, 'interpreter', 'latex');
text(1.4, 6.4e-7, '$D_t$', 'fontsize', 24, 'interpreter', 'latex');

print -painters -dpdf -r200 ~/Desktop/thesis/figures/qp.pdf
