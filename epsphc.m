load epsphc_data.mat;

% for some reason pml creeps into nonpml region
figure;
pcolor(epsphcR(2:end-1, 2:end-1)-epsphcold(2:end-1, 2:end-1)); axis equal; 
axis off; shading interp; piyg; cb = colorbar; set(cb, 'location', 'east');
caxis( 0.7*[-1, 1]);
text(90.3, 86.5, '$\delta \varepsilon$', 'fontsize', 20, 'interpreter', 'latex');

figure;
pcolor(epsphcI(2:end-1, 2:end-1)); axis equal; axis off; shading interp; 
piyg; cb = colorbar; set(cb, 'location', 'east'); caxis( 0.05*[-1, 1]);
text(90.3, 86.5, '$\delta \varepsilon$', 'fontsize', 20, 'interpreter', 'latex');