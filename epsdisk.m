load epsdisk_data.mat;

figure;
pcolor(epsdiskR(2:end-1, 2:end-1)-epsdiskold(2:end-1, 2:end-1)); axis equal; 
axis off; shading interp; piyg; cb = colorbar; set(cb, 'location', 'east');
caxis( 0.025*[-1, 1]);
text(96.9, 100.5, '$\delta \varepsilon$', 'fontsize', 20, 'interpreter', 'latex');

figure;
pcolor(epsdiskI(2:end-1, 2:end-1)); axis equal; axis off; shading interp; 
piyg; cb = colorbar; set(cb, 'location', 'east'); 
text(96.9, 100.5, '$\delta \varepsilon$', 'fontsize', 20, 'interpreter', 'latex');
caxis( 1e-3*0.12*[-1, 1]);