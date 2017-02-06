Ns = [25, 30, 35, 40, 45, 50:10:100, 111:10:211, 231:10:291, 301];
% removed data point at 221, splitting too low
L = 3.2;


% to print out a row of seven digits, do
% grep "#0" tally | sed s/.*"w = "// | sed s/" +".*//
% then format long, then do "single( [ %paste into matlab ... ].')"

w0 = [3.6504741 3.7351287 3.793846  3.7834386 3.7881012 3.8127954   3.8527558   3.8312275   3.8484561   3.8510289   3.8591478   3.8644125   3.8641462   3.8632023   3.8680258   3.8659501   3.8751945   3.8729823   3.8709669   3.8750095     3.8717651   3.8725128   3.8717320   3.8747873   3.8733115   3.8761141   3.8745012   3.8747311   3.8743253   3.8752811];
w1 = [3.5476631 3.6654981 3.7447976  3.761753  3.8146588 3.8413866   3.8487830   3.8443699   3.8554006   3.8639562   3.8639083   3.8662591   3.8685060   3.8682823   3.8772192   3.8731735   3.8812177   3.8767133   3.8752673   3.8759880     3.8735201   3.8753359   3.8764963   3.8806953   3.8747623   3.8771136   3.8771636   3.8751304   3.8754592   3.8768506];
% doing real parts only for simplicity
% Npml numbers up to 301 given in circa notes090515
% 35: 3, 50: 4
% for second round of data:
% 311: 25 ...  370: 30 ... 450: 35 ... 540 40

NsB = [311 321 330 350 370 390 420 450 480 510 540];
% put last one at multiple of 100 for neater xlim
w0B = [3.8764421 3.8763336 3.8747996 3.8759104 3.8769372 3.8770653 ...
    3.874636 3.8770011 3.8763646 3.8780576 3.8779151];
w1B = [3.8795585 3.878284 3.8757862 3.8768909 3.877065 3.8773772 ...
    3.8765516 3.8775743 3.8770989 3.8782522 3.8780536];

Ns = [Ns, NsB];
w0 = [w0, w0B];
w1 = [w1, w1B];
hs = L ./ (Ns-1);

set(0, 'DefaultAxesFontSize',20);

splitting = abs(w0 - w1);

f=figure('Visible','off');
loglog(1./hs, splitting, 'b-', 1./hs, hs.^2, 'r--' );
% note 0 and 1 in passive could've been switched, so assume wcos > wsin always or vice versa always

%xlim([min(Ns.^2), max(Ns.^2)]);
%ylim([min(splitting), max(splitting)]);

xlabel('Resolution ($1/h$)', 'interpreter', 'latex');
ylabel('$\omega_2 - \omega_1$', 'interpreter', 'latex');

xlim([min(1./hs), max(1./hs)]);
% otherwise a lot of whitespace

annotation('textbox', [0.4, 0.2, 0.3, 0.2],'fontsize', 20, ...
    'linestyle', 'none', 'interpreter', 'latex', ...
    'string', '$O(h^{-2})$' );

%can't use pastecylinder here because y axis is in log units
ratio = get(gca, 'plotboxaspectratio');
w = 0.16;
h = w * ratio(1)/ratio(2);
xrel = 0.6; yrel = 0.63;
annotation('ellipse', [xrel yrel w h], 'facecolor', [1 1 1]*0.8);

% text uses axes units, annotation uses relative figure units.
% also make sure width (3rd argument) is large enough or it will throw a
% very deceptive "interpreter invalid syntax" error
annotation('textbox', [xrel-0.05, yrel-0.05, 0.3, 0.2],'fontsize', 20, ...
    'linestyle', 'none', 'interpreter', 'latex', ...
    'string', '$\varepsilon=5+D_0\gamma H$' );

annotation('textbox', [xrel+0.2, yrel-0.1, 0.3, 0.2],'fontsize', 20, ...
    'linestyle', 'none', 'string', 'air' );


print -painters -dpdf cylinder_resolution.pdf
