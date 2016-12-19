
np=round(Nx/(64/2));
cl = get(gca,'colororder');

figNo=1;

nRow=2;nCol=2;
figure(figNo);clf;
% pos1 = [440   378   920   380];
% set(gcf,'Position',pos1);
subplot(nRow,nCol,3)
pcolor2(X,Y,zeta/etarms);shading interp;
axis square;colorbar;
axis([0 x(end) 0 y(end)]);
title('$\zeta(x,y,t)/\eta_{\rm rms}$','fontsize',18,'interpreter','latex');
xlabel('$x/L$','fontsize',18,'interpreter','latex');
ylabel('$y/L$','fontsize',18,'interpreter','latex');
set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10])
clim=get(gca,'CLim');caxis([-1 1]*.805*min(abs(clim)))

subplot(nRow,nCol,4)
pcolor2(X,Y,zetaMean/etarms);shading interp;
axis square;colorbar;
axis([0 x(end) 0 y(end)]);
title('$\bar{\zeta}(x,y,t)/\eta_{\rm rms}$','fontsize',18,'interpreter','latex');
xlabel('$x/L$','fontsize',18,'interpreter','latex');
ylabel('$y/L$','fontsize',18,'interpreter','latex');
set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10])
clim=get(gca,'CLim');caxis([-1 1]*.805*min(abs(clim)))

subplot(nRow,nCol,2)
pcolor2(X,Y,psi/(ell^2*etarms));shading interp;
hold on;
quiver(X(1:np:end,1:np:end),Y(1:np:end,1:np:end),u(1:np:end,1:np:end),v(1:np:end,1:np:end),2,'k','linewidth',0.5)
hold off;
axis square;colorbar;
axis([0 x(end) 0 y(end)]);
title('$\psi(x,y,t)/(\ell_\eta^2 \eta_{\rm rms})$','fontsize',18,'interpreter','latex');
xlabel('$x/L$','fontsize',18,'interpreter','latex');
ylabel('$y/L$','fontsize',18,'interpreter','latex');
set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10])
hmax=max(abs(eta(:)));
levh=linspace(-hmax,hmax,21);
clim=get(gca,'CLim');caxis([-1 1]*.805*min(abs(clim)))

subplot(nRow,nCol,1)
pcolor2(X,Y,eta);shading interp;
axis square;colorbar;
axis([0 x(end) 0 y(end)]);
title('$\eta(x,y)/\eta_{\rm rms}$','fontsize',18,'interpreter','latex');
xlabel('$x/L$','fontsize',18,'interpreter','latex');
ylabel('$y/L$','fontsize',18,'interpreter','latex');
set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10])
clim=get(gca,'CLim');caxis([-1 1]*.805*min(abs(clim)))

figNo=2;

nRow=2;nCol=2;
figure(figNo);
subplot(nRow,nCol,2)
plot(mu*T(1:it),Qpsi(1:it),'color',cl(1,:),'linewidth',2);
hold on
plot(mu*T(1:it),beta*Ut(1:it),'color',cl(2,:),'linewidth',2);
hold off;
xlim(mu*[0 T(it)])
title('b: $Q_{\psi}$~,~r: $Q_U$','fontsize',18,'interpreter','latex');
xlabel('$\mu t$','fontsize',18,'interpreter','latex');
ylabel('$Q$','fontsize',18,'interpreter','latex');

subplot(nRow,nCol,4)
plot(mu*T(1:it),formstress(1:it)/F,'color',cl(1,:),'linewidth',2);
xlim(mu*[0 T(it)])
hold on;plot(mu*[T(1) tfin],[0 0],'--k');hold off;
xlabel('$\mu t$','fontsize',18,'interpreter','latex');
ylabel('$\langle\psi\eta_x\rangle/F$','fontsize',18,'interpreter','latex');
hold off;

subplot(nRow,nCol,3)
plot(mu*T(1:it),mu*Ut(1:it)/F,'color',cl(2,:),'linewidth',2);
hold on;plot(mu*[T(1) tfin],[0 0],'--k');hold off;
xlim([0 mu*T(it)]);
xlabel('$\mu t$','fontsize',18,'interpreter','latex');
ylabel('$\mu U/F$','fontsize',18,'interpreter','latex');
hold off;

subplot(nRow,nCol,1)
plot(mu*T(1:it),Epsi(1:it),'color',cl(1,:),'linewidth',2);
hold on
plot(mu*T(1:it),EU(1:it),'color',cl(2,:),'linewidth',2);
hold off;
title('b: $E_{\psi}$~,~r: $E_U$','fontsize',18,'interpreter','latex');
xlabel('$\mu t$','fontsize',18,'interpreter','latex');
ylabel('$E$','fontsize',18,'interpreter','latex');
xlim(mu*[0 T(it)]);
drawnow;
