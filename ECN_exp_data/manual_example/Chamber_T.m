figure(1)
h=slice(X,Y,Z,T,[],[],30);
set(h, 'EdgeColor', 'none');
colorbar
caxis([200 350])
zlabel('Z axis [mm]')
hold on
h2=slice(X,Y,Z,T,0,[],[]);
set(h2, 'EdgeColor', 'none');
colorbar
caxis([200 350])
ylabel('Y axis [mm]')
hold on
h3=slice(X,Y,Z,T,[],0,[]);
set(h3, 'EdgeColor', 'none');
colorbar
caxis([200 350])
xlabel('X axis [mm]')
set(gcf,'color','w')