ic8_l_2e4=load('V:\PLVF2\iso-octane\G2\primary\liquid_data_2e-4_thr.mat');
ic8_l_2e3=load('V:\PLVF2\iso-octane\G2\primary\liquid_data_2e-3_thr.mat');
ic8_v=load('R:\Spray G\19_05_Schlieren\iso-octane\normal\G2\vapor_data.mat');

tc_l_2e4=load('V:\PLVF2\3-components\G2\primary\liquid_data_2e-4_thr.mat');
tc_l_2e3=load('V:\PLVF2\3-components\G2\primary\liquid_data_2e-3_thr.mat');
tc_v=load('R:\Spray G\19_05_Schlieren\3-components\normal\G2\vapor_data.mat');

figure(1)
plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,16),'--','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,16),':','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,12),'color',[0 0 1.0000],'linewidth',3)
hold on
plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,16),'--','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,16),':','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,12),'color',[0 0.5 0],'linewidth',3)
hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,16)+ic8_l_2e4.liquid_data_low_thr(:,18),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,16)-ic8_l_2e4.liquid_data_low_thr(:,18),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,16)+ic8_l_2e3.liquid_data_high_thr(:,18),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,16)-ic8_l_2e3.liquid_data_high_thr(:,18),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,16)+tc_l_2e4.liquid_data_low_thr(:,18),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,16)-tc_l_2e4.liquid_data_low_thr(:,18),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,16)+tc_l_2e3.liquid_data_high_thr(:,18),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,16)-tc_l_2e3.liquid_data_high_thr(:,18),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,12)+ic8_v.vapor_data(:,14),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,12)-ic8_v.vapor_data(:,14),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,12)+tc_v.vapor_data(:,14),':','color',[0 0.5 0],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,12)-tc_v.vapor_data(:,14),':','color',[0 0.5 0],'linewidth',2)
hold on

legend('ic8(l)0.2e-3','ic8(l)2e-3','ic8(v)','3c(l)0.2e-3','3c(l)2e-3','3c(v)','Location','northeastoutside')
xlabel('Time aSOI [ms]')
ylabel('Length [mm]')
set(gcf,'color','w');
grid on
set(gca,'FontSize',18)
set(gca,'linewidth',1.5)
box on

figure(2)
plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,2),'--','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,2),':','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,2),'color',[0 0 1.0000],'linewidth',3)
hold on
plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,2),'--','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,2),':','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,2),'color',[0 0.5 0],'linewidth',3)
hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,2)+ic8_l_2e4.liquid_data_low_thr(:,5),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,2)-ic8_l_2e4.liquid_data_low_thr(:,5),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,2)+ic8_l_2e3.liquid_data_high_thr(:,5),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,2)-ic8_l_2e3.liquid_data_high_thr(:,5),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,2)+tc_l_2e4.liquid_data_low_thr(:,5),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,2)-tc_l_2e4.liquid_data_low_thr(:,5),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,2)+tc_l_2e3.liquid_data_high_thr(:,5),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,2)-tc_l_2e3.liquid_data_high_thr(:,5),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,2)+ic8_v.vapor_data(:,4),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,2)-ic8_v.vapor_data(:,4),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,2)+tc_v.vapor_data(:,4),':','color',[0 0.5 0],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,2)-tc_v.vapor_data(:,4),':','color',[0 0.5 0],'linewidth',2)
hold on

legend('ic8(l)0.2e-3','ic8(l)2e-3','ic8(v)','3c(l)0.2e-3','3c(l)2e-3','3c(v)','Location','northeastoutside')
xlabel('Time aSOI [ms]')
ylabel('Width @15mm [mm]')
set(gcf,'color','w');
grid on
set(gca,'FontSize',18)
set(gca,'linewidth',1.5)
box on

figure(3)
plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,6),'--','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,6),':','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,5),'color',[0 0 1.0000],'linewidth',3)
hold on
plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,6),'--','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,6),':','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,5),'color',[0 0.5 0],'linewidth',3)
hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,6)+ic8_l_2e4.liquid_data_low_thr(:,9),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,6)-ic8_l_2e4.liquid_data_low_thr(:,9),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,6)+ic8_l_2e3.liquid_data_high_thr(:,9),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,6)-ic8_l_2e3.liquid_data_high_thr(:,9),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,6)+tc_l_2e4.liquid_data_low_thr(:,9),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,6)-tc_l_2e4.liquid_data_low_thr(:,9),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,6)+tc_l_2e3.liquid_data_high_thr(:,9),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,6)-tc_l_2e3.liquid_data_high_thr(:,9),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,5)+ic8_v.vapor_data(:,7),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,5)-ic8_v.vapor_data(:,7),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,5)+tc_v.vapor_data(:,7),':','color',[0 0.5 0],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,5)-tc_v.vapor_data(:,7),':','color',[0 0.5 0],'linewidth',2)
hold on

legend('ic8(l)0.2e-3','ic8(l)2e-3','ic8(v)','3c(l)0.2e-3','3c(l)2e-3','3c(v)','Location','northeastoutside')
xlabel('Time aSOI [ms]')
ylabel('Width @25mm [mm]')
set(gcf,'color','w');
grid on
set(gca,'FontSize',18)
set(gca,'linewidth',1.5)
box on

figure(4)
plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,10),'--','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,10),':','color',[0.3010 0.7450 0.9330],'linewidth',3)
hold on
plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,8),'color',[0 0 1.0000],'linewidth',3)
hold on
plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,10),'--','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,10),':','color',[0.4660 0.6740 0.1880],'linewidth',3)
hold on
plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,8),'color',[0 0.5 0],'linewidth',3)
hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,10)+ic8_l_2e4.liquid_data_low_thr(:,13),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e4.liquid_data_low_thr(:,1),ic8_l_2e4.liquid_data_low_thr(:,10)-ic8_l_2e4.liquid_data_low_thr(:,13),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,10)+ic8_l_2e3.liquid_data_high_thr(:,13),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(ic8_l_2e3.liquid_data_high_thr(:,1),ic8_l_2e3.liquid_data_high_thr(:,10)-ic8_l_2e3.liquid_data_high_thr(:,13),':','color',[0.3010 0.7450 0.9330],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,10)+tc_l_2e4.liquid_data_low_thr(:,13),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e4.liquid_data_low_thr(:,1),tc_l_2e4.liquid_data_low_thr(:,10)-tc_l_2e4.liquid_data_low_thr(:,13),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,10)+tc_l_2e3.liquid_data_high_thr(:,13),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(tc_l_2e3.liquid_data_high_thr(:,1),tc_l_2e3.liquid_data_high_thr(:,10)-tc_l_2e3.liquid_data_high_thr(:,13),':','color',[0.4660 0.6740 0.1880],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,8)+ic8_v.vapor_data(:,10),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(ic8_v.vapor_data(:,1),ic8_v.vapor_data(:,8)-ic8_v.vapor_data(:,10),':','color',[0 0 1.0000],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,8)+tc_v.vapor_data(:,10),':','color',[0 0.5 0],'linewidth',2)
% hold on
% plot(tc_v.vapor_data(:,1),tc_v.vapor_data(:,8)-tc_v.vapor_data(:,10),':','color',[0 0.5 0],'linewidth',2)
hold on

legend('ic8(l)0.2e-3','ic8(l)2e-3','ic8(v)','3c(l)0.2e-3','3c(l)2e-3','3c(v)','Location','northeastoutside')
xlabel('Time aSOI [ms]')
ylabel('Width @35mm [mm]')
set(gcf,'color','w');
grid on
set(gca,'FontSize',18)
set(gca,'linewidth',1.5)
box on