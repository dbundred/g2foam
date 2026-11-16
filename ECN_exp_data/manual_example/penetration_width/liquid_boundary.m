load('C:\Users\joohwan\Desktop\ECN7_data_release\data\ic8\G1\liquid_length\liquid_boundary_2e-4_womax.mat');
figure(1)

for i=1:1:400
   imshow(outline1_wo_outline{1,i})
   time=(i-20)*1/67200*1000;
   time=['aSOI:',num2str(time),'ms'];
   text(20,20,time,'Fontsize',16,'color','w'); 
   pause(0.1)
end