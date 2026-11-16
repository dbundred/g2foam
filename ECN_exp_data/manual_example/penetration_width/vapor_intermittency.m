load('C:\Users\joohwan\Desktop\ECN7_data_release\data\ic8\G1\vapor_length\Schlieren_probability.mat');
figure(1)

for i=1:1:400
   imagesc(Schlieren_process{1,i})
   axis image;
   colorbar
   caxis([0 100]);
   time=(i-21)*1/67200*1000;
   time=['aSOI:',num2str(time),'ms'];
   text(20,20,time,'Fontsize',16,'color','w'); 
   pause(0.1)
end