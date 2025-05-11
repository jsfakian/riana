function   file_exists=Save_data_check_if_exists(material,Ep,DATA,status, stri)


if status==1
    
 if ~(strcmp(material,'Si'))==1 
     
    cd (material)
    
 else
     cd (material)
     cd (material)
    
 end
    
 
              file_exists=0;
%               save(['RESULTS_' material '_material_substrate=' DATA.substrate '_wavelength=' ...
%                 num2str(1000*DATA.wavelength_nm) 'nm_delay=' ... 
%                '_Ep=' num2str(DATA.fluence_J_cm2) '_thickness=' num2str(DATA.film_thickness_nm) ...
%                '_Date=' datestr(now,'dd-mmm-yyyy-HHMMSS') 'nm.mat'],...
%                     'DATA','-mat');
%                 save(['RESULTS_dimension=' int2str(dimen) '_Ep=' num2str(Ep) '_'...
%                     material '_Date=' datestr(now,'dd-mmm-yyyy-HHMMSS')  '.mat'],...
%                     'DATA','-mat');

a=pwd;





cd ..


reflectivity_scaled=scaledata(DATA.intensity, min(DATA.reflectivity), max(DATA.reflectivity));
transmissivity_scaled=scaledata(DATA.intensity, min(DATA.transmissivity), max(DATA.transmissivity));
absorbed_scaled=scaledata(DATA.intensity, min(DATA.absorpt), max(DATA.absorpt));
intensity_scaled=scaledata(DATA.intensity, min(DATA.Te/1000), max(DATA.Te/1000));

cd (a)

if strcmp(DATA.substrate,'SiO2')==1
    title_substrate='SiO_2';
else 
    title_substrate=DATA.substrate;
    
end 



%%% Refletivity evolution 
fig=figure('Visible', stri);;
H1=axes;
plot(DATA.t,DATA.reflectivity,'r-',  'LineWidth', 2.5 );hold on;
if strcmp(DATA.material,'Mo')==1 || (strcmp(DATA.material,'Steel')==1 && ~(DATA.wavelength==1.026)) 
     legend('Reflectivity', 'FontSize', 14, 'Location', 'best')
else 
     plot(DATA.t,reflectivity_scaled,'k--','LineWidth', 1.5);
     legend('Reflectivity', 'Laser Pulse (a.u)', 'FontSize', 14, 'Location', 'best')
end 

xlabel('Time [sec]',  'FontName','Arial','FontSize',16); 
ylabel ('Reflectivity', 'FontName','Arial','FontSize',16)
  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')
  title([DATA.material '/' title_substrate])
xlim([-20*DATA.t(2),(1+0.1)*max(DATA.t)])
set(H1,'FontSize',15,'FontName','Arial');
set(gcf,'Color',[1,1,1])
saveas(fig, 'Reflectivity_png.png');
close(figure)
% saveas(fig, 'Reflectivity_pdf.pdf');

% 
% print(gcf, 'Reflectivity_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'Reflectivity_jpeg.jpeg', '-djpeg', '-r300')


%%% Transmissivity evolution 
fig=figure('Visible', stri);;
H1=axes;
plot(DATA.t,DATA.transmissivity,'r-',  'LineWidth', 2.5 );hold on;

if strcmp(DATA.material,'Mo')==1 || (strcmp(DATA.material,'Steel')==1 && ~(DATA.wavelength==1.026)) 
legend('Transmissivity', 'FontSize', 14, 'Location', 'best')
else 
plot(DATA.t,transmissivity_scaled,'k--','LineWidth', 1.5);
legend('Transmissivity', 'Laser Pulse (a.u)', 'FontSize', 14, 'Location', 'best')
end 



xlabel('Time [sec]',  'FontName','Arial','FontSize',16); 
ylabel ('Transmissivity', 'FontName','Arial','FontSize',16)
  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')
    title([DATA.material '/' title_substrate])
xlim([-20*DATA.t(2),(1+0.1)*max(DATA.t)])
set(H1,'FontSize',15,'FontName','Arial');
set(gcf,'Color',[1,1,1])
saveas(fig, 'Transmissivity_png.png');
close(figure)
% saveas(fig, 'Transmissivity_pdf.pdf');

% print(gcf, 'Transmissivity_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'Transmissivity_jpeg.jpeg', '-djpeg', '-r300')


%%% Absorptivity evolution 
fig=figure('Visible', stri);;
H1=axes;
plot(DATA.t,DATA.absorpt,'r-',  'LineWidth', 2.5 );hold on;


if strcmp(DATA.material,'Mo')==1 || (strcmp(DATA.material,'Steel')==1 && ~(DATA.wavelength==1.026)) 

legend('Absoprtivity', 'FontSize', 14, 'Location', 'best')
else 
plot(DATA.t,absorbed_scaled,'k--','LineWidth', 1.5);
legend('Absoprtivity', 'Laser Pulse (a.u)', 'FontSize', 14, 'Location', 'best')
end 




xlabel('Time [sec]',  'FontName','Arial','FontSize',16); 
ylabel ('Absoprtivity', 'FontName','Arial','FontSize',16)
  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')
    title([DATA.material '/' title_substrate])
xlim([-20*DATA.t(2),(1+0.1)*max(DATA.t)])
set(H1,'FontSize',15,'FontName','Arial');
set(gcf,'Color',[1,1,1])
saveas(fig, 'Absorptivity_png.png');
close(figure)
% saveas(fig, 'Absorptivity_pdf.pdf');
% print(gcf, 'Absorptivity_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'Absorptivity_jpeg.jpeg', '-djpeg', '-r300')



%%% Te and TL evolution on the first pixel 
fig=figure('Visible', stri);;
H1=axes;
plot(DATA.t,DATA.Te/1000,'r-',  'LineWidth', 2.5 );hold on;
plot(DATA.t, DATA.TL/1000, 'b', 'LineWidth', 1.5 );hold on;
plot(DATA.t,intensity_scaled,'k--','LineWidth', 1.5);
legend('T_e', 'T_L', 'Laser Pulse (a.u)', 'FontSize', 14, 'Location', 'best')
xlabel('Time [sec]',  'FontName','Arial','FontSize',16); 
ylabel ('Temperature [10^3 K]', 'FontName','Arial','FontSize',16)
  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')
    title([DATA.material '/' title_substrate ' (Surface of Stack)'])
xlim([-20*DATA.t(2),(1+0.1)*max(DATA.t)])
set(H1,'FontSize',15,'FontName','Arial');
set(gcf,'Color',[1,1,1]);
saveas(fig, 'Te_Tl_t_png.png');
close(figure)
% saveas(fig, 'Te_Tl_t_pdf.pdf');

% print(gcf, 'Te_Tl_t_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'Te_Tl_t_jpeg.jpeg', '-djpeg', '-r300')


% 

%%% Electron temmperature inside the first material
fig=figure('Visible', stri);;

H1=axes;
imagesc(DATA.t,DATA.L1_first, DATA.Te_t_z_upper/1000);colorbar;
 %set(gca,'YDir','normal')
a=colorbar;
  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')
set(H1,'FontSize',15,'FontName','Arial');
 L=ylabel(a,'Electron Temperature [10^3 K]');
 set(L,'FontSize',15);
 set(a,'FontSize',15);
% set(gca,'YDir','normal')
   title(['T_e in ' DATA.material])
xlabel('Time [s]', 'FontName','Arial','FontSize',16);
ylabel('Depth [m]','FontName','Arial','FontSize',16);
set(gcf,'Color',[1,1,1])
saveas(fig, 'Te_t_z_png.png');
close(figure)
% saveas(fig, 'Te_t_z_pdf.pdf');

% 
% print(gcf, 'Te_t_z_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'Te_t_z_jpeg.jpeg', '-djpeg', '-r300')

%%% Lattice temmperature inside the stack material
fig=figure('Visible', stri);;

H1=axes;
if max(DATA.TL_t_z(:)>1000)
imagesc(DATA.t,DATA.L_stack,DATA.TL_t_z/1000);colorbar;
else 
    imagesc(DATA.t,DATA.L_stack,DATA.TL_t_z);colorbar;

end 
 %set(gca,'YDir','normal')
a=colorbar;

  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')

set(H1,'FontSize',15,'FontName','Arial');

if max(DATA.TL_t_z(:)>1000)
 L=ylabel(a,'Lattice Temperature [10^3 K]');
 
else 
     L=ylabel(a,'Lattice Temperature [K]');
    
end
 set(L,'FontSize',15);
 set(a,'FontSize',15);
% set(gca,'YDir','normal')
   title(['T_L in ' DATA.material '/' title_substrate])
xlabel('Time [s]', 'FontName','Arial','FontSize',16);
ylabel('Depth [m]','FontName','Arial','FontSize',16);
hold on;plot([DATA.t(1) DATA.t(end)],[DATA.L1 DATA.L1],'w--', 'LineWidth', 2)
set(gcf,'Color',[1,1,1])

saveas(fig, 'TL_t_z_stack_png.png');
close(figure)
% saveas(fig, 'TL_t_z_stack_pdf.pdf');

% print(gcf, 'TL_t_z_stack_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'TL_t_z_stack_jpeg.jpeg', '-djpeg', '-r300')


%%% TL_upper_layer

fig=figure('Visible', stri);;

H1=axes;
if max(DATA.TL_t_z(:)>1000)
imagesc(DATA.t, DATA.L1_first, DATA.TL_t_z_upper/1000);colorbar;

else 
    imagesc(DATA.t, DATA.L1_first, DATA.TL_t_z_upper);colorbar;
end 
 %set(gca,'YDir','normal')
a=colorbar;
  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')

set(H1,'FontSize',15,'FontName','Arial');


if max(DATA.TL_t_z(:)>1000)
 L=ylabel(a,'Lattice Temperature [10^3 K]');
 
else 
     L=ylabel(a,'Lattice Temperature [K]');
    
end
 set(L,'FontSize',15);
 set(a,'FontSize',15);
% set(gca,'YDir','normal')
 
xlabel('Time [s]', 'FontName','Arial','FontSize',16);
ylabel('Depth [m]','FontName','Arial','FontSize',16);
  title(['T_L in ' DATA.material])
set(gcf,'Color',[1,1,1])


saveas(fig, 'TL_t_z_upper_png.png');
close(figure)
% saveas(fig, 'TL_t_z_upper_pdf.pdf');

% print(gcf, 'TL_t_z_upper_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'TL_t_z_upper_jpeg.jpeg', '-djpeg', '-r300')






%%%% Lattice temperature inside the depth of upper material-contour plot

fig=figure('Visible', stri);;
 H1=axes;
 if max(DATA.TL_t_z_upper(:))>1000
 [C,h]=contourf(DATA.t, DATA.L1_first, DATA.TL_t_z_upper/1000);colorbar;
 else 
      [C,h]=contourf(DATA.t, DATA.L1_first, DATA.TL_t_z_upper);colorbar;
 end 
 
 set(gca,'YDir','reverse');
  clabel(C, h, 'FontSize', 10, 'Color', 'w', 'LabelSpacing', 100);
 %set(gca,'YDir','normal')
a=colorbar;

  set(H1,'XMinorTick','On');set(H1,'YMinorTick','On')

set(H1,'FontSize',15,'FontName','Arial');

% hold on;plot([DATA.t(1) DATA.t(end)],[DATA.L1 DATA.L1],'w--', 'LineWidth', 2)
if max(DATA.TL_t_z_upper(:))>1000
 L=ylabel(a,'Lattice Temperature [10^3 K]');
 
else 
     L=ylabel(a,'Lattice Temperature [K]');
end 
 set(L,'FontSize',15);
 set(a,'FontSize',15);
% set(gca,'YDir','normal')
   title(['T_L in ' DATA.material])
xlabel('Time [s]', 'FontName','Arial','FontSize',16);
ylabel('Depth [m]','FontName','Arial','FontSize',16);
set(gcf,'Color',[1,1,1])


saveas(fig, 'TL_t_z_upper_layer_contour_png.png');
close(figure)
% saveas(fig, 'TL_t_z_upper_layer_contour_pdf.pdf');

% print(gcf, 'TL_t_z_upper_layer_contour_pdf.pdf', '-dpdf');  % -dpdf = device type PDF
% print(gcf,'TL_t_z_upper_layer_contour_jpeg.jpeg', '-djpeg', '-r300')

%%%% save in an XL file
% filename = ['Simulation_results ' material '.xlsx'];
% 
% writematrix(DATA.material, filename, 'Sheet', 'Material');
% writematrix(DATA.substrate, filename, 'Sheet', 'Substrate');
% writematrix(DATA.L1, filename, 'Sheet', 'Thickness of First Layer (m)');
% writematrix(DATA.L2, filename, 'Sheet', 'Thickness of Substrate (m)');
% writematrix(DATA.delay, filename, 'Sheet', 'Pulse delay (s)');
% writematrix(DATA.wavelength, filename, 'Sheet', 'Wavelength (m)');
% writematrix(DATA.Ep, filename, 'Sheet', 'Fluence (J per cm2)');
% writematrix(DATA.t, filename, 'Sheet', 'Time (sec)');
% writematrix([DATA.t, DATA.intensity], filename, 'Sheet', 'Time_Intensity');
% writematrix([DATA.t, DATA.Te,DATA.TL], filename, 'Sheet', 'Time_Te_TL');
% writematrix(DATA.L1_first, filename, 'Sheet', 'First Upper Layer depth (m)');
% writematrix(DATA.L_stack, filename, 'Sheet', 'Stack depth  (m)');
% writematrix(DATA.Te_t_z_upper, filename, 'Sheet', 'Te_upper_layer (t,z) (K)');
% writematrix(DATA.TL_t_z_upper, filename, 'Sheet', 'TL_upper_layer (t,z) (K)');
% writematrix(DATA.TL_t_z, filename, 'Sheet', 'TL_stack (t,z) (K)');


% Open a text file for writing
fileID = fopen('Conditions_materials.txt', 'w');
% writematrix(DATA.material, filename, 'Sheet', 'Material');
% writematrix(DATA.substrate, filename, 'Sheet', 'Substrate');
% writematrix(DATA.L1, filename, 'Sheet', 'Thickness of First Layer (m)');
% writematrix(DATA.L2, filename, 'Sheet', 'Thickness of Substrate (m)');
% writematrix(DATA.delay, filename, 'Sheet', 'Pulse delay (s)');
% writematrix(DATA.wavelength, filename, 'Sheet', 'Wavelength (m)');
% writematrix(DATA.Ep, filename, 'Sheet', 'Fluence (J per cm2)');

% Write the values to the file
fprintf(fileID, 'Material :  %s\n', DATA.material);
fprintf(fileID, 'Substrate :  %s\n', DATA.substrate);
fprintf(fileID, 'Thickness of First Layer (m) = %d\n', DATA.L1);
fprintf(fileID, 'Thickness of Substrate (m) = %d\n', DATA.L2);
fprintf(fileID, 'Pulse delay (s) = %d\n', DATA.delay);
fprintf(fileID, 'Wavelength (m) = %d\n', DATA.wavelength);
fprintf(fileID, 'Fluence (J/cm2) = %d\n',  DATA.fluence_J_cm2);
fprintf(fileID, 'Melting Temperature (Tmelt) (K) = %d\n',  DATA.Tmelt);
if   DATA.molten>0
fprintf(fileID, 'Thickness of Molten material (i.e. if Lattice Temperature in the layer is larger than Tmelt) (m) = %d\n',  DATA.molten);
end 
    
% Close the file
fclose(fileID);




A1=[DATA.t, DATA.intensity];
save('Time_Intensity.txt','A1','-ascii'); A1=[];
A1=[DATA.t, DATA.Te,DATA.TL];
save('Time_Te_TL.txt','A1','-ascii'); A1=[];
A1=DATA.L_stack;
save('Stack depth_(m).txt', 'A1', '-ascii');A1=[];
A1=DATA.Te_t_z_upper;
save('DATA.Te_t_z_upper.txt', 'A1', '-ascii');A1=[];
A1=DATA.TL_t_z_upper;
save('Te_upper_layer_(t,z)_(K).txt', 'A1', '-ascii');A1=[];
A1=DATA.TL_t_z_upper;
save('TL_upper_layer_(t,z)_(K).txt', 'A1', '-ascii');A1=[];
A1=DATA.TL_t_z;
save('TL_stack_(t,z)_(K).txt', 'A1', '-ascii');A1=[];

A1=[DATA.t, DATA.intensity, (DATA.transmissivity)'];
save('Time_Intensity_tranmissivity.txt','A1','-ascii'); A1=[];

A1=[DATA.t, DATA.intensity, (DATA.absorpt)'];
save('Time_Intensity_absorptivity.txt','A1','-ascii'); A1=[];

A1=[DATA.t, DATA.intensity, (DATA.reflectivity)'];
save('Time_Intensity_reflectivity.txt','A1','-ascii'); A1=[];



end