clear all;
warning off;
global n_2 k_2 n_1 k_1 stri


%% INPUT FOR CODE
%%%%%%%%%% THE FOLLOWING WILL BE ENETERED FROM USER
 Ep1=0.1; % J/cm^2;                                %%%%%%%%%%   TO BE INSERTED FROM USER
 wavelength1=1026; % nm;     %%%%%%%%%%   TO BE INSERTED FROM USER
 tp1=170;%fs                      %%%%%%%%%%   TO BE INSERTED FROM USER
 t_delay1=0; %fs          %%%%%%%%%%   TO BE INSERTED FROM USER
 t_max1=20; %ps               %%%%%%%%%%   TO BE INSERTED FROM USER
 L1=50; %nm                 %%%%%%%%%%   TO BE INSERTED FROM USER
 material='Cr';              % Ni, Ti, Au, Cu, Ag, Al, Steel, Cr, W, Pt, Mo    %%%%%%%%%%   TO BE INSERTED FROM USER
 material_substrate='SiO2';  % SiO2, Si, Soda               %%%%%%%%%%   TO BE INSERTED FROM USER
 %%%%%%%% END OF INPUT VARIABLES FROM USER

 stri='on'; % if you dont wish the plots to be displayed, it should be 'off' on the interface
 
 
 
 %% SUBSTRATE 
 
 if strcmp(material_substrate,'Si')==1 

      if wavelength1/1000==0.800
         n_2=3.6750; %Green
         k_2=0.0054113;
         
     elseif wavelength1/1000==1.026
         n_2=3.5632; %Green
         k_2=0.00027806;
         
     elseif wavelength1/1000==0.515
         
         n_2=4.2170; %Green
         k_2=0.037;
         
   
     elseif wavelength1/1000==0.248
         
         n_2=1.57;%Aspnes
         k_2=3.5650;  
         
      else 
          
          n_2=1;%%%%%%%%%%   TO BE INSERTED FROM USER
          k_2=1;%%%%%%%%%%%%%%% TO BE INSERTED FROM USER
         
      end 
     
 elseif strcmp(material_substrate,'Soda')==1   

        n_2=1.5130-0.003169*(wavelength1/1000)^2+0.003962/(wavelength1/1000)^2;
        k_2=0;
    
       
 elseif strcmp(material_substrate,'SiO2')==1
    
     wavelength=wavelength1*1e-3; % in mum
    term1=0.6961663*wavelength^2/(wavelength^2-0.0684043^2);
    term2=0.4079426*wavelength^2/(wavelength^2-0.1162414^2);
   term3=0.8974794*wavelength^2/(wavelength^2-9.896161^2);
   n_2=sqrt(1+term1+term2+term3);
   k_2=0;
      
   
          
 elseif strcmp(material_substrate,'Air')==1
   
   n_2=1;
   k_2=0;
     
 end 
 

 
 
 %% UPPER MATERIAL (METAL)
         if strcmp(material,'Mo')==1 && wavelength1==1026
          n_1= 2.4357; 
           k_1=4.1672;  
           
         elseif strcmp(material,'Mo')==1 && ~(wavelength1==1026)
         
          n_1= 2.4357;                  %%%%%%%%%%   TO BE INSERTED FROM USER
          k_1=4.1672;                  %%%%%%%%%%   TO BE INSERTED FROM USER
        
         end 
 
 
         if strcmp(material,'Steel')==1   
         
             if ~(wavelength1==1026)
                n_1=2.9602;                  %%%%%%%%%%   TO BE INSERTED FROM USER
                k_1=3.8990;                  %%%%%%%%%%   TO BE INSERTED FROM USER

             end 
         
         end 
         
         
         
   %% RUN THE CODES      
 
 if  ((strcmp(material,'Au')==1|| strcmp(material,'Cu')==1||strcmp(material,'Cr')==1 ...
         || strcmp(material,'Ag')==1  || strcmp(material,'W')==1 || strcmp(material,'Al')==1 || strcmp(material,'Mo')==1)&& t_delay1<=5000)
     

        Code_for_all_metals_but_two(Ep1,wavelength1, tp1, t_delay1, t_max1, material, material_substrate, L1);
  
     
 elseif (strcmp(material,'Ni')==1 || strcmp(material,'Steel')==1 || strcmp(material,'Ti')==1 || strcmp(material,'Pt')==1 )||...
         ((strcmp(material,'Au')==1|| strcmp(material,'Cu')==1||strcmp(material,'Cr')==1 ...
         || strcmp(material,'Ag')==1  || strcmp(material,'W')==1 || strcmp(material,'Al')==1 || strcmp(material,'Mo')==1)&& t_delay1>5000)
     
     
  
        Code_for_two(Ep1,wavelength1, tp1, t_delay1, t_max1, material, material_substrate, L1);
        
 else 
     
 end 



