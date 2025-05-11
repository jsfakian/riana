function Code_for_all_metals_but_two(Ep1,wavelength1, tp1, t_delay1, t_max1, material, material_substrate, L1)

%%

global wavelength1 Ae BL n_2 k_2 stri Tmelt n_1 k_1 alpha1 break_code T_cr


%% Define Energy pulse, laser wavelength, pulse duration, thickness of the
%%% material

% Ep1=0.3; % J/cm^2;
% wavelength1=1026; % nm; 
% tp1=170;%fs 
% t_delay1=0*300; %fs 
% t_max1=20; %ps 
% material='Ti'; 
% material_substrate='SiO2'; 
L1=L1; % nm
 L1=L1*1e-9 ;%m  
L2=100e-9; %m
L=L1+L2; %m



%% Laser parameters

Ep = Ep1*1e4 ;  % Energy density (J/m^2)
wavelength1=wavelength1; %nm Laser wavelength
wavelength=wavelength1*1e-9; % in meters
tp = tp1*1e-15;  % Pulse duration (s)
t_delay=t_delay1*1e-15; % Pulse delay(s)
t_max = t_max1*1e-12;  % Simulation time (s)
break_code=0; % label to ensure that if Te>50000 it does not continue to run 



%% Spatial and temporal discretization



 dx1=0.5e-9;
 l1=0:dx1:L1;
 N1=round(length(l1));


dx1 = L1/(N1-1);% Spatial step 
% dx2=  L2/(N2-1);  % Spatial step
dx2=dx1;

second_layer=0:dx2:L2;
N2=length(second_layer);
num_x=N1+N2;


DEPTH=[0:dx1:(N1-1)*dx1,(N1-1)*dx1+dx2:dx2:L2+L1];


%%% boundaries of the two materials
first_material=0:dx1:(N1-1)*dx1;
second_material=(N1-1)*dx1+dx2:dx2:L1+L2;


if t_max1<=1000
    number_t=500;
else 
        number_t=t_max1;
       
end 
 tspan=linspace(0, t_max,number_t);
 % tspan=[0, t_max];

%% Initialize Temperature Field

Te=300*ones(num_x,1);  % Electron Temperature profile
TL=300*ones(num_x,1);  % Lattice Temperature profile
T0=[Te;TL];          % Temperatures
%%%%%%%%%%%%%%%%%%




%% Solve using ode15s or ode23tb (stiff solver)
% 
%    options = odeset('OutputFcn', @myOutputFcn);
% 
%    options = odeset('RelTol',1e-8,'AbsTol',1e-10);



% Use an anonymous wrapper to pass extra args
options = odeset('Events', @stop_condition);

  if t_max/tp<3000
    
[t,y]=ode23tb(@(t,y) heat_diffusion(t,y,Ep,dx1,dx2, num_x, ...
    material,material_substrate, N1, N2, wavelength,L1,L2, tp, t_delay), tspan, T0, options);
 else 
 [t,y]=ode15s(@(t,y) heat_diffusion(t,y,Ep,dx1,dx2, num_x, ...
     material,material_substrate, N1, N2, wavelength,L1,L2, tp, t_delay), tspan, T0, options);
  end 
  
 

  
  
  
 %% Extract Solution
 % row: time, column: position 
 Te_sol_1=y(:,1:N1);
 Te_sol_2=y(:,N1+1:N1+N2); % 0=L1
 TL_sol_1=y(:,N1+N2+1:N1+N2+N1); %[0 L1]
 TL_sol_2=y(:,N1+N2+N1+1:2*(N1+N2)); %[L1 L2]
 
 A=[TL_sol_1(:,1:end)';TL_sol_2(:,2:end)'];

%%%% to check wherther Te>50000 K. if this the case,it will be mentioned in the 'Conditions' file  
     if size(Te_sol_1,1)< length(tspan)
                    break_code=1;
         
     end 
 %%%%
 
 %%% calculation of depth of molten volume 
 [r,c]=find(TL_sol_1>=Tmelt);
 if length(r)>0
 maximum_depth_pixel=max(c);
depthofmoltenpart=first_material(maximum_depth_pixel);
 else
     depthofmoltenpart=0;
 end 
 
  %%% calculation of depth of ablated volume 
 [r,c]=find(TL_sol_1>=T_cr);
 if length(r)>0
 maximum_depth_pixel_ablated=max(c);
depthofablatedpart=first_material( maximum_depth_pixel_ablated);
 else
     depthofablatedpart=0;
 end 
 

%%  save data in an excel file

% Calculate Intensity at all t
Icoeff = 1/2*exp(-4 * log(2) * (t-3*tp).^2 / tp^2)+...
   1/2 * exp(-4 * log(2) * (t-3*tp-t_delay).^2 / tp^2);

% calculate Optical parameters at all t
 [n_2,k_2, alpha2]=substrate_optical_parameters(material_substrate,  n_2, k_2);
 
 
 N_1=[];K_1=[];
 
 if ~(strcmp(material,'Mo')==1)
for i=1:length(t)
 [n_1,k_1, alpha1]=material_optical_parameters(material,Te_sol_1(i,1),TL_sol_1(i,1), Ae, BL, N1, wavelength);
 N_1=[N_1,n_1];
  K_1=[K_1,k_1];
end 

 else 
     for i=1:length(t)
 
 N_1=[N_1,n_1];
  K_1=[K_1,k_1];
     end 
 end 


[absorpt, reflectivity, transmissivity]=optical_parameters_plot(N_1,K_1, n_2, k_2, L1, wavelength);




                        DATA.material=material;
                        DATA.substrate=material_substrate;
%                          DATA.refractive_index_material=material;
%                         DATA.refractive_index_substrate=material_substrate;
                        DATA.fluence_J_cm2=Ep1;
                        DATA.wavelength=wavelength;
                        DATA.pulse_duration_fs=tp1;
                        DATA.delay=t_delay;
                        DATA.dz_nm=dx1;
                        DATA.intensity=Icoeff;
                        DATA.Te=Te_sol_1(:,1);
                        DATA.TL=TL_sol_1(:,1);
                        DATA.t=t;
                        DATA.L1=L1;
                        DATA.L2=L2;
                        DATA.L1_first=0:dx1:(N1-1)*dx1;
                        DATA.L_stack=[0:dx1:(N1-1)*dx1,(N1-1)*dx1+dx2:dx2:L2+L1];
                        DATA.Te_t_z_upper=Te_sol_1';
                        DATA.TL_t_z=A;
                        DATA.TL_t_z_upper=TL_sol_1(:,1:end)';
                        DATA.reflectivity=reflectivity;
                        DATA.absorpt=absorpt;
                        DATA.transmissivity=transmissivity;
                        DATA.molten=depthofmoltenpart;
                        DATA.Tmelt=Tmelt;
                        DATA.Tablation=0.95*T_cr;
                        DATA.ablated=depthofablatedpart;

% save('RESULTS.mat','DATA','-mat')

 [status,msg,msgID]=mkdir(material);


 
 file_exists=Save_data_check_if_exists_test(material,Ep1,DATA,status, stri);
         
      cd('..')


      
end 

%% 

function [dTdt]=heat_diffusion(t,y,Ep,dx1,dx2, num_x, material, material_substrate, ...
    N1, N2, wavelength,L1,L2, tp, t_delay)
global Ae BL n_2 k_2 Tmelt n_1 k_1 T_cr

dTdt=zeros(2*num_x,1);



% extract electron and lattice temperatures 
Te1=y(1:N1,1);
Te2=y(N1+1:N1+N2,1);
TL1=y(N1+N2+1:N1+N2+N1,1);
TL2=y(N1+N2+N1+1:N1+N2+N1+N2,1);





%%%%%%%%%%%%%%%%%
%%%% Material thermophysical properties and electron phonon coupling and
%%%% electron heat capaicty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
[rho_s,rho_l,CL_s,CL_l,Tmelt,Ae,BL,ke0_1,T_cr] = ...
    material_parameters_constant(material);
 ke0=ke0_1; %J/(m*sec*K)
 CL_s=CL_s*rho_s;
 CL_l=CL_l*rho_l;
 
 
 AB_new=Ae/BL*ones(size(Te1,1),1); %1/K 
 
 

if strcmp(material,'Ni')==1 

 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_ni(1); % calculation of Ce and G for material
 
elseif strcmp(material,'Au')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Au(1); % calculation of Ce and G for material
 
 elseif strcmp(material,'Steel')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Fe(1); % calculation of Ce and G for material
 
  elseif strcmp(material,'Ti')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_ti(1); % calculation of Ce and G for material
 
  elseif strcmp(material,'W')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_W(1); % calculation of Ce and G for material
 
   elseif strcmp(material,'Cu')==1
    

 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_cu(1); % calculation of Ce and G for material
 
   elseif strcmp(material,'Ag')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Ag(1); % calculation of Ce and G for material
 
    elseif strcmp(material,'Al')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Al(1); % calculation of Ce and G for material
 
     elseif strcmp(material,'Pt')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Pt(1); % calculation of Ce and G for material
 
      elseif strcmp(material,'Mo')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Mo(1); % calculation of Ce and G for material
 
 
   elseif strcmp(material,'Cr')==1
 
 [p4,T_max_e]=ce_coupl_cr(1); % calculation of Ce and G for material
 
end


    
   
 %%%%%


 %%%%%%%%%%%%%%%%%% second layer
 
 if strcmp(material_substrate,'SiO2')==1 || strcmp(material_substrate,'Si')==1 ||  strcmp(material_substrate,'Soda')==1   
    
 [rho_s2,rho_l2,CL_s2,CL_l2,Tmelt2,Ae2,BL2,ke0_12,T_cr2] = material_parameters_constant(material_substrate);
 kL2=ke0_12; %J/(m*sec*K)
 CL_s2=CL_s2*rho_s2;
 CL_l2=CL_l2*rho_l2; % not conditions that lead to phase transitions
 
 elseif  strcmp(material_substrate,'Air')==1
 kL2=0;
 CL_s2=100000000000000;
  CL_l2=100000000000000;
 

 end

 
 % optical parameters for the substrate
 
 [n_2,k_2, alpha2]=substrate_optical_parameters(material_substrate,  n_2, k_2);

%

if ~(strcmp(material,'Mo'))==1

[n_1,k_1, alpha1]=material_optical_parameters(material,Te1,TL1, Ae, BL, N1, wavelength);

else 
    n_1=n_1;
    k_1=k_1;
 alpha1=(4*pi)*k_1/wavelength*ones(size(Te1)); %1/m

end 

   %Calculate the optical paramters


 %%%%%%%%%%%%%%%%
 %%%% COMPUTATION OF REFLECTIVITY AND TRANSMISSIVITY FOR THIN FILMS: ONE
 %%%% THIN FILM IS SANDWICHED BETWEEN AIR AND SILICON/SILICA
 

 
 REF_0=1+0*1j; %complex refractive index for material : air
 REF_1=n_1+k_1*1j; %complex refractive index for material 1
 REF_2=n_2+k_2*1j; %complex refractive index for material 2 (substrate)
 
  % A. Reflectivity

 r_01=(REF_1-REF_0)/(REF_1+REF_0);
 r_12=(REF_2-REF_1)/(REF_2+REF_1);
 

 beta_distance=2*pi*L1/wavelength*(REF_1);
 
 reflectivity_complex=(r_01+r_12*exp(2*beta_distance*1j))/...
     (1+r_01*r_12*exp(2*beta_distance*1j));
 
 reflectivity=reflectivity_complex*conj(reflectivity_complex);
 

 %%%%%%%%%%%%%%%%%%%%%
 % B. Transmissivity  
 
t_01=(2*REF_0)/(REF_1+REF_0);
 t_12=(2*REF_1)/(REF_2+REF_1);
 
 transmissivity_complex=(t_01*t_12*exp(1*beta_distance*1j))/...
     (1+r_01*r_12*exp(2*beta_distance*1j));
 
 transmissivity=transmissivity_complex*conj(transmissivity_complex)*real(REF_2)/real(REF_0);
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Energy absoprtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

absorpt=1-reflectivity-1*transmissivity;

if strcmp(material,'Ni')==1 |  strcmp(material,'Ti')==1 | strcmp(material,'Steel')==1 | strcmp(material,'Cr')==1  | strcmp(material,'Pt')==1 
 
    ballistic=1/alpha1(1);
  
    
else 
    
    if strcmp(material,'Cu')==1
        
            ballistic=70e-3*1e-6; %m
        
    elseif strcmp(material,'Au')==1
          ballistic=100e-3*1e-6; %m 
        
    elseif strcmp(material,'Al')==1
           ballistic=46e-3*1e-6; %m
        
    elseif strcmp(material,'Ag')==1
        
           ballistic=142e-3*1e-6; %m 
           
             elseif strcmp(material,'W')==1
        
           ballistic=0e-3*1e-6; %m 
           
                  elseif strcmp(material,'Mo')==1
        
           ballistic=20e-3*1e-6; %m 
    end 
    
end 



enhanced=1/(1/alpha1(1)+1*ballistic);
coefent=1-1*exp(-L1*enhanced);

   
% Calculate the heat source that attenuates
Intensity=zeros(2*(N1+N2),1);


Icoeff = 2/ (tp * sqrt(pi)) * sqrt(log(2)) * Ep/2 * exp(-4 * log(2) * (t-3*tp).^2 / tp^2)+...
    2/ (tp * sqrt(pi)) * sqrt(log(2)) * Ep/2 * exp(-4 * log(2) * (t-3*tp-t_delay).^2 / tp^2);

%for i=1:N1
 %   Intensity(i)=exp(-alpha1(1)*(i-1)*dx1);
    
%end



Source1=enhanced*exp(-enhanced*(0:dx1:L1)')/coefent*Icoeff*absorpt;



% second film

second_material_matrix=(N1-1)*dx1:dx2:L1+L2;


 %Source2=1*Source1(end)/enhanced*exp(alpha2(1)*(L1-second_material_matrix))*alpha2(1);



if ~(strcmp(material_substrate,'Air')==1)
Source2=1*transmissivity*Icoeff*exp(alpha2(1)*(L1-second_material_matrix))*alpha2(1);
else
    Source2=0*transmissivity*Icoeff*exp(alpha2(1)*(L1-second_material_matrix))*alpha2(1);
    
end 








        % heat capacity lattice first material
         CL_L= CL_l*ones(size(TL1));
        CL_S= CL_s*ones(size(TL1));
        CL_S(find(TL1>=Tmelt))=CL_l;
        
        
        % 2 Heat conductivity of electron first material
      
        if ~(strcmp(material,'Cr')==1)
        ke1=ke0*Te1./(AB_new.*Te1.^2+TL1);
        else 
            ke1=ke0*ones(size(Te1));
        end 
        
        % 2 Heat conductivity of lattice first material
         kL1=1*0.01*ke1;
       
        % 3 Electron-phono coupling and heat capacity of electron -first
        % material
        coupl_const1=ones(N1,1);
         CE=ones(N1,1);
         
         if ~(strcmp(material,'Cr')==1)
          q1=find(Te1<=threshold_ce_coupl);
         q2=find(Te1>threshold_ce_coupl);
                  if length(q1)>0
                    coupl_const(q1)=polyval(p3,Te1(q1));% to convert from W/(m^3*K) to  W/(microns^3*K)
                       CE(q1)=polyval(p1,Te1(q1)); % to convert from J/(m^3*K) to  J/(microns^3*K);
                  end
   
                   if q2>0
                    coupl_const(q2)=polyval(p4,Te1(q2));% to convert from W/(m^3*K) to  W/(microns^3*K)
                    CE(q2)=polyval(p2,Te1(q2)); % to convert from J/(m^3*K) to  J/(microns^3*K);
                    end
         elseif strcmp(material,'Cr')==1
              CE=5.8e4/300*Te1; 
              coupl_const=polyval(p4,Te1);
         end 
        
 % compute TTM       
         % Finite Difference: Second Order Central Scheme
    dTe1_dt = zeros(N1,1);
       dTL1_dt = zeros(N1,1);
            dTe2_dt = zeros(N2,1);
       dTL2_dt = zeros(N2,1);

            
       
numerator=kL1(N1-1)*TL1(N1-1)/dx1+kL2*TL2(2)/dx2;
denominator=kL1(N1-1)/dx1+kL2/dx2;
interface_temperature=numerator/denominator;
TL2(1)=interface_temperature;
TL1(N1)=interface_temperature;
      

        
        
if ~(strcmp(material,'Cr')==1)
    ke1(N1)=ke0*Te1(N1)./(AB_new(1).*Te1(N1).^2+TL1(N1));
else 
    ke1(N1)=ke0;
end 
         kL1(N1)=1*0.01*ke1(N1);
         
         
         if ~(strcmp(material,'Cr')==1)
          q1=find(Te1(N1)<=threshold_ce_coupl);
         q2=find(Te1(N1)>threshold_ce_coupl);
                  if length(q1)>0
                    coupl_const(q1)=polyval(p3,Te1(q1));% to convert from W/(m^3*K) to  W/(microns^3*K)
                       CE(q1)=polyval(p1,Te1(q1)); % to convert from J/(m^3*K) to  J/(microns^3*K);
                  end
   
                   if q2>0
                    coupl_const(q2)=polyval(p4,Te1(q2));% to convert from W/(m^3*K) to  W/(microns^3*K)
                    CE(q2)=polyval(p2,Te1(q2)); % to convert from J/(m^3*K) to  J/(microns^3*K);
                    end
         elseif strcmp(material,'Cr')==1
              CE(N1)=5.8e4/300*Te1(N1); 
              coupl_const(N1)=polyval(p4,Te1(N1));
         end 
         
         
         
         
         
         
         
               
           
 %%%%% midpoint calculations
      
       % Interior Points  % material 1
    
         for   i=2:N1-1

              ke_iph = 1/2 * (ke1(i)) + ke1(i+1);
              ke_imh = 1/2 * (ke1(i)) + ke1(i-1);
              
        kL_iph = 1/2 * (kL1(i)) + kL1(i+1);
         kL_imh = 1/2 * (kL1(i)) + kL1(i-1);
        
            
       dTe1_dt(i)=((ke_iph*(Te1(i+1)-Te1(i))-ke_imh*(Te1(i)-Te1(i-1)))/(dx1)^2-...
           coupl_const(i).*(Te1(i)-TL1(i))+...
           Source1(i))./CE(i);
       dTL1_dt(i)=((kL_iph*(TL1(i+1)-TL1(i))-kL_imh*(TL1(i)-TL1(i-1)))/(dx1)^2 +...
           coupl_const(i).*(Te1(i)-TL1(i)))./CL_S(i);
      
        end 
        
        
       % upper material, upper surface , at i=1 
       i=1;
             
       ke_iph = 1/2 * (ke1(i)) + ke1(i+1);
        kL_iph = 1/2 * (kL1(i)) + kL1(i+1);
        
dTe1_dt(1)=(ke_iph*(Te1(2)-Te1(i))/dx1^2- coupl_const(i).*(Te1(i)-TL1(i))+...
           Source1(i))./CE(i);
dTL1_dt(1)=(kL_iph*(TL1(2)-TL1(i))/dx1^2+ coupl_const(i).*(Te1(i)-TL1(i)))./CL_S(i);;
%  if strcmp(material,'Ni')==1 || strcmp(material,'Steel')==1
%    dTe1_dt(1)=dTe1_dt(2);
%         dTL1_dt(1)=dTL1_dt(2);
%  end 
%           
          %upper material, back surface, at i=N1
          
          % to ensure that theta_Ke/Theta_e=0;theta_Te/Theta_e=0
      
          i=N1;
           ke_imh = 1/2 * (ke1(i) + ke1(i-1));

           
           
         dTe1_dt(i)=(-ke_imh*(Te1(i)-Te1(i-1))/dx1^2-...
            coupl_const(i).*(Te1(i)-TL1(i))+...
           Source1(i))./CE(i);
       
    
       
   
        
       
           % Interior Points  % material 2
 for i=2:N2-1
    
     
       dTL2_dt(i)=1*(kL2*(TL2(i+1)+1*TL2(i-1)-2*TL2(i))/dx2^2+1*Source2(i))/CL_s2;
 end 
       
        % material 2: back end
        i=N2;
        
         

 dTL2_dt(i)=(-kL2*(TL2(i)-TL2(i-1))/dx2^2+1*Source2(i))./CL_s2;
       
       


       
       
% %        % interface between the two materials



i=N1;

%     A= ((kL1(i)+kL1(i-1))/2*(TL1(i) - TL1(i-1))/dx1 - 0*kL2*(TL2(2) - TL2(1))/dx2) / ((dx1 + dx2)/2);
%    A=A+1*coupl_const(i).*(Te1(i)-TL1(i));
%   
% 
 dTL1_dt(N1)=(dTL1_dt(N1-1)+dTL2_dt(2))/2;
%      dTL1_dt(N1)=A/CL_S(N1);
 dTL2_dt(1)=dTL1_dt(N1);

       
       % temporal gradient matrix
       dTdt=[ dTe1_dt;dTe2_dt;dTL1_dt;dTL2_dt];
       
 


end 


  %%
function  [n_2,k_2, alpha2]=substrate_optical_parameters(material_substrate,  n_2, k_2)
global wavelength1

wavelength=wavelength1*1e-3; % in mum
if strcmp(material_substrate,'SiO2')==1
    
    
    term1=0.6961663*wavelength^2/(wavelength^2-0.0684043^2);
    term2=0.4079426*wavelength^2/(wavelength^2-0.1162414^2);
   term3=0.8974794*wavelength^2/(wavelength^2-9.896161^2);
   n_2=sqrt(1+term1+term2+term3);
   k_2=0;
%  
% 
  elseif  strcmp(material_substrate,'Si')==1
      
        n_2=n_2;
      k_2=k_2;

  elseif  strcmp(material_substrate,'Si')==1
      
      n_2=n_2;
      k_2=k_2;

     
 elseif strcmp(material_substrate,'Soda')==1

     n_2=n_2;
     k_2=k_2;
       
end
 alpha2=4*pi*k_2/(wavelength*1e-6); % 1/m

end 
%%
function [n_1,k_1, alpha]=material_optical_parameters(material,Te1,TL1,Ae,BL,N1, wavelength)



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Material optical parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(material,'Ni')==1 | strcmp(material,'Au')==1|...
   strcmp(material,'Steel')==1| strcmp(material,'Ti')==1 ...
        | strcmp(material,'Cu')==1 | strcmp(material,'Ag')==1 | strcmp(material,'Cr')==1 ...
    | strcmp(material,'Al')==1   | strcmp(material,'W')==1 | strcmp(material,'Pt')==1
 
     
AB_new=Ae/BL*ones(N1,1); %1/K 


if strcmp(material,'Ni')==1 | strcmp(material,'Au')==1| ...
        strcmp(material,'Ti')==1| strcmp(material,'Cu')==1 | strcmp(material,'Ag')==1 ...
        | strcmp(material,'Cr')==1 | strcmp(material,'Al')==1 | strcmp(material,'W')==1  | strcmp(material,'Pt')==1

[reflect1,alpha,e1,e2]= LD_Tin_new(wavelength,material,'LD',Te1,TL1,AB_new, BL);  

alpha=alpha/1e-6; %1/m


%reflectivity=reflect1(2); %reflectivity of material;
penetration_depth=1/alpha(1);% in m  


e_1=e1(1);
e_2=e2(1);

 n_1=sqrt((e_1+sqrt(e_1^2+e_2^2))/2);
  k_1=sqrt((-e_1+sqrt(e_1^2+e_2^2))/2);
 

elseif  strcmp(material,'Steel')==1
    
      [pAaa,pBaa]=reflectivity_absorption_Fe(1);

    reflectivity=polyval(pAaa,Te1);
     penetration_depth=1e-3*polyval(pBaa,Te1)*1e-6;%m

     alpha=1./penetration_depth;
 

     k_1=alpha(1)*wavelength/(4*pi);

    trion_A=reflectivity(1)-1;
    trion_B=2*reflectivity(1)+2;
    trion_C=reflectivity(1)-1-k_1^2+reflectivity(1)*k_1^2;
    n_1=(-trion_B+sqrt(trion_B^2-4*trion_A*trion_C))/(2*trion_A);



end   
 
end
end 

function [absorpt, reflectivity, transmissivity]=optical_parameters_plot(N_1,K_1, n_2, k_2, L1, wavelength)

 REF_0=1+0*1j; %complex refractive index for material : air
 REF_1=N_1+K_1*1j; %complex refractive index for material 1
 REF_2=n_2+k_2*1j; %complex refractive index for material 2 (substrate)
 
  % A. Reflectivity
 
 r_01=(REF_1-REF_0)./(REF_1+REF_0);
 r_12=(REF_2-REF_1)./(REF_2+REF_1);

 beta_distance=2*pi*L1./wavelength*(REF_1);
 
 reflectivity_complex=(r_01+r_12.*exp(2*beta_distance*1j))./...
     (1+r_01.*r_12.*exp(2*beta_distance*1j));
 
 reflectivity=reflectivity_complex.*conj(reflectivity_complex);


 %%%%%%%%%%%%%%%%%%%%%
 % B. Transmissivity  
 
t_01=(2*REF_0)./(REF_1+REF_0);
 t_12=(2*REF_1)./(REF_2+REF_1);
 
 transmissivity_complex=(t_01.*t_12.*exp(1*beta_distance*1j))./...
     (1+r_01.*r_12.*exp(2*beta_distance*1j));
 
 transmissivity=transmissivity_complex.*conj(transmissivity_complex)*real(REF_2)/real(REF_0);
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Energy absoprtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

absorpt=1-reflectivity-1*transmissivity;


end 

function [value, isterminal, direction] = stop_condition(t, y)


    threshold = 50000;
    value = y(:,1) - threshold;   % Stop when y(1) crosses 5000

      isterminal=1;    % Stop the integration
    

    direction = 0;              % Detect all directions
end