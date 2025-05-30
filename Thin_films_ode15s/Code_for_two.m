function Code_for_two(Ep1,wavelength1, tp1, t_delay1, t_max1, material, material_substrate, L1)

%%

global  Ae BL n_2 k_2 n_1 k_1 stri Tmelt  alpha1 break_code T_cr


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
 L1d=L1; %nm
 L1=L1*1e-9 ;%m  
L2=100e-9; %m
L=L1+L2; %m

%% initialize parameters 

Te1_matrix=[];
TL1_matrix=[];
TL2_matrix=[];
t=0; % start time
Ta=[];
TE=[];
TL=[];
AA=0;
REFL=[];
ABSO=[];
TRANS=[];
 break_code=0; % label to ensure that if Te>50000 it does not continue to run 


%% Laser parameters

Ep = Ep1*1e4 ;  % Energy density (J/m^2)
wavelength1=wavelength1; %nm Laser wavelength
wavelength=wavelength1*1e-9; % in meters
tp = tp1*1e-15;  % Pulse duration (s)
t_delay=t_delay1*1e-15; % Pulse delay(s)
t_max = t_max1*1e-12;  % Simulation time (s)


%% Spatial and temporal discretization


% discretization size
if L1d<=20
 dx1=0.5e-9;
else 
dx1=1e-9;
end

 l1=0:dx1:L1;
 N1=round(length(l1));


dx1 = L1/(N1-1);% Spatial step 
% dx2=  L2/(N2-1);  % Spatial step 
dx2=dx1;
second_layer=0:dx2:L2;
N2=length(second_layer);



DEPTH=[0:dx1:(N1-1)*dx1,(N1-1)*dx1+dx2:dx2:L2+L1];



%%% boundaries of the two materials
first_material=0:dx1:(N1-1)*dx1;

second_material=(N1-1)*dx1+dx2:dx2:L1+L2;



if t_max1<=1000
    number_t=1000;
else 
        number_t=t_max1;
end 
 tspan=linspace(0, t_max,number_t);
 % tspan=[0, t_max];

%% Initialize Temperature Field

        Te1=300*ones(N1+1,1);  % Electron Temperature profile
TL1=300*ones(N1+1,1);  % Lattice Temperature profile
TL2=300*ones(N2+1,1);  % Lattice Temperature profile

%%%%%%%%%%%%%%%%%%




%% Solve using manual model
% 

    %%%%%%%%%%%%%%%%%
%%%% Material thermophysical properties and electron phonon coupling and
%%%% electron heat capaicty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
[rho_s,rho_l,CL_s,CL_l,Tmelt,Ae,BL,ke0_1,T_cr] = material_parameters_constant(material);
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
 
   elseif strcmp(material,'Cu')==1
    

 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_cu(1); % calculation of Ce and G for material
 
      elseif strcmp(material,'Pt')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Pt(1); % calculation of Ce and G for material
  
    elseif strcmp(material,'Al')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Al(1); % calculation of Ce and G for material
 
   elseif strcmp(material,'Ag')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Ag(1); % calculation of Ce and G for material
 
   elseif strcmp(material,'W')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_W(1); % calculation of Ce and G for material
 
       elseif strcmp(material,'Mo')==1
 
 [threshold_ce_coupl,p1,p2,p3,p4,T_max_e]=ce_coupl_Mo(1); % calculation of Ce and G for material
 
 
   elseif strcmp(material,'Cr')==1
 
 [p4,T_max_e]=ce_coupl_cr(1); % calculation of Ce and G for material
 
end






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
    
 
 %%%% start a loop to compute temperatures for electron and latice system.
 %%%% Similarly for the substate

while t<=t_max
    AA=AA+1;
    
    
    if t==0
    
  Te1=300*ones(N1+1,1);  % Electron Temperature profile
TL1=300*ones(N1+1,1);  % Lattice Temperature profile
TL2=300*ones(N2+2,1);  % Lattice Temperature profile
        
    else
      %%%% first layer
            Te1=[Te1(1,1); Te1];
            TL1=[TL1(1,1); TL1];
            
            if max(Te1(:))>=50000
                    break_code=1;
            break; 
            end 
 %%%% second layer

            TL2=[TL2;TL2(end)];
    end 
    
    
    
     if t<=6*tp
      vis=50;% this is used to save or save plot at times that is 
               % multiple of 40 10 
    else
         vis=100;% this is used to save or save plot at times that is 
               % multiple of 40 50
    end


%%% optical parameters of first material
if ~(strcmp(material,'Mo'))==1
[n_1,k_1, alpha1]=material_optical_parameters(material,Te1,TL1, Ae, BL, N1, wavelength);
else 
    n_1=n_1;
    k_1=k_1;
 alpha1=(4*pi)*k_1/wavelength*ones(size(Te1)); %1/m

end 

%%%% optical parameters of substrate

 alpha2=4*pi*k_2/(wavelength); % 1/m


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



 
      


if  strcmp(material,'Ti')==1 | strcmp(material,'Steel')==1 | strcmp(material,'Pt')==1 | strcmp(material,'W')==1
 
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
           
    elseif strcmp(material,'Mo')==1
        
           ballistic=20e-3*1e-6; %m 
           
         elseif strcmp(material,'Cr')==1 
             
              ballistic=14e-3*1e-6; %m 
              
                   elseif strcmp(material,'Ni')==1 
             
              ballistic=11e-3*1e-6; %m 
    end 
    
end 


enhanced=1/(1/alpha1(1)+1*ballistic);
coefent=1-1*exp(-L1*enhanced);


% Calculate the heat source that attenuates



Icoeff = 2/ (tp * sqrt(pi)) * sqrt(log(2)) * Ep/2 * exp(-4 * log(2) * (t-3*tp).^2 / tp^2)+...
    2/ (tp * sqrt(pi)) * sqrt(log(2)) * Ep/2 * exp(-4 * log(2) * (t-3*tp-t_delay).^2 / tp^2);



Source1=enhanced*exp(-enhanced*(0:dx1:L1)')/coefent*Icoeff*absorpt;


% second film

second_material_matrix=(N1-1)*dx1:dx2:L1+L2;


 %Source2=1*Source1(end)/enhanced*exp(alpha2(1)*(L1-second_material_matrix))*alpha2(1);
if ~(strcmp(material_substrate,'Air')==1)
Source2=1*transmissivity*Icoeff*exp(alpha2(1)*(L1-second_material_matrix))*alpha2(1);
else
    Source2=0*transmissivity*Icoeff*exp(alpha2(1)*(L1-second_material_matrix))*alpha2(1);;
    
end 


        % heat capacity lattice first material
         CL_L= CL_l*ones(size(TL1));
        CL_S= CL_s*ones(size(TL1));
        CL_S(find(TL1>=Tmelt))=CL_l;
        
        
        % 2 Heat conductivity of electron first material
      
        ke1=ke0*Te1./(AB_new.*Te1.^2+TL1);
        
        % 2 Heat conductivity of lattice first material
         kL1=1*0.01*ke1;
       
        % 3 Electron-phono coupling and heat capacity of electron -first
        % material
        coupl_const=ones(N1+1,1);
         CE=ones(N1+1,1);
         
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
coupl_const=coupl_const(:);


     %%%
  
     dTe1_dt=zeros(size(Te1));
     dTL1_dt=zeros(size(TL1));
       dTL2_dt=zeros(size(TL2));
    
      
     
          if t<=1e-15
              dt=0.8*1/2*(dx1^2*min(CE(:))/max(ke1(:)));
             
              
%         
         else 
           dt=0.8*1/2*(dx1^2*min(CE(:))/max(ke1(:)));
            end 
             
   
           Te1b=Te1;
   
           
i=2:size(dTe1_dt,1)-1;% 

A=0; 
  

% %    % ke*theta^2 Te/theta z^2
%    A=A+(ke0*Te1(i)./((AB_new(i).*Te1(i).^2+TL1(i))*dx1^2)).*(Te1(i+1)+Te1(i-1)-2*Te1(i));
% 
%   % theta ke/theta Te*(theta Te/theta z)^2
%     A1=ke0*(AB_new(i).*Te1(i).^2+TL1(i))-2*ke0*AB_new(i).*Te1(i).^2;
%     A1=A1./(AB_new(i).*Te1(i).^2+TL1(i)).^2;
%   A=A+A1.*((Te1(i+1)-Te1(i-1))/(2*dx1)).^2 ;
%   
%    % theta ke/theta TL*(theta Te/theta z)*(theta TL/theta z)
%   
%   A=A+(-ke0*Te1(i)./(AB_new(i).*Te1(i).^2+TL1(i)).^2).*(Te1(i+1)-Te1(i-1))/(2*dx1).*...
%       (TL1(i+1)-TL1(i-1))/(2*dx1);
%   
j1=1;

 
 
dTe1_dt(i,j1)=((ke1(i+1,j1)+ke1(i,j1))./(2*dx1).*(Te1(i+1,j1)-Te1(i,j1))/(dx1)-...
    (ke1(i,j1)+ke1(i-1,j1))./(2*dx1).*(Te1(i,j1)-Te1(i-1,j1))/(dx1)-...
    coupl_const(i,j1).*(Te1(i,j1)-TL1(i,j1))+Source1(1:end-1))./CE(i,j1);

    
%end 
 
%     dTe1_dt(i)=(A-coupl_const(i).*(Te1(i)-TL1(i))+Source1)./CE(i);
           



Te1(i)=Te1(i)+dt*dTe1_dt(i);
% Te1(end)=Te1(end-1);
% 
% Te1b(end)=Te1b(end-1);


  i=2:size(dTL1_dt,1)-1; % 

 
          
dTL1_dt(i)=((kL1(i+1)+kL1(i))./(2*dx1).*(TL1(i+1)-TL1(i))/(dx1)-...
    (kL1(i)+kL1(i-1))./(2*dx1).*(TL1(i)-TL1(i-1))/(dx1)+...
    coupl_const(i).*(Te1b(i)-TL1(i)))./CL_S(i);
       


TL1(i)=TL1(i)+dt*dTL1_dt(i);



i=2:size(dTL2_dt,1)-1 ;


dTL2_dt(i)=(kL2*(TL2(i+1)+1*TL2(i-1)-2*TL2(i))/dx2^2+1*Source2(1:end)')/CL_s2;

TL2(i)=TL2(i)+dt*dTL2_dt(i);




%       
Te1=Te1(2:end-1,1);
TL1=TL1(2:end-1,1);
TL2=TL2(2:end-1,1);



   
        ke1a=ke0*Te1(end)./(AB_new(1,1).*Te1(end).^2+TL1(end));  
        % 2 Heat conductivity of lattice first material
         kL1a=1*0.01*ke1a;


numerator=kL1a*TL1(end)/dx1+kL2*TL2(1)/dx2;
denominator=kL1a/dx1+kL2/dx2;
TL12=numerator/denominator;   %%% on the interface


Te1=[Te1;Te1(end)];
TL1=[TL1;TL12];
TL2=[TL12;TL2];

    if mod(AA,vis)==0 | AA==1
Ta=[Ta,t];
TE=[TE,Te1(1)];
TL=[TL,TL1(1)];
REFL=[REFL,reflectivity];
ABSO=[ABSO,absorpt];
TRANS=[TRANS,transmissivity];



 Te1_matrix=[Te1_matrix;Te1'];
 TL1_matrix=[TL1_matrix;TL1'];
  TL2_matrix=[TL2_matrix;TL2'];


    end 

t=t+dt;


end 

  Te_sol_1=Te1_matrix;
  TL_sol_1=TL1_matrix; %[0 L1]
  TL_sol_2=TL2_matrix; %[L1 L2]
  


   A=[TL_sol_1(:,1:end)';TL_sol_2(:,2:end-1)'];
   
 
 [r,c]=find(TL_sol_1>=Tmelt );
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
t=Ta(:);
clear TA;
Icoeff = 1/2*exp(-4 * log(2) * (t-3*tp).^2 / tp^2)+...
   1/2 * exp(-4 * log(2) * (t-3*tp-t_delay).^2 / tp^2);


% 


                        DATA.material=material;
                        DATA.substrate=material_substrate;
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
                        DATA.reflectivity=REFL;
                        DATA.absorpt=ABSO;
                        DATA.transmissivity=TRANS;
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
function [n_1,k_1, alpha]=material_optical_parameters(material,Te1,TL1,Ae,BL,N1, wavelength)

global n_1 k_1 
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Material optical parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(material,'Ni')==1 | strcmp(material,'Au')==1|...
   strcmp(material,'Steel')==1| strcmp(material,'Ti')==1 ...
        | strcmp(material,'Cu')==1 | strcmp(material,'Ag')==1 | strcmp(material,'Cr')==1 | strcmp(material,'Al')==1 | strcmp(material,'W')==1 | strcmp(material,'Pt')==1
 
 
    
    
     
AB_new=Ae/BL*ones(N1+1,1); %1/K 


if strcmp(material,'Ni')==1 | strcmp(material,'Au')==1| strcmp(material,'Ti')==1| strcmp(material,'Cu')==1 ...
        | strcmp(material,'Ag')==1 | strcmp(material,'Cr')==1 | strcmp(material,'Al')==1 | strcmp(material,'W')==1 | strcmp(material,'Pt')==1


[reflect1,alpha,e1,e2]= LD_Tin_new(wavelength,material,'LD',Te1,TL1,AB_new, BL);  


alpha=alpha/1e-6; %1/m

%reflectivity=reflect1(2); %reflectivity of material;
penetration_depth=1/alpha(1);% in m  


e_1=e1(1);
e_2=e2(1);

 n_1=sqrt((e_1+sqrt(e_1^2+e_2^2))/2);
  k_1=sqrt((-e_1+sqrt(e_1^2+e_2^2))/2);
 


elseif  strcmp(material,'Steel')==1
    
 
    
    if wavelength==1.026*1e-6
    
      [pAaa,pBaa]=reflectivity_absorption_Fe(1);

    reflectivity=polyval(pAaa,Te1);
     penetration_depth=1e-3*polyval(pBaa,Te1)*1e-6;%m

     alpha=1./penetration_depth;
 

     k_1=alpha(1)*wavelength/(4*pi);

    trion_A=reflectivity(1)-1;
    trion_B=2*reflectivity(1)+2;
    trion_C=reflectivity(1)-1-k_1^2+reflectivity(1)*k_1^2;
    n_1=(-trion_B+sqrt(trion_B^2-4*trion_A*trion_C))/(2*trion_A);

    else 
        n_1=n_1;
        k_1=k_1;
       
          alpha=(4*pi)*k_1/wavelength*ones(size(Te1)); %1/m
         
    end 
%      
%     reflectivity=reflectivity(1); %reflectivity of material;



end   
 
end
end 
