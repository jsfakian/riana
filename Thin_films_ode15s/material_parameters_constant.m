function [ rho_s,rho_l,CL_s,CL_l,Tmelt,Ae,BL,ke0,T_cr] = material_parameters_constant(material)
% The function material_parameters yields the thermophysical properties of
% the material. 
%


if strcmp(material,'Ni')==1|strcmp(material,'Ni_xuv')==1
 rho_l=7500; % in Kgr/m^3 (in liquid phase)
 rho_s=8900; % in Kgr/m^3 (in solid phase)
 CL_s=457; % in J/(Kgr*K)
 CL_l=734.16; % in J/(Kgr*K) (molten phase)
 Tmelt=1728  ;% K
 Ae=0.59*1e7; %s-1K=2; 
 BL=1.4e11 ; %s-1K-1;
 ke0=90; %Jm-1s-1K-1; 90/(8900*734)
 T_cr=9460; %K
 
elseif strcmp(material,'Au')==1
    
 rho_l=17310; % in Kgr/m^3 (in liquid phase)
 rho_s=19300; % in Kgr/m^3 (in solid phase)
 CL_s=129; % in J/(Kgr*K)
 CL_l=129; % in J/(Kgr*K) (molten phase) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Tmelt=1337  ;% K
 Ae=1.18*1e7; %s-1K=2; 
 BL=1.25e11 ; %s-1K-1;
 ke0=318; %Jm-1s-1K-1; 318(19300*129)
 T_cr=6250; %K
 
 
 elseif strcmp(material,'Cu')==1
    
 rho_l=8020; % in Kgr/m^3 (in liquid phase)
 rho_s=8960; % in Kgr/m^3 (in solid phase)
 CL_s=380; % in J/(Kgr*K)
 CL_l=380; % in J/(Kgr*K) (molten phase) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Tmelt=1357  ;% K
 Ae=1.28*1e7; %s-1K=2; 
 BL=1.23e11 ; %s-1K-1;
 ke0=401; %Jm-1s-1K-1;
 T_cr=6250; %K   %%%%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFY
 Tcr=7700;   % set boiling point
 
  elseif strcmp(material,'Cr')==1
    
 rho_l=6300; % in Kgr/m^3 (in liquid phase)
 rho_s=7190; % in Kgr/m^3 (in solid phase)
  CL_s=449.23; % in J/(Kgr*K)
 CL_l=449.23; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=7.9*1e7; %s-1K=2;    CHEN
 BL=13.4e11 ; %s-1K-1;    CHEN
 ke0=93.9; %Jm-1s-1K-1;
 T_cr=2944; %K  %%%%%%%%%%%%%%%%SPECIFY--- this i sthe boiling point
 Tmelt=2180  ;% K
 
 
 
 
 elseif strcmp(material,'Steel')==1
    
 rho_l=7000; % in Kgr/m^3 (in liquid phase) 
 rho_s=7750; % in Kgr/m^3 (in solid phase)
 CL_s=475; % in J/(Kgr*K)
 CL_l=748; % in J/(Kgr*K) (molten phase) 
 Ae=0.98*1e7; %s-1K=2; 
 BL=2.8e11 ; %s-1K-1;
 ke0=46.6; %Jm-1s-1K-1;
 T_cr=8500; %K
 Tmelt=1811  ;% K
 
 
 elseif strcmp(material,'Ti')==1
    
 rho_l=4110; % in Kgr/m^3 (in liquid phase) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 rho_s=4506; % in Kgr/m^3 (in solid phase)
 CL_s=522; % in J/(Kgr*K)
 CL_l=522; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=1*1e7; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%approximately?
 BL=1.5e11 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% approximatekly?
 ke0=21.9; %Jm-1s-1K-1;  21.9/(2506*522)
 T_cr=15500; %K
 Tmelt=1941  ;% K
 Tboil=3287+273;
%  T_cr=Tboil;
 
    
  elseif strcmp(material,'Pt')==1
    
 rho_l=19770; % in Kgr/m^3 (in liquid phase) 
 rho_s=21452; % in Kgr/m^3 (in solid phase)
 CL_s=0.13*1000; % in J/(Kgr*K)
 CL_l=0.13*1000; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=1*1e7; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%approximately?
 BL=1.5e11 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% approximatekly?
 ke0=72; %Jm-1s-1K-1;  21.9/(2506*522)
 T_cr=13000; %K
 Tmelt=1772+273  ;% K
 Tboil=3827+273;
%  T_cr=Tboil;
 
 
  elseif strcmp(material,'Mo')==1

 rho_s=10223; % in Kgr/m^3 (in liquid phase) 
 rho_l=9330; % in Kgr/m^3 (in solid phase)
 CL_s=0.25*1000; % in J/(Kgr*K)
 CL_l=0.25*1000; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=1e7; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%approximately?
 BL=1.5e11 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% approximatekly?
 ke0=138; %Jm-1s-1K-1;
 T_cr=9450; %K
 Tmelt=2896  ;% K
  
        
 
 
 elseif strcmp(material,'W')==1
    
 rho_l=17600; % in Kgr/m^3 (in liquid phase) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 rho_s=19254; % in Kgr/m^3 (in solid phase)
 CL_s=0.13*1000; % in J/(Kgr*K)
 CL_l=0.13*1000; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=1*1e7; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%approximately?
 BL=1.5e11 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% approximatekly?
 ke0=21.9; %Jm-1s-1K-1;
 T_cr=15500; %K
 Tmelt=1941  ;% K
 
  elseif strcmp(material,'Ag')==1
    
 rho_l=9320; % in Kgr/m^3 (in liquid phase)
 rho_s=10490; % in Kgr/m^3 (in solid phase)
 CL_s=238; % in J/(Kgr*K)
 CL_l=238; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=0.932*1e7; %s-1K=2;    
 BL=1.02e11 ; %s-1K-1;    
 ke0=428; %Jm-1s-1K-1;
 T_cr=7700; %K 
 Tmelt=1234  ;% K
  Tboil=2435;% https://periodictable.com/Elements/047/data.html
%  T_cr=Tboil;
 
   elseif strcmp(material,'Al')==1

 rho_l=2350; % in Kgr/m^3 (in liquid phase)
 rho_s=2700; % in Kgr/m^3 (in solid phase)
 CL_s=921; % in J/(Kgr*K)
 CL_l=921; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY%%%%%%%%%%%
 Ae=0.376*1e7; %s-1K=2;    
 BL=3.9e11 ; %s-1K-1;    
 ke0=235; %Jm-1s-1K-1; 235/(2700*921)
 T_cr=6700; %K  The Critical Temperature of Aluminum 
 Tmelt=933.47  ;% K

 
 
 
    elseif strcmp(material,'GaN')==1
    
 rho_l=6100; % in Kgr/m^3 (in liquid phase)  %%%%% specifiy
 rho_s=6100; % in Kgr/m^3 (in solid phase)
 CL_s=726; % in J/(Kgr*K)
 CL_l=726; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 Ae=0; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT applied
 BL=0 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT applied
 ke0=1.3*100; %Jm-1s-1K-1;  %%%KL 130/(726*6100)
 T_cr=1600+273; %K  %%%%%%%%%%%%%% spcify  
 Tmelt=1600+273  ;% K
    
 
 
 
 
 
 
 
    elseif strcmp(material,'SiO2')==1
    
 rho_l=2050; % in Kgr/m^3 (in liquid phase)
 rho_s=2203; % in Kgr/m^3 (in solid phase)
 CL_s=726; % in J/(Kgr*K)
 CL_l=726; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 Ae=0; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT applied
 BL=0 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT applied
 ke0=1.3; %Jm-1s-1K-1;  %%%KL
 T_cr=7000; %K  
 Tmelt=1738  ;% K
    
    
     elseif strcmp(material,'Si')==1
    
 rho_l=2570; % in Kgr/m^3 (in liquid phase)
 rho_s=2329; % in Kgr/m^3 (in solid phase)
 CL_s=706; % in J/(Kgr*K)
 CL_l=706; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 Ae=0; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT applied
 BL=0 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT applied
 ke0=149; %Jm-1s-1K-1;  %%%KL
 T_cr=5159; %K  
 Tmelt=1687  ;% K
    
 
     elseif strcmp(material,'Soda')==1
    
 rho_l=2500; % in Kgr/m^3 (in liquid phase)%%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 rho_s=2500; % in Kgr/m^3 (in solid phase)
 CL_s=870; % in J/(Kgr*K)
 CL_l=870; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 Ae=0; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT applied
 BL=0 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT applied
 ke0=1.06; %Jm-1s-1K-1;  %%%KL
 T_cr=5159; %K %%%%%%%%%%%%%%% not applied  
 Tmelt=1273  ;% K
 
      elseif strcmp(material,'Air')==1
    
 rho_l=1.1455	; % in Kgr/m^3 (in liquid phase)%%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 rho_s=1.1455	; % in Kgr/m^3 (in solid phase)
 CL_s=1000; % in J/(Kgr*K)
 CL_l=1000; % in J/(Kgr*K) %%%%%%%%%%%%SPECIFY %%%%%%%%%%%
 Ae=0; %s-1K=2;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT applied
 BL=0 ; %s-1K-1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT applied
 ke0=0.0257; %Jm-1s-1K-1;  %%%KL
 T_cr=5159; %K %%%%%%%%%%%%%%% not applied  
 Tmelt=1273  ;% K %%%%%%%%%%%%%% not applied  
 
 
end







