function [p1,p2]=reflectivity_absorption_Fe(q)
warning off;
% this script helps to fit the curves for the capacity and 
% coupling constant data for Ti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fe 
% 
% reflectivity
%
A=load('reflectivity_iron.txt');

x1=(A(:,1));y1=A(:,2);


%   figure;plot(x1(1:end)/1e4,y1(1:end)/1e5,'k');
% dfgdf

p1=polyfit(x1,y1,7);
y1a=polyval(p1,x1);
%   figure;plot(x1,y1,x1,y1a)

% absorption
%
A=load('skin_depth_iron.txt');

x1=(A(:,1));y1=A(:,2);



p2=polyfit(x1,y1,7);
y1a=polyval(p2,x1);
%  figure;plot(x1,y1,x1,y1a)
% 
% bhjgh

end

