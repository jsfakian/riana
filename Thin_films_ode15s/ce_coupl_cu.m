function [thresh1,p1,p2,p3,p4,Xaa]=ce_coupl_ti(q)
warning off;
% this script helps to fit the curves for the capacity and 
% coupling constant data for Ti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cu 
% 
% capacity
%
A=load('Ce_Cu.dat');

x=(1e4*A(:,1))';y=(1e5*A(:,2))';
Xaa=max(x); %maximujm temperature used in the simulations fitting
 %figure;plot(x(1:10:end)/1e4,y(1:10:end)/1e5,'k.');
x1=x(find(x<=1e4*0.205));
y1=y(find(x<=1e4*0.205));

p1=polyfit(x1,y1,1);
y1=polyval(p1,x1);

x2=[x(find(x>=1e4*0.205))];
y2=[ y(find(x>=1e4*0.205))];

p2=polyfit(x2,y2,11);
y2=polyval(p2,x2);
thresh1=1e4*0.205;

 % hold on;plot([x1/1e4 x2/1e4],[y1/1e5 y2/1e5],'k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coupling constant
%%%%%%%%%%%%%%%%%%%%%%%


A=load('G_Cu.dat');

x=(1e4*A(:,1))';y=(1e17*A(:,2))';
%figure;plot(x(1:10:end)/1e4,y(1:10:end)/1e17,'k.');
x1=x(find(x<=1e4*0.205));
y1=y(find(x<=1e4*0.205));

p3=polyfit(x1,y1,5);
y1=polyval(p3,x1);

x2=[x(find(x>=1e4*0.205))];
y2=[ y(find(x>=1e4*0.205))];



p4=polyfit(x2,y2,11);
y2=polyval(p4,x2);
%hold on;plot([x1/1e4 x2/1e4],[y1/1e17 y2/1e17],'k')

