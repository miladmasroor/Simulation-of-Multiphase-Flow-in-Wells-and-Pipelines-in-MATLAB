clc;
clear;
disp('shoma darid az raveshe homogenus estefade mikonid.');
disp('etelaat dade shode dar sorat soal ra vared konid.');
Vsl=input('Vsl ra vared kon(ft/sec):');
MUl=input('MUl ra vared kon(cp):');
MUl1=MUl/1488.16; %viscozity ra az cp be lb/ft.sec tabdil vahed kardim.
p=input('p ra vared kon(psig):');
p1=p+14.7; %feshar ra az psig be psia tabdil vahed mikonim.
Vsg=input('Vsg ra vared kon(ft/sec):');
MUg=input('MUg ra vared kon(cp):');
MUg1=MUg/1488.16; %viscozity ra az cp be lb/ft.sec tabdil vahed kardim.
d=input('d ra vared kon(in):');
d1=d/12; %ghotr ra az in be ft tabdil vahed kardim.
ROl=input('ROl ra vared kon(lbm/ft^3):');
ZG=input('ZG ra vared kon(dynes/cm):');
ZG1=ZG/453.632; %zigma ra az dynes/cm be lbm/sec^3 tabdil vahed kardim.
ROg=input('ROg ra vared kon(lbm/ft^3):');
g=32.174;
x=(Vsg*ROg)/(Vsg*ROg+Vsl*ROl);
%hala mirim soraghe mohasebeye oft feshar.
gc=32.174;
MUs=(MUl1*(1-x))+MUg1*x; 
Vm=Vsg+Vsl;
fg=Vsg/Vm;
ros=((1-x)/ROl)+x/ROg;
ROs=1/ros;
Gm=(Vsg*ROg)+(Vsl*ROl);
Rem=(Gm*d1)/MUs;
fm=0.32*(Rem)^(-0.25);
dpdzf=(fm*((Vm)^2)*ROs)/(2*gc*d1);%vahed=psf/ft.
dpdzH=ROs;%vahed=psf/ft.
dpdztotal=dpdzf+dpdzH; %vahed=psf/ft.
dpdztotal1=dpdztotal/148.73; %vahed=psi/ft.
disp(['(dp/dz)total ba estefade az raveshe homogeneus barabar ast ba(psi/ft)=' num2str(dpdztotal1)]);












