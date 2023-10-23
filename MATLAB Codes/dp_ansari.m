clc;
clear;
close all;
a=input('input deviation angle from horizontal(degrees):');
a1=(a*pi)/180;
MUg1=input('input gas viscozity(cp):');
MUg=MUg1*0.001;%pa.s
MUo1=input('input oil viscozity(cp):');
MUo=MUo1*0.001;%pa.s
MUw1=input('input water viscozity(cp):');
MUw=MUw1*0.001;%pa.s
d=input('input diameter(in):');
d2=d/12;
d1=d*0.0254;%m
WOR=input('input water oil ratio:');
qo_prim=input('input oil rate(STB/day):');
qw_prim=input('input water rate (STB/day):');
qg=input('input gas rate(ft^3/s):');
ROo1=input('input oil density(lbm/ft^3):');
ROo=ROo1*16.018;%kg/m^3
ROg1=input('input gas density(lbm/ft^3):');
ROg=ROg1*16.018;%kg/m^3
ROw1=input('input water density(lbm/ft^3):');
ROw=ROw1*16.018;%kg/m^3
ZG_o1=input('input oil surface tension(dynes/cm):');
ZG_o=ZG_o1*0.001;%N/m
ZG_w1=input('input water surface tension(dynes/cm):');
ZG_w=ZG_w1*0.001;%N/m
Bo=input('input oil FVF(bbl/STB):');
Bw=input('input water FVF(bbl/STB):');
g=9.81;
qo=qo_prim*Bo;%bbl/day
qw=qw_prim*Bw;%bbl/day
ql=(qo+qw)*(5.617/(24*3600));%ft^3/s
S=(pi/4)*d2^2;
Vsg1=qg/S;
Vsg=Vsg1*0.3048;%m/s
Vsl1=ql/S;
Vsl=Vsl1*0.3048;%m/s
fo=1/(1+WOR);
fw=1-fo;
ROl=ROo*fo+ROw*fw;%kg/m^3
MUl=MUo*fo+MUw*fw;
ZG=ZG_o*fo+ZG_w*fw;
Vm=Vsl+Vsg;
Vb=1.53*(((g*(ROl-ROg)*ZG)/(ROl)^2)^0.25); %manzor az Vb Vbinahayat mibashad.
VbT=(0.35*sin(a1)+0.54*cos(a1))*((g*d1*(ROl-ROg)/ROl)^0.5);%VbT haman Vbinahayat teylor bubble mibashad.
VTB=1.2*Vm+VbT; 
fgls=Vsg/(1.208*Vm+1.41*((g*(ROl-ROg)*ZG/(ROl)^2)^0.25)*sqrt(sin(a1)));
Vgls=1.08*Vm+1.4*((g*ZG*(ROl-ROg)/(ROl)^2)^0.25)*sqrt(sin(a1));
VLLS=(Vm-Vgls*fgls)/(1-fgls);
syms Fltb
F=9.916*sqrt(g*d1)*(1-sqrt(1-Fltb))^0.5-(((VTB-VLLS)*(1-fgls))/Fltb)+VTB;
dif_F=diff(F);
F=inline(F);
dif_F=inline(dif_F);
error=1;
fltb0=0.005;
while abs(error)>1e-5
    fltb=fltb0-feval(F,fltb0)/feval(dif_F,fltb0);
    error=fltb-fltb0;
    fltb0=fltb;
end

VLTB=9.916*(sqrt(g*d1)*((1-sqrt(1-fltb))^0.5));
%A=Lls/Lsu
A=(Vsl+(VTB-VLLS)*(1-fgls)-VTB*fltb)/(VTB*((1-fgls)-fltb));
ROls=ROl*(1-fgls)+ROg*fgls;
MUls=MUl*(1-fgls)+MUg*fgls;
NRels=(ROls*Vm*d1)/MUls;
if NRels<2100
    fls=64/(NRels);
elseif NRels>50000
    fls=0.184*(NRels)^(-0.2);
else
    fls=0.32*(NRels)^(-0.25);
end
dpdz1=A*ROls*g*sin(a1)+A*fls*ROls*(Vm^2)/2*d1;
dpdz=dpdz1*0.0063664;
disp(['meghdar oft feshar barabar ast ba(psf/ft):' num2str(dpdz)]);