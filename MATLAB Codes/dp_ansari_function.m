function [dpdz]=dp_ansari_function(a1,d1,MUg,ROg,Vsg,Vsl,ROl,MUl,ZG,g,gc)
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
dpdz=(A*ROls*(g/gc)*sin(a1)+A*fls*ROls*(Vm^2)/(2*d1*gc))/148.73;
%disp(['meghdar oft feshar barabar ast ba(psi/ft):' num2str(dpdz)]);
end