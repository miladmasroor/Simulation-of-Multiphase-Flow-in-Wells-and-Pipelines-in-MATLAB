function[dp_dz]=stratified_flow_model(ROl,ROg,Vsl,Vsg,d,MUl,MUg,a,HL)
% ROl & ROg=kg/m^3
% Vsl & Vsg=m/s
% d=in
% MUl & MUg=pa.s
% ZG=N/m
a1=(a*pi)/180;
d1=d/12;%ft
d2=d*0.0254;%m
g=9.81;
A=(pi/4)*(d2^2);%m^2
teta=2*acos(1-HL/d2);
EL=(teta-sin(teta))/(2*pi);
Vs=Vsl+Vsg;%Vs=Vm=mixture velocity
AL=HL*pi*d2/4;
AG=((pi/4)*(d2^2))-AL;
SL=HL*pi;
SG=(pi*d)-SL;
Si=sqrt(4*(d2-(HL)^2));
DL=4*AL/SL;
DG=4*AG/(SG+Si);
qL=Vsl*A;
qg=Vsg*A;
VL=qL/AL;
Vg=qg/AG;
 NReL=ROl*VL*DL/MUl;
if NReL<2100
    fWL=64/(NReL);
elseif NReL>50000
    fWL=0.184*(NReL)^(-0.2);
else
    fWL=0.32*(NReL)^(-0.25);
end
  tWL=fWL*ROl*(VL^2)/2;
  NReg=ROg*Vg*DG/MUg;
if NReg<2100
    fWg=64/(NReg);
elseif NReg>50000
    fWg=0.184*(NReg)^(-0.2);
else
    fWg=0.32*(NReg)^(-0.25);
  end
tWg=fWg*ROg*(Vg^2)/2;
dp_dz=(tWL*SL+tWg*SG)+((AL/A)*ROl+(AG/A)*ROg)*g*sin(a1);
end