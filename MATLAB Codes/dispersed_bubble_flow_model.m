function[dp_dz]=dispersed_bubble_flow_model(ROl,ROg,Vsl,Vsg,d,MUl,MUg,a)
% ROl & ROg=kg/m^3
% Vsl & Vsg=m/s
% d=in
% MUl & MUg=pa.s
% ZG=N/m
a1=(a*pi)/180;
d1=d/12;%ft
d2=d*0.0254;%m
g=9.81;
A=(pi/4)*(d2^2);
Vm=Vsl+Vsg;
EL=Vsl/Vm;
ROm=ROl*EL+ROg*(1-EL);
MUm=MUl*EL+MUg*(1-EL);
NRem=ROm*Vm*d2/MUm;
  if NRem<2100
    fm=64/(NRem);
elseif NRem>50000
    fm=0.184*(NRem)^(-0.2);
else
    fm=0.32*(NRem)^(-0.25);
  end
  dp_dz=(2*fm*ROm*(Vm^2)/d2)+ROm*g*sin(a1);
end