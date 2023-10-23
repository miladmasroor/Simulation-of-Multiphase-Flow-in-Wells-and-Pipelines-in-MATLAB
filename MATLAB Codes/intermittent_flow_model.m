function [dp_dz]=intermittent_flow_model(ROl,ROg,Vsl,Vsg,d,MUl,MUg,ZG,a,HL)
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
 teta=2*acos((1-2*HL/d2)^2);
 Ef=(teta-sin(teta))/(2*pi)
if d<1.5
    Ls=30*d2;%slug length=Ls
else 
    Ls=exp((-26.6)+28.5*((log(d2)+3.67)^0.1))%slug length=Ls
end
  Vs=Vsl+Vsg;%Vs=Vm=mixture velocity
  Es=1/(((Vs/8.66)^1.39)+1);%Es=liquid holdup in slug Body
  Vb=1.2*Vs+(1.53*(ZG*g*(ROl-ROg)/(ROl^2))^0.25)*(Es^0.1)*sin(a1);%Vb=The velocity of dispersed
  ROs=Es*ROl+(1-Es)*ROg;
  MUs=Es*MUl+(1-Es)*MUg;
  NRes=ROs*Vs*d2/MUs;
  if NRes<2100
    fs=64/(NRes);
elseif NRes>50000
    fs=0.184*(NRes)^(-0.2);
else
    fs=0.32*(NRes)^(-0.25);
  end
ts=fs*ROs*(Vs^2)/2;
if NRes<2100 %laminar
    C=2;
else
    C=1.2;
end
Vt=C*Vs+0.35*(sqrt(g*d2))*sin(a1)+0.54*(sqrt(g*d2))*cos(a1);
VL=(Vs-((1-Es)*Vb))/Es;
EL=(Vt*Es+Vb*(1-Es)-Vsg)/Vt;
Vf=Vt-((Vt-VL)*Es/Ef);
Vg=(Vs-Vf*Ef)/(1-Ef);
Lu=Ls*(((VL*Es)-(Vf*Ef))/(Vsl-(Vf*Ef)))
AL=HL*pi*d2/4;
AG=((pi/4)*(d2^2))-AL;
SL=HL*pi;
SG=(pi*d2)-SL;
Si=sqrt((((d2)^2)/4)-(((d2/2)-HL)^2));
DL=4*AL/SL;
DG=4*AG/(SG+Si);
NRef=ROl*Vf*DL/MUl;
if NRef<2100
    ff=64/(NRef);
elseif NRef>50000
    ff=0.184*(NRef)^(-0.2);
else
    ff=0.32*(NRef)^(-0.25);
end
  tf=ff*ROl*(Vf^2)/2;
NReg=ROg*Vg*DG/MUg;
if NReg<2100
    fg=64/(NReg);
elseif NReg>50000
    fg=0.184*(NReg)^(-0.2);
else
    fg=0.32*(NReg)^(-0.25);
  end
tg=fg*ROg*(Vg^2)/2;
ROu=EL*ROl+(1-EL)*ROg;
Lf=(EL*Lu-Es*Ls)/Ef
dp_dz=ROu*g*sin(a1)+(1/Lu)*((ts*pi*d2*Ls/A)+((tf*SL+tg*SG)*Lf)/A);
end

