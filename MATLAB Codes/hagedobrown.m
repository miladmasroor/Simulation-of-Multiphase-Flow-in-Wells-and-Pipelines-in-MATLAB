clc;
clear;
disp('shoma darid az raveshe hagedobrown estefade mikonid.');
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
%game badi mohasebate adade bedone bood mibashad.
Nlv=Vsl*(ROl/(g*ZG1))^0.25;
Nl=MUl1*(1/(ROl*(ZG1)^3))^0.25;
Nd=d1*((g*ROl)/ZG1)^0.5;
Ngv=Vsg*(ROl/(g*ZG1))^0.25;

%hala mikhahim parametr haye say,CNl,fl/say ra be dast biyavarim.
B=(Ngv*(Nl)^0.38)/(Nd)^2.14;
W=(1.0886-69.9473*B+2334.349*(B)^2-12896.638*(B)^3)/(1-53.4401*B+1517.936*(B)^2-8419.8115*(B)^3);

CNl=(0.0019+0.032*Nl-0.6642*(Nl)^2+4.995*(Nl)^3)/(1-10.01*Nl+33.869*(Nl)^2+277.28*(Nl)^3);

H=(Nlv/(Ngv)^0.575)*((p1/14.7)^0.1)*(CNl/Nd);
flW=((0.0047+1123.32*H+429489.64*(H)^2)/(1+1097.1566*H+722153.95*(H)^2))^0.5;

%hala bayad liquid holdup ra mohasebe konim.
fl=flW*W;
%hala MUs va ROs ra mohasebe mikonim.
ROs=(ROl*fl)+(ROg*(1-fl));
MUs=(MUl1^fl)*(MUg1^(1-fl)); %vahed MUs=lb/ft.sec mibashad.

%hala nobate mohasebeye reynoldz va firiction factor mibashad.
Vm=Vsg+Vsl;
Re=(ROs*Vm*d1)/MUs;
f=0.32*(Re)^(-0.25); %Blasius

%dar nahayat be mohasebate ofte feshar mipardazim.
gc=32.174;
dpdzf=(f*ROs*(Vm)^2)/2*gc*d1; %vahed=psf/ft.
dpdzH=(g/gc)*ROs; %vahed=psf/ft.
dpdztotal=dpdzf+dpdzH; %vahed=psf/ft.
dpdztotal1=dpdztotal/148.73; %vahed=psi/ft.
disp(['(dp/dz)total ba estefade az raveshe hagedobrown barabar ast ba(psi/ft)=' num2str(dpdztotal1)]);






























