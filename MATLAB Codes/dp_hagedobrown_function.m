function[dpdztotal]=dp_hagedobrown_function(Vsl,MUl,Pavg,Vsg,MUg,ROl,ZG,d1,g,gc,ROg)
%game badi mohasebate adade bedone bood mibashad.
Nlv=Vsl*(ROl/(g*ZG))^0.25;
Nl=MUl*(1/(ROl*(ZG)^3))^0.25;
Nd=d1*((g*ROl)/ZG)^0.5;
Ngv=Vsg*(ROl/(g*ZG))^0.25;

%hala mikhahim parametr haye say,CNl,fl/say ra be dast biyavarim.
B=(Ngv*(Nl)^0.38)/(Nd)^2.14;
W=(1.0886-69.9473*B+2334.349*(B)^2-12896.638*(B)^3)/(1-53.4401*B+1517.936*(B)^2-8419.8115*(B)^3);

CNl=(0.0019+0.032*Nl-0.6642*(Nl)^2+4.995*(Nl)^3)/(1-10.01*Nl+33.869*(Nl)^2+277.28*(Nl)^3);

H=(Nlv/(Ngv)^0.575)*((Pavg/14.7)^0.1)*(CNl/Nd);
flW=((0.0047+1123.32*H+429489.64*(H)^2)/(1+1097.1566*H+722153.95*(H)^2))^0.5;

%hala bayad liquid holdup ra mohasebe konim.
fl=flW*W;
%hala MUs va ROs ra mohasebe mikonim.
ROs=(ROl*fl)+(ROg*(1-fl));
MUs=(MUl^fl)*(MUg^(1-fl)); %vahed MUs=lb/ft.sec mibashad.

%hala nobate mohasebeye reynoldz va firiction factor mibashad.
Vm=Vsg+Vsl;
Re=(ROs*Vm*d1)/MUs;
f=0.32*(Re)^(-0.25); %Blasius

%dar nahayat be mohasebate ofte feshar mipardazim.
%gc=32.174;
dpdzf=(f*ROs*(Vm)^2)/2*gc*d1; %vahed=psf/ft.
dpdzH=(g/gc)*ROs; %vahed=psf/ft.
dpdztotal=(dpdzf+dpdzH)/148.73; %vahed=psi/ft.
%dpdztotal=dpdztotal/148.73; %vahed=psi/ft.
%disp(['(dp/dz)total ba estefade az raveshe hagedobrown barabar ast ba(psi/ft)=' num2str(dpdztotal)]);
end