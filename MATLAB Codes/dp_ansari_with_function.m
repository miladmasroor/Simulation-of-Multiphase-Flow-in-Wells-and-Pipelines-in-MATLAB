  clc;
clear;
close all;
GLR=input('input gas liquid ratio(scf/STB):');
Pbh1=input('input bottom hole pressure(psig):');
Pbh=Pbh1+14.7;%psia
a=input('input deviation angle from horizontal(degrees):');
a1=(a*pi)/180;
d=input('input diameter(in):');
d1=d/12;%ft
WOR=input('input water oil ratio:');
qo_prim=input('input oil rate(STB/day):');
qw_prim=input('input water rate (STB/day):');
%qg=input('input gas rate(ft^3/s):');
ZG_o1=input('input oil surface tension(dynes/cm):');
ZG_o=ZG_o1/453.632;%lb/s^3
ZG_w1=input('input water surface tension(dynes/cm):');
ZG_w=ZG_w1/453.632;%lb/s^3
API=input('input API:');
gamma_g=input('input gas specific gravity: ');
L=input('input length of well(ft):');
N=input('input Number of segments:');
dL=L/N;
Tsc=520;%R
Psc=14.7;%psia
R=10.732;
ROwsc=67;%lbm/ft^3 
Pt(1)=Pbh;
for i=1:N;
%i=1;
   T1=120+0.017*(dL *(i-1));
   T2=120+0.017*(dL*i);
   Tavg=(T1+T2)/2+460;%R
   Pavg=Pbh;%psia
   
   %compressibility factor
   if gamma_g<0.75 
    Tc=168+325*gamma_g-12.5*(gamma_g^2);
    Pc=667+15*gamma_g-37.5*(gamma_g^2);
   else 
    Tc=187+330*gamma_g-71.5*(gamma_g^2);
    Pc=706-51.7*gamma_g-11.1*(gamma_g^2);
   end
    Pr=Pavg/Pc;
    Tr=Tavg/Tc;
    z=zfactor(Tr,Pr);
    %Bg
    Bg=(Psc/Tsc)*(z*Tavg/Pavg);
    %ROg
    ROg=28.97*Pavg*gamma_g/(z*R*Tavg);
    %MUg
    ROg1=0.0433*gamma_g*Pavg/(z*Tavg);%mg/cm^3
   Mg=gamma_g*28.97;
   MUg1=VisLeeGonzalez(Mg,Tavg,ROg1);
   MUg=MUg1/1488.16;%lbm/ft.s
   %Solution Gas_Oil Ratio
   gamma_o=141.5/(131.5+API);
   Tavg1=Tavg-460;%F
   Rs=gamma_g*(((0.055*Pavg+1.4)*(10^(0.0125*API))/(10^(0.00091*Tavg1)))^1.205);%scf/STB
   %oil formation volume factor(Bo)
   F=Rs*((gamma_g/gamma_o)^0.5)+1.25*Tavg1;
   Bo=0.972+0.000147*(F^1.175);%bbl/STB
   %ROo
   ROo=((62.4*gamma_o)+(0.0136*gamma_g*Rs))/Bo;%lbm/ft^3
   %MUo
  AA=((Rs/gamma_g)^0.83)*(10^(0.00091*Tavg1-0.0125*API));
  PB=18.2*(AA-1.4);
  if Pavg<=PB
     ZZ=3.0324-0.02023*API;
     yy=10^ZZ;
     XX=yy*(Tavg1^(-1.163));
     MUoD=(10^XX)-1;%dead oil viscosity
     AAA=10.715*((Rs+100)^(-0.515)) ;
     BBB=5.44*((Rs+150)^(-0.338));
     MUo1=AAA*(MUoD^BBB);
  else
     ZZ=3.0324-0.02023*API;
     yy=10^ZZ;
     XX=yy*(Tavg1^(-1.163));
     MUoD=(10^XX)-1;%dead oil viscosity
     AAA=10.715*((Rs+100)^(-0.515)) ;
     BBB=5.44*((Rs+150)^(-0.338));
     MUoB=AAA*(MUoD^BBB);
     m=2.6*(Pavg^1.187)*exp(-11.513-(8.98*10^(-5))*Pavg) ;
     MUo1=MUoB*((Pavg/PB)^m);
  end
  MUo=MUo1/1488.16;%lbm/ft.s
  %MUw
  MUw1=exp(1.003-(1.479*(10^(-2))*Tavg1)+(1.982*(10^(-5))*(Tavg1^2)));
  MUw=MUw1/1488.16;%lbm/ft.s
  %Bw
  TX=Tavg1-60;
  Bw=1+1.2*(10^(-4))*TX+1*(10^(-6))*(TX^2)-3.33*(10^(-6))*Pavg;
  %ROw
  ROw=ROwsc/Bw;
qg=((GLR-Rs)*qo_prim*Bg)/(24*3600);%ft^3/s
g=32.174;
gc=32.174;
qo=qo_prim*Bo;%bbl/day
qw=qw_prim*Bw;%bbl/day
ql=(qo+qw)*(5.617/(24*3600));%ft^3/s
S=(pi/4)*d1^2;%ft^2
Vsg=qg/S;%ft/s
Vsl=ql/S;%ft/s
fo=1/(1+WOR);
fw=1-fo;
ROl=ROo*fo+ROw*fw;%lbm/ft^3
MUl=MUo*fo+MUw*fw;
ZG=ZG_o*fo+ZG_w*fw;
dpdz=dp_ansari_function(a1,d1,MUg,ROg,Vsg,Vsl,ROl,MUl,ZG,g,gc);
Pt(i+1)=Pavg-dpdz*dL;
Pbh=Pt(i+1);
end
H=(N:-1:0).*dL;
plot(Pt,H);
xlabel('Pressure(psia)');
ylabel('Measured Depth(ft)');