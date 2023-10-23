clc;
clear;
close all;
ROwsc=67;%lbm/ft^3 
API=16.01;
R=10.732;
Pbh=1014.7;%psia
d1=1.995;%in
d=d1/12;%ft
qo_prime=400;%stb/day
qw_prime=600;%stb/day
%qg_prime=100000000;%scf/day
L=2000;%ft
n=10;
dL=L/n;%ft
WOR=1.5;
GLR=500;%scf/STB
gamma_g=0.65;
Tsc=520;%R
Psc=14.7;%psia
%Mg=gamma_g*28.97;


for i=1:n

%i=1;
   T1=120+0.017*(dL *(i-1));
   T2=120+0.017*(dL*i);
   Tavg=(T1+T2)/2+460;%R
   %Tavg=80+460;%R
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
   MUg=VisLeeGonzalez(Mg,Tavg,ROg1);
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
     MUo=AAA*(MUoD^BBB);
  else
     ZZ=3.0324-0.02023*API;
     yy=10^ZZ;
     XX=yy*(Tavg1^(-1.163));
     MUoD=(10^XX)-1;%dead oil viscosity
     AAA=10.715*((Rs+100)^(-0.515)) ;
     BBB=5.44*((Rs+150)^(-0.338));
     MUoB=AAA*(MUoD^BBB);
     m=2.6*(Pavg^1.187)*exp(-11.513-(8.98*10^(-5))*Pavg) ;
     MUo=MUoB*((Pavg/PB)^m);
  end
  %MUw
  MUw=exp(1.003-(1.479*(10^(-2))*Tavg1)+(1.982*(10^(-5))*(Tavg1^2)));
  %Bw
  TX=Tavg1-60;
  Bw=1+1.2*(10^(-4))*TX+1*(10^(-6))*(TX^2)-3.33*(10^(-6))*Pavg;
  
end

