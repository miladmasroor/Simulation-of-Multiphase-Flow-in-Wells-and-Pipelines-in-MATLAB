clc;
clear;
close all;
GLR=input('input gas liquid ratio(scf/STB):');
Pbh1=input('input bottom hole pressure(psig):');
Pbh=Pbh1+14.7;
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
i=1;
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
%MUg1=input('input gas viscozity(cp):');
%MUg=MUg1/1488.16;%lbm/ft.s
%MUo1=input('input oil viscozity(cp):');
%MUo=MUo1/1488.16;%lbm/ft.s
%MUw1=input('input water viscozity(cp):');
%MUw=MUw1/1488.16;%lbm/ft.s
%d=input('input diameter(in):');
%d1=d/12;%ft
%WOR=input('input water oil ratio:');
%qo_prim=input('input oil rate(STB/day):');
%qw_prim=input('input water rate (STB/day):');
%qg=input('input gas rate(ft^3/s):');
%ROo=input('input oil density(lbm/ft^3):');
%ROg=input('input gas density(lbm/ft^3):');
%ROw=input('input water density(lbm/ft^3):');
%ZG_o1=input('input oil surface tension(dynes/cm):');
%ZG_o=ZG_o1/453.632;%lb/s^3
%ZG_w1=input('input water surface tension(dynes/cm):');
%ZG_w=ZG_w1/453.632;%lb/s^3
%Bo=input('input oil FVF(bbl/STB):');
%Bw=input('input water FVF(bbl/STB):');
%R=10.732;
qg=((GLR-Rs)*qo_prim*Bg*5.617)/(24*3600);%ft^3/s
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

landa=ql/(ql+qg);
Vm=Vsl+Vsg;%ft/s
NFr=(Vm^2)/(g*d1);

L1=316*(landa^0.302);
L2=0.0009252*(landa^(-2.4684));
L3=0.10*(landa^(-1.4516));
L4=0.5*(landa^(-6.738));

if landa<0.01 && NFr<L1 || landa>=0.01 && NFr<L2
    disp('we have segregated flow.');
elseif landa>=0.01 &&  L2<=NFr && NFr<=L3
    disp('we have transition flow.');
elseif landa>=0.01 && landa<0.4 &&  L3<NFr && NFr<=L1 || landa>=0.4 && L3<NFr && NFr<=L4
     disp('we have intermittent flow.');
else
     disp('we have distributed flow.');
end
n=input('what is your flow regime(1.intermittent / 2.segregated / 3.distributed /4.transition)=');
switch n
    case 1
Nlv=1.938*Vsl*((ROl/ROg)^0.25);
HL0=0.845*(landa^0.5351)/(NFr^0.0173);
m=input('what is your flow direction(1.up / 2.down )=');
switch m
    case 1
C=(1-landa)*(log((2.69*(landa^0.305)*(NFr^0.0978))/(Nlv^0.4473)));
    case 2
C=(1-landa)*(log((4.7*(Nlv^0.1244))/((NFr^0.5056)*(landa^0.3692))));
end    
Z=1+C*(sin(a1)-(1/3)*((sin(a1))^3));%z=inclination correction factor.
HLa=HL0*Z;

ROn=ROl*landa+ROg*(1-landa);
MUn=MUl*landa+MUg*(1-landa);
NRe=ROn*Vm*d1/MUn;
fn=1/((2*log10(NRe/(4.5223*(log10(NRe))-3.8215)))^2);
y=landa/(HLa^2);
if y>1 && y<1.2
    s=log(2.2*y-1.2);
else
    s=(log(y))/(3.182*(log(y))-0.8728*((log(y))^2)+0.01853*((log(y))^4)-0.0523);
end
ftpfn=exp(s);
ftp=fn*ftpfn;

ROtp=ROl*HLa+ROg*(1-HLa);
dpdz=((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1);
disp(['dpdz(psf/ft):' num2str(dpdz)]);

    case 2
Nlv=1.938*Vsl*((ROl/ROg)^0.25);
HL0=0.98*(landa^0.4846)/(NFr^0.0868);
m=input('what is your flow direction(1.up / 2.down )=');
switch m
    case 1
C=(1-landa)*(log((0.011*(Nlv^3.539))/((NFr^1.614)*(landa^3.768))));
    case 2
C=(1-landa)*(log((4.7*(Nlv^0.1244))/((NFr^0.5056)*(landa^0.3692))));
end    
Z=1+C*(sin(a1)-(1/3)*((sin(a1))^3));%z=inclination correction factor.
HLa=HL0*Z;  

ROn=ROl*landa+ROg*(1-landa);
MUn=MUl*landa+MUg*(1-landa);
NRe=ROn*Vm*d1/MUn;
fn=1/((2*log10(NRe/(4.5223*(log10(NRe))-3.8215)))^2);
y=landa/(HLa^2);
if y>1 && y<1.2
    s=log(2.2*y-1.2);
else
    s=(log(y))/(3.182*(log(y))-0.8728*((log(y))^2)+0.01853*((log(y))^4)-0.0523);
end
ftpfn=exp(s);
ftp=fn*ftpfn;

ROtp=ROl*HLa+ROg*(1-HLa);
dpdz=((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1);
disp(['dpdz(psf/ft):' num2str(dpdz)]);
    case 3
Nlv=1.938*Vsl*((ROl/ROg)^0.25);
HL0=1.065*(landa^0.5824)/(NFr^0.0609);
m=input('what is your flow direction(1.up / 2.down )=');
switch m
    case 1
C=0;
    case 2
C=(1-landa)*(log((4.7*(Nlv^0.1244))/((NFr^0.5056)*(landa^0.3692))));
end    
Z=1+C*(sin(a1)-(1/3)*((sin(a1))^3));%z=inclination correction factor.
HLa=HL0*Z;

ROn=ROl*landa+ROg*(1-landa);
MUn=MUl*landa+MUg*(1-landa);
NRe=ROn*Vm*d1/MUn;
fn=1/((2*log10(NRe/(4.5223*(log10(NRe))-3.8215)))^2);
y=landa/(HLa^2);
if y>1 && y<1.2
    s=log(2.2*y-1.2);
else
    s=(log(y))/(3.182*(log(y))-0.8728*((log(y))^2)+0.01853*((log(y))^4)-0.0523);
end
ftpfn=exp(s);
ftp=fn*ftpfn;

ROtp=ROl*HLa+ROg*(1-HLa);
dpdz=((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1);
disp(['dpdz(psf/ft):' num2str(dpdz)]);
    case 4
        A=(L3-NFr)/(L3-L2);
        B=1-A;
        Nlv=1.938*Vsl*((ROl/ROg)^0.25);
        HL0int=0.845*(landa^0.5351)/(NFr^0.0173);
        HL0seg=0.98*(landa^0.4846)/(NFr^0.0868);
        i=input('what is your flow direction(intermittent)(1.up / 2.down )=');
switch i
    case 1
C=(1-landa)*(log((2.69*(landa^0.305)*(NFr^0.0978))/(Nlv^0.4473)));
    case 2
C=(1-landa)*(log((4.7*(Nlv^0.1244))/((NFr^0.5056)*(landa^0.3692))));
end   
w=input('what is your flow direction(segregated)(1.up / 2.down )=');
switch w
    case 1
C=(1-landa)*(log((0.011*(Nlv^3.539))/((NFr^1.614)*(landa^3.768))));
    case 2
C=(1-landa)*(log((4.7*(Nlv^0.1244))/((NFr^0.5056)*(landa^0.3692))));
end
Z=1+C*(sin(a1)-(1/3)*((sin(a1))^3));%z=inclination correction factor.
HLaint=HL0int*Z;
HLaseg=HL0seg*Z;
HLatrans=A*HLaseg+B*HLaint;

ROn=ROl*landa+ROg*(1-landa);
MUn=MUl*landa+MUg*(1-landa);
NRe=ROn*Vm*d1/MUn;
fn=1/((2*log10(NRe/(4.5223*(log10(NRe))-3.8215)))^2);
y=landa/(HLatrans^2);
if y>1 && y<1.2
    s=log(2.2*y-1.2);
else
    s=(log(y))/(3.182*(log(y))-0.8728*((log(y))^2)+0.01853*((log(y))^4)-0.0523);
end
ftpfn=exp(s);
ftp=fn*ftpfn;

ROtp=ROl*HLatrans+ROg*(1-HLatrans);
dpdz=((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1);
disp(['dpdz(psf/ft):' num2str(dpdz)]);

end


















