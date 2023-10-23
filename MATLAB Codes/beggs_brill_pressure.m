function [dpdz]=beggs_brill_pressure(a1,d1,ROg,qg,g,gc,Vsg,Vsl,ROl,MUl,ql,MUg,m,j,w)
landa=ql/(ql+qg);
Vm=Vsl+Vsg;%ft/s
NFr=(Vm^2)/(g*d1);

L1=316*(landa^0.302);
L2=0.0009252*(landa^(-2.4684));
L3=0.10*(landa^(-1.4516));
L4=0.5*(landa^(-6.738));

if landa<0.01 && NFr<L1 || landa>=0.01 && NFr<L2
    %disp('we have segregated flow.');
    n=2;
elseif landa>=0.01 &&  L2<=NFr && NFr<=L3
    %disp('we have transition flow.');
    n=4;
elseif landa>=0.01 && landa<0.4 &&  L3<NFr && NFr<=L1 || landa>=0.4 && L3<NFr && NFr<=L4
     %disp('we have intermittent flow.');
     n=1;
else
     %disp('we have distributed flow.');
     n=3;
end
%n=input('what is your flow regime(1.intermittent / 2.segregated / 3.distributed /4.transition)=');
switch n
    case 1
Nlv=1.938*Vsl*((ROl/ROg)^0.25);
HL0=0.845*(landa^0.5351)/(NFr^0.0173);
%m=input('what is your flow direction(1.up / 2.down )=');
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
dpdz=(((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1))/(148.73);
disp(['dpdz(psi/ft):' num2str(dpdz)]);

    case 2
Nlv=1.938*Vsl*((ROl/ROg)^0.25);
HL0=0.98*(landa^0.4846)/(NFr^0.0868);
%m=input('what is your flow direction(1.up / 2.down )=');
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
dpdz=(((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1))/(148.73);
disp(['dpdz(psi/ft):' num2str(dpdz)]);
    case 3
Nlv=1.938*Vsl*((ROl/ROg)^0.25);
HL0=1.065*(landa^0.5824)/(NFr^0.0609);
%m=input('what is your flow direction(1.up / 2.down )=');
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
dpdz=(((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1))/(148.73);
disp(['dpdz(psi/ft):' num2str(dpdz)]);
    case 4
        A=(L3-NFr)/(L3-L2);
        B=1-A;
        Nlv=1.938*Vsl*((ROl/ROg)^0.25);
        HL0int=0.845*(landa^0.5351)/(NFr^0.0173);
        HL0seg=0.98*(landa^0.4846)/(NFr^0.0868);
        %i=input('what is your flow direction(intermittent)(1.up / 2.down )=');
switch j
    case 1
C=(1-landa)*(log((2.69*(landa^0.305)*(NFr^0.0978))/(Nlv^0.4473)));
    case 2
C=(1-landa)*(log((4.7*(Nlv^0.1244))/((NFr^0.5056)*(landa^0.3692))));
end   
%w=input('what is your flow direction(segregated)(1.up / 2.down )=');
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
dpdz=(((g/gc)*ROtp*sin(a1))+ftp*ROn*(Vm^2)/(2*gc*d1))/(148.73);
%disp(['dpdz(psi/ft):' num2str(dpdz)]);

end
end