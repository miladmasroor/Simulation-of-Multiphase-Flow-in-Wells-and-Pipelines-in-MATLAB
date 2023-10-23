clc;
clear;
disp('shoma darid az raveshe lockhartmartinelli estefade mikonid.');
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

%dar in marhale mikhahim meghdar X ra be dast biyavarim.
x=(Vsg*ROg)/(Vsg*ROg+Vsl*ROl);
X=((((1-x)/x)^1.8)*(ROg/ROl)*((MUl/MUg)^0.2))^0.5;

%hala mikhahim meghdar liquid holdup ra bedast biyavarim.
fg=(1+(X)^0.8)^(-0.378);
fl=1-fg;
%hala mikhahim ba estefade az liquid holdup ROs ra be dast biyavarim.
ROs=(ROl*fl)+(ROg*(1-fl));

%dar in marhale mirim soraghe mohasebeye oft feshar.
Rel=(ROl*Vsl*d1)/MUl1;
Reg=(ROg*Vsg*d1)/MUg1;
fliquid=0.32*(Rel)^(-0.25); %Blasius
if Rel>1000 && Reg>1000
    c=20;
elseif Rel>1000 && Reg<1000;
    c=10;
else
    c=5;
end
gc=32.174;
PHIl2=(1+(c/X)+1/(X)^2);
dpdzf=((fliquid*ROl*(Vsl)^2)/(2*gc*d1))*PHIl2; %vahed=psf/ft.
dpdzH=ROs; %vahed=psf/ft.
dpdztotal=dpdzf+dpdzH; %vahed=psf/ft.
dpdztotal1=dpdztotal/148.73; %vahed=psi/ft.
disp(['(dp/dz)total ba estefade az raveshe lockhartmartinelli barabar ast ba(psi/ft)=' num2str(dpdztotal1)]);


















