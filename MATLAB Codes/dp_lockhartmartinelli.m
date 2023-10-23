function[dpdztotal1]=dp_lockhartmartinelli(Vsl,MUl,Vsg,MUg,ROl,d1,gc,ROg)
%dar in marhale mikhahim meghdar X ra be dast biyavarim.
x=(Vsg*ROg)/(Vsg*ROg+Vsl*ROl);
X=((((1-x)/x)^1.8)*(ROg/ROl)*((MUl/MUg)^0.2))^0.5;

%hala mikhahim meghdar liquid holdup ra bedast biyavarim.
fg=(1+(X)^0.8)^(-0.378);
fl=1-fg;
%hala mikhahim ba estefade az liquid holdup ROs ra be dast biyavarim.
ROs=(ROl*fl)+(ROg*(1-fl));

%dar in marhale mirim soraghe mohasebeye oft feshar.
Rel=(ROl*Vsl*d1)/MUl;
Reg=(ROg*Vsg*d1)/MUg;
fliquid=0.32*(Rel)^(-0.25); %Blasius
if Rel>1000 && Reg>1000
    c=20;
elseif Rel>1000 && Reg<1000;
    c=10;
else
    c=5;
end
%gc=32.174;
PHIl2=(1+(c/X)+1/(X)^2);
dpdzf=((fliquid*ROl*(Vsl)^2)/(2*gc*d1))*PHIl2; %vahed=psf/ft.
dpdzH=ROs; %vahed=psf/ft.
dpdztotal=dpdzf+dpdzH; %vahed=psf/ft.
dpdztotal1=dpdztotal/148.73; %vahed=psi/ft.
disp(['(dp/dz)total ba estefade az raveshe lockhartmartinelli barabar ast ba(psi/ft)=' num2str(dpdztotal1)]);
end