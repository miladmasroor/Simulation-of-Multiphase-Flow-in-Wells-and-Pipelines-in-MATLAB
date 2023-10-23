function[dpdztotal1]=dp_homogeneus_function_Q1(Vsl,MUl,Vsg,MUg,ROl,d1,gc,ROg,g,a1)
x=(Vsg*ROg)/(Vsg*ROg+Vsl*ROl);
MUs=(MUl*(1-x))+MUg*x; 
Vm=Vsg+Vsl;
fg=Vsg/Vm;
ros=((1-x)/ROl)+x/ROg;
ROs=1/ros;
Gm=(Vsg*ROg)+(Vsl*ROl);
Rem=(Gm*d1)/MUs;
fm=0.32*(Rem)^(-0.25);
dpdzf=(fm*((Vm)^2)*ROs)/(2*gc*d1);%vahed=psf/ft.
dpdzH=ROs*g*sin(a1)/gc;%vahed=psf/ft.
dpdztotal=dpdzf+dpdzH; %vahed=psf/ft.
dpdztotal1=dpdztotal/148.73; %vahed=psi/ft.
disp(['(dp/dz)total ba estefade az raveshe homogeneus barabar ast ba(psi/ft)=' num2str(dpdztotal1)]);
end