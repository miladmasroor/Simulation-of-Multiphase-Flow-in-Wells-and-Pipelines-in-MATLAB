clc;
clear;
a=input('input deviation angle from horizontal(degrees):');
a1=(a*pi)/180;
MUg1=input('input gas viscozity(cp):');
MUg=MUg1*0.001;%pa.s
MUo1=input('input oil viscozity(cp):');
MUo=MUo1*0.001;%pa.s
MUw1=input('input water viscozity(cp):');
MUw=MUw1*0.001;%pa.s
d=input('input diameter(in):');
d2=d/12;
d1=d*0.0254;%m
WOR=input('input water oil ratio:');
qo_prim=input('input oil rate(STB/day):');
qw_prim=input('input water rate (STB/day):');
qg=input('input gas rate(ft^3/s):');
ROo1=input('input oil density(lbm/ft^3):');
ROo=ROo1*16.018;%kg/m^3
ROg1=input('input gas density(lbm/ft^3):');
ROg=ROg1*16.018;%kg/m^3
ROw1=input('input water density(lbm/ft^3):');
ROw=ROw1*16.018;%kg/m^3
ZG_o1=input('input oil surface tension(dynes/cm):');
ZG_o=ZG_o1*0.001;%N/m
ZG_w1=input('input water surface tension(dynes/cm):');
ZG_w=ZG_w1*0.001;%N/m
Bo=input('input oil FVF(bbl/STB):');
Bw=input('input water FVF(bbl/STB):');
g=9.81;
qo=qo_prim*Bo;%bbl/day
qw=qw_prim*Bw;%bbl/day
ql=(qo+qw)*(5.617/(24*3600));%ft^3/s
S=(pi/4)*d2^2;
Vsg1=qg/S;
Vsg=Vsg1*0.3048;%m/s
Vsl1=ql/S;
Vsl=Vsl1*0.3048;%m/s
fo=1/(1+WOR);
fw=1-fo;
ROl=ROo*fo+ROw*fw;%kg/m^3
MUl=MUo*fo+MUw*fw;
ZG=ZG_o*fo+ZG_w*fw;
E=1-exp((-0.125)*((10^4*Vsg*MUg/ZG)*((ROg/ROl)^0.5)-1.5));
Elc=((Vsl*E)/(Vsg+(Vsl*E)));
ROc=ROl*Elc+ROg*(1-Elc);
MUc=MUl*Elc+MUg*(1-Elc);
NResl=(ROl*Vsl*d1)/(MUl);
NReF=(ROl*Vsl*d1*(1-E))/(MUl);
if NResl<2100
    fsl=64/(NResl);
elseif NResl>50000
    fsl=0.184*(NResl)^(-0.2);
else
    fsl=0.32*(NResl)^(-0.25);
end
if NReF<2100
    fF=64/(NReF);
elseif NReF>50000
    fF=0.184*(NReF)^(-0.2);
else
    fF=0.32*(NReF)^(-0.25);
end
dpdzsl=(fsl*ROl*(Vsl)^2)/(2*d1); 
Vsc=Vsg+Vsl*E;
NRec=(ROc*Vsc*d1)/MUc;
if NRec<2100
    fc=64/(NRec);
elseif NRec>50000
    fc=0.184*(NRec)^(-0.2);
else
    fc=0.32*(NRec)^(-0.25);
end
dpdzc=(fc*ROc*(Vsc)^2)/(2*d1);

xm=((1-E)^2)*(fF/fsl)*(dpdzsl/dpdzc);%manzoram az xm, xm^2 mibashad.
ym=((ROl-ROg)*g)/dpdzc;
%tayin flf 
%estefade az raveshe newton_rabson
%fi/fc=1+24*(pl/pg)^1/3*deltahad
%flf=1-(1-2*deltahad)^2
Flf0 = 0.5 ; %hadse avaliye
    error = 1;
    while abs(error)>=1e-5
        delta = (1-sqrt(1-Flf0))/2;
        fifc=1+24*(ROl/ROg)^(1/3)*delta; 
        syms x
        f=ym-fifc/((1-x)^2.5*x)+xm/x^3;
        dif_f=diff(f);
        f=inline(f);
        dif_f=inline(dif_f);
        Flf=Flf0-feval(f,Flf0)/feval(dif_f,Flf0);
        error=Flf-Flf0;
        Flf0=Flf;
    end



Vma=Vsl+Vsg; %manzor az Vma Vmactual mibashad.
NRe=(ROl*Vma*d1)/MUl;
if NRe<2100
    ff=64/(NRe);
elseif NRe>50000
    ff=0.184*(NRe)^(-0.2);
else
    ff=0.32*(NRe)^(-0.25);
end
syms Vm
B=2*Vm^1.2*(ff/(2*d1))^0.4*(ROl/ZG)^0.6*sqrt(0.4*ZG/(g*(ROl-ROg)))-0.725-4.15*sqrt(Vsg/Vm);
dif_B=diff(B);
B=inline(B);
dif_B=inline(dif_B);
Vm0=Vma;
    while abs(feval(B,Vm0))>=1e-5
        Vm_calc=Vm0-feval(B,Vm0)/feval(dif_B,Vm0);
        Vm0=Vm_calc;
    end
Vb=1.53*(((g*(ROl-ROg)*ZG)/(ROl)^2)^0.25);%manzor az Vb Vbinahayat mibashad.
VbT=(0.35*sin(a1)+0.54*cos(a1))*((g*d1*(ROl-ROg)/ROl)^0.5);

if Vsg>3.1*((ZG*g*(ROl-ROg)/(ROg)^2)^0.25)
    if (Flf+Elc*(1-Flf))<0.12
        disp('we have annular flow.');
    else
        if Vma>Vm_calc
            if Vsg<1.08*Vsl
                disp('we have dispersed flow.');
            else
                if Vsg<=(0.429*Vsl+0.357*Vb)*sin(a1)
                    if d1>19.1*((ZG/(g*(ROl-ROg)))^0.5)
                        disp('we have bubbly flow. ');
                else
                    if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end                      
                    end
                else
                    if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end 
             
                end
            end
        else
            if Vsg<=(0.429*Vsl+0.357*Vb)*sin(a1)
                    if d1>19.1*((ZG/(g*(ROl-ROg)))^0.5)
                        disp('we have bubbly flow. ');     
                else
                    if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end  
                    end
            else
                if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                end 
                    
            end
        
        end
    end
else 
    if Vma>Vm_calc
            if Vsg<1.08*Vsl
                disp('we have dispersed flow.');
            else
                if Vsg<=(0.429*Vsl+0.357*Vb)*sin(a1)
                   if d1>19.1*((ZG/(g*(ROl-ROg)))^0.5)
                       disp('we have bubbly flow. ');
                else
                   if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end 
                   end
                else
                   if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end 
                end
            end
    else
        if Vsg<=(0.429*Vsl+0.357*Vb)*sin(a1)
                   if d1>19.1*((ZG/(g*(ROl-ROg)))^0.5)
                       disp('we have bubbly flow. ');
                else
                    if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end 
                   end
                else
                   if Vsg>12.19*(1.2*Vsl+VbT)
                        disp('we have churn flow.');
                    else
                        disp('we have slug flow. ');   
                    end 
        end
    end
end


        
        
        

















   
