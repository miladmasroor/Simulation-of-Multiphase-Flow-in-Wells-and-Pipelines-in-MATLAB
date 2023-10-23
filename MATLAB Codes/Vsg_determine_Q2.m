clc;
clear;
close all;

n=100;
a=linspace(0,44,n)';
x=zeros(n,1);
for i=1:n
    f=@ (x)  ((-35.5896)*x^2)+159.0244*x+21653+8593.6*sind(a(i));
    x(i)=fzero(f,5); 
    
end

figure;
plot(a,x,'LineWidth',2);
xlabel('deviation angle from horizontal');
ylabel(' superficial gas velocity (Vsg,m/s)');