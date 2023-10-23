function [MUg1]=VisLeeGonzalez(Mg,Tavg,ROg1)
K=(9.4+0.02*Mg)*(Tavg^1.5)/(209+19*Mg+Tavg);
X=3.5+(986/Tavg)+0.01*Mg;
y=2.4-(0.2*X);
MUg1=K*(10^(-4))*exp(X*(ROg1^y));
end