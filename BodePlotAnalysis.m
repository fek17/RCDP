clear
clc
omega=logspace(-1,2,1000);
magnitude=0.9686./sqrt(power(omega,2)*power(2.056,2)+1);
argument=-(atan(0.2056*omega)+16.33*omega);
semilogx(omega,magnitude);
hold on
semilogx(omega,argument);
%wcg=11.09
%mag=0.043
