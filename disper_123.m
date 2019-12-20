%B.DispersionrelationinMATLAB
%B.1.AngleresolvedsimulationbyBruggemanandMaxwell-GarnettModel
%AngleresolvedsimulationbyBruggemanModel

%% section 1
e3=1;
eair=1;
e1=2.25;
d=120;
ed=eair;
c=3*10^17;%veloflightinnm/s
g1=42*pi/180;
g2=60*pi/180;
s=1*pi/180;
theta=g1:s:g2;
hold on
%fid=fopen('e_per_lambda.txt','a+');
epsinfinity=1.53;
lambdap=145;
gammap=17000;
A1=.94;
phi1=-pi/4;
lambda1=468;
gamma1=2300;
A2=1.36;
phi2=-pi/4;
lambda2=331;
gamma2=940;
for ii=1:340
    lambdaa(ii)=200+ii*(1400-200)/340;
    lambda=lambdaa(ii);
    epsilon(ii)=epsinfinity-1/lambdap^2/(1/lambda^2+sqrt(-1)/gammap/lambda)+A1/lambda1*((exp(sqrt(-1)*phi1)/(1/lambda1-1/lambda-sqrt(-1)/gamma1)...
        +(exp(-sqrt(-1)*phi1)/(1/lambda1+1/lambda+sqrt(-1)/gamma1))))...
        +A2/lambda2*((exp(sqrt(-1)*phi2)/(1/lambda2-1/lambda-sqrt(-1)/gamma2)...
        +(exp(-sqrt(-1)*phi2)/(1/lambda2+1/lambda+sqrt(-1)/gamma2))));
    e2=epsilon(ii);

    w=2*pi*c/lambda;%omega
    f=0.5;%f=1means100%puregold
    M=(3*f-1)*e2+(2-3*f)*e3;
    v=(M^2)+8*e2*e3;
    vol=imag(v^0.5);
    %mol=imag(M);
    if vol<0
    e2brug=(M-(v^0.5))/4;
    else
    e2brug=(M+(v^0.5))/4;
    end
    kx=(w/c)*((e1)^0.5)*sin(g1);%+ve
    kz1=(((w^2/c^2)*e1)-kx^2)^0.5;%(isrealnumber)
    kz2=(((w^2/c^2)*e2brug)-kx^2)^0.5;%(isaComplexno.)
    kz3=(((w^2/c^2)*e3)-kx^2)^0.5;%(shouldbepurelyimaginary)
    A=(1-((kz2*e1)/(kz1*e2brug)))-(1+((kz2*e1)/(kz1*e2brug)))*(kz3*e2brug-kz2*e3)/(kz3*e2brug+kz2*e3)*exp(2i*kz2*d);
    B=(1+((kz2*e1)/(kz1*e2brug)))-(1-((kz2*e1)/(kz1*e2brug)))*(kz3*e2brug-kz2*e3)/(kz3*e2brug+kz2*e3)*exp(2i*kz2*d);
    Hry=A/B;
    y=abs(Hry);
    x=lambda;
%fprintf(fid,'%E\n',x);
    plot(x,y);
end%fclose(fid);

% %% section 2 
% 
% d=60;%thicknessoflayerinnmc=3*10^17;%veloflightinnm/se1=2.25;%epsilonforglass
% eair=1;e3=eair;ed=eair;%epsilont1=42*pi/180;
% t2=57*pi/180;s=1*pi/180;
% 
% 
% theta=t1;%:s:t2;
% f=0.50;%percentageofair
% hold on
% epsinfinity=1.53;
% lambdap=145;
% gammap=17000;
% A1=.94;
% phi1=-pi/4;
% lambda1=468;
% gamma1=2300;
% A2=1.36;
% phi2=-pi/4;
% lambda2=331;
% gamma2=940;
% for ii=1:540
% lambdaa(ii)=496.28+ii*(2000-496.28)/540;
% lambda=lambdaa(ii);
% epsilon(ii)=epsinfinity-1/lambdap^2/(1/lambda^2+sqrt(-1)/gammap/lambda)+...
% A1/lambda1*((exp(sqrt(-1)*phi1)/(1/lambda1-1/lambda-sqrt(-1)/gamma1)+(exp(-sqrt(-1)*phi1)/(1/lambda1+1/lambda+sqrt(-1)/gamma1))))+...
% A2/lambda2*((exp(sqrt(-1)*phi2)/(1/lambda2-1/lambda-sqrt(-1)/gamma2)+(exp(-sqrt(-1)*phi2)/(1/lambda2+1/lambda+sqrt(-1)/gamma2))));egold=epsilon(ii);
% e_ll=f*eair+(1-f)*egold;e_per=egold*(((1+f)*eair+(1-f)*egold)/((1-f)*eair+(1+f)*egold));
% w=2*pi*c/lambda;%omegaforg=1:length(theta);
% g1=theta(g);
% kx=(w/c)*((e1)^0.5)*sin(g1);kz3=(((w^2/c^2)*ed)-kx^2)^0.5;kz2=(((w^2/c^2)*e_per)-(e_per/e_ll)*kx^2)^0.5;kz1=(((w^2/c^2)*e1)-kx^2)^0.5;%(isrealnumber)
% A=(1-((kz2*e1)/(kz1*e_per)))-(1+((kz2*e1)/(kz1*e_per)))*(kz3*e_per-kz2*e3)/(kz3*e_per+kz2*e3)*exp(2i*kz2*d);
% B=(1+((kz2*e1)/(kz1*e_per)))-(1-((kz2*e1)/(kz1*e_per)))*(kz3*e_per-kz2*e3)/(kz3*e_per+kz2*e3)*exp(2i*kz2*d);
% Hry=A/B;x=lambda;y=abs(Hry);
% plot(x,y);
% end
% 
% 
% w=2*pi*c/lambda;
% f=0.5;%f=1means100%puregold
% M=(3*f-1)*e2+(2-3*f)*e4;v=(M^2)+8*e2*e4;vol=imag(v^0.5);
% if vol<0
%     e2brug=(M-(v^0.5))/4;
% else
%     e2brug=(M+(v^0.5))/4;
% end
% 
% %Definingvariables
% d=120;%thicknessoflayerinnm
% c=3*10^17;%veloflightinnm/s%e1=2.25;%epsilonforGee3=1;%epsilonvalueofSi12.25
% e4=1;
% hold on
% %fid=fopen('Au60min_x2.txt','a+');
% epsinfinity=1.53;lambdap=145;gammap=17000;A1=.94;phi1=-pi/4;lambda1=468;
% gamma1=2300;A2=1.36;phi2=-pi/4;lambda2=331;gamma2=940;
% for ii=1:300
% lambdaa(ii)=400+ii*(1900-400)/300;
% lambda=780;%lambdaa(ii);
% epsilon(ii)=epsinfinity-
% 1/lambdap^2/(1/lambda^2+sqrt(-1)/gammap/lambda)+
% A1/lambda1*((exp(sqrt(-1)*phi1)/(1/lambda1-1/lambda-sqrt(-1)/gamma1)+(exp(-sqrt(-1)*phi1)/(1/lambda1+1/lambda+sqrt(-1)/gamma1))))+
% A2/lambda2*((exp(sqrt(-1)*phi2)/(1/lambda2-1/lambda-sqrt(-1)/gamma2)+(exp(-sqrt(-1)*phi2)/(1/lambda2+1/lambda+sqrt(-1)/gamma2))));e2=complex(real(epsilon(ii)),imag(epsilon(ii)));
% 
% A=(e2brug*e3)/(e3+e2brug);kx=(w/c)*(A^0.5)*10^3;
% x=real(kx)
% y=1240.7/lambda;
% %fprintf(fid,'%E\n',x);
% plot(x,y);end
% %fclose(fid);
% WavevectorbyMaxwell-GarnettModel
% eair=1;
% f=0.5;%percentageofair
% ed=eair;
% c=3*10^17;%veloflightinnm/shold on
% %fid=fopen('e_per_lambda.txt','a+');
% epsinfinity=1.53;lambdap=145;gammap=17000;A1=.94;phi1=-pi/4;lambda1=468;
% gamma1=2300;A2=1.36;phi2=-pi/4;lambda2=331;gamma2=940;
% for ii=1:340
% lambdaa(ii)=200+ii*(1900-200)/340;
% lambda=780;%lambdaa(ii);
% epsilon(ii)=epsinfinity-
% 1/lambdap^2/(1/lambda^2+sqrt(-1)/gammap/lambda)+
% A1/lambda1*((exp(sqrt(-1)*phi1)/(1/lambda1-1/lambda-sqrt(-1)/gamma1)+(exp(-sqrt(-1)*phi1)/(1/lambda1+1/lambda+sqrt(-1)/gamma1))))+
% A2/lambda2*((exp(sqrt(-1)*phi2)/(1/lambda2-1/lambda-sqrt(-1)/gamma2)+(exp(-sqrt(-1)*phi2)/(1/lambda2+1/lambda+sqrt(-1)/gamma2))));egold=epsilon(ii);
% e_ll=f*eair+(1-f)*egold;e_per=egold*(((1+f)*eair+(1-f)*egold)/((1-f)*eair+(1+f)*egold));
% w=2*pi*c/lambda;%omegaA=ed-e_per;
% 100
% 
% B=ed-(e_ll*e_per/ed);
% kx=(w/c)*((e_ll*A/B)^0.5);kzd=(((w^2/c^2)*ed)-kx^2)^0.5;kzm=(((w^2/c^2)*e_per)-(e_per/e_ll)*kx^2)^0.5;
% %y=real(e_per);y=1240.7/lambda;x=real(kx)
% %fprintf(fid,'%E\n',x);plot(x,y);
% end%fclose(fid);
% B.3.TransmissionCalculation
% %Definingvariablesclearall
% d=120;%thicknessoflayerinnmc=3*10^17;%veloflightinnm/se1=1;%epsilonforglass
% e3=1;%epsilonvalueofairt1=0;
% theta=t1;g1=theta;%:s:t2;hold on
% fid=fopen('Trans_Au46_y.txt','a+');
% epsinfinity=1.53;lambdap=145;gammap=17000;A1=.94;phi1=-pi/4;
% lambda1=468;gamma1=2300;A2=1.36;phi2=-pi/4;lambda2=331;gamma2=940;for ii=1:700
% lambdaa(ii)=300+ii*(1800-300)/700;
% lambda=lambdaa(ii);
% epsilon(ii)=epsinfinity-1/lambdap^2/(1/lambda^2+sqrt(-1)/gammap/lambda)+
% A1/lambda1*((exp(sqrt(-1)*phi1)/(1/lambda1-1/lambda-sqrt(-1)/gamma1)+(exp(-sqrt(-1)*phi1)/(1/lambda1+1/lambda+sqrt(-1)/gamma1))))+
% A2/lambda2*((exp(sqrt(-1)*phi2)/(1/lambda2-1/lambda-sqrt(-1)/gamma2)+(exp(-sqrt(-1)*phi2)/(1/lambda2+1/lambda+sqrt(-1)/gamma2))));e2=complex(real(epsilon(ii)),imag(epsilon(ii)));
% 
% %forjj=1:100
% %theta(jj)=42+jj*(52-42)/100;
% %forg=1:length(theta);w=2*pi*c/lambda;%omega
% f=0.5;%f=1means100%puregoldM=(3*f-1)*e2+(2-3*f)*e3;
% v=(M^2)+8*e2*e3;vol=imag(v^0.5);
% %mol=imag(M);ifvol<0
% e2brug=(M-(v^0.5))/4;
% else
% e2brug=(M+(v^0.5))/4;end
% RE=(real(e2brug))^2;IE=(imag(e2brug))^2;
% UN=RE+IE;
% n=(0.5*(UN^0.5+RE^0.5))^0.5;s=(0.5*(UN^0.5-RE^0.5))^0.5;t1=2*3.5/(3.5+n);
% t3=2*n/(1+n);
% g1=0;%theta(g);kx=(w/c)*((e1)^0.5)*sin(g1);%+ve
% %kz1=(((w^2/c^2)*e1)-kx^2)^0.5;%(isrealnumber)
% kz2=(((w^2/c^2)*e2brug)-kx^2)^0.5;%(isaComplexno.)
% %kz3=(((w^2/c^2)*e3)-kx^2)^0.5;%(shouldbepurelyimaginary)Int=(t1*t3*exp(-imag(kz2)*d))^2;
% x=lambda;%g1;
% y=Int;
% %axis([300180000.2])%fprintf(fid,'%E\n',y);
% plot(x,y);%end
% end%fclose(fid);