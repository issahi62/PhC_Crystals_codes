%Matlab code that calculates the angular dependence of Rayleigh anomalies in diffraction gratings:
%Diffraciton grating
clc
clear all
%input
a_g=(700)*10^-9; % gratting latice constant in m (Warning! only use 1 or 2 steps if double length(v_f)>1 only give more than one value when length(theta=1))
n=1;%refractive index
la=0.4:0.001:1; %wavelength range in um
%thetag2=45*pi/180;
%wavelength convert
lambda=la*10^-6; %m
%Rayleigh 2D
i=0;
for int=-1:2:1
    for v_g=-4:1:4
            i=i+1;
thetar(i,:)=asin(int-lambda/a_g*v_g)/pi*180;
    end
end
%for int=-1:2:1
    %for v_g=-3:1:3
%i=i+1;
%thetar(i,:)=asin(int-lambda/(sqrt(a_g^2+a_g^2))*v_g/cos(the- tag2))/pi*180;
%end %end
figure;
GraphTitle = ('diffraction dispersion'); 
plot(lambda*10^6,real(thetar), 'Linewidth', 2.5)
xlabel ('Wavelength (um)','FontSize',20); 
ylabel ('angle (degrees)','FontSize',20); title(GraphTitle,'FontSize',20);
%saveas (gcf,GraphTitle,'jpg');
k=0;
 for i=1:2:2*length(thetar(:,1)')
    j=i+1;
    k=k+1;
    print1(:,i)=lambda';
    print2(:,j)=real(thetar(k,:)');
end