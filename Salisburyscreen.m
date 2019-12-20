%Matlab code that calculates the wavelength of destructive interference in a Salisbury screen for a given spacer width:
% ********************************************************* 
%Testing now:
% ********************************************************* 
%write
% ********************************************************* 
% Salisbury screen
% ********************************************************* 
clc
clear all;
% ********************************************************* 
% Code options
% ********************************************************* 
test=1;
set=1;% 1 change ds, 2 change dig, 3 ns
datac=1;%1 enables data save 
nig=0;%refractive index of ion gel
ns=1.8;%refractive index of spacer 
dig=0;%thickness of ion gel in um 
ds=0:0.01:0.462;%thickness of ion gel in um 
nn=1:2:6;%order
lrange=1:0.1:2;%wavelength range in um
% ********************************************************* % Output file names
% *********************************************************
if set==1
    nig=0;%refractive index of ion gel
    ns=1.8;%refractive index of spacer
    dig=0;%thickness of ion gel in um
    l=zeros(length(ds),length(nn));
    for i=1:length(ds)
        for j=1:length(nn)
            l(i,j)=((4*ns*ds(i)+4*nig*dig)/nn(j));
        end
    end
figure;
GraphTitle = ('Salisbury screen central wavelength 1'); plot(l,ds)
xlabel ('Wavelength (um)','FontSize',20);
ylabel ('Spacer thickness (um)','FontSize',20); title(GraphTitle,'FontSize',20);
saveas (gcf,GraphTitle,'fig');
%cd('data')
    if datac==1
        data=table(l,ds');
        writetable(data)% save data
    end
elseif set==2
    nig=1.42;%refractive index of ion gel
    ns=1.42;%refractive index of spacer
    ds=0.12;
    l=zeros(length(dig),length(nn));
    for i=1:length(dig)
        for j=1:length(nn)
            l(i,j)=((4*ns*ds+4*nig*dig(i)))/nn(j);
        end
    end
figure;
GraphTitle = ('Salisbury screen central wavelength 2'); plot(l,dig)
xlabel ('Wavelength (um)','FontSize',20);
ylabel ('Ion gel thickness (um)','FontSize',20); title(GraphTitle,'FontSize',20);
%saveas (gcf,GraphTitle,'fig');
    if datac==1
        data=table(l,dig');
        writetable(data)% save data
    end
    else
    for i=1:length(ns)
        for j=1:length(nn)
            l(i,j)=((4*ns(i)*ds+4*nig*dig))/nn(j);
        end
    end
figure;
GraphTitle = ('Salisbury screen central wavelength 3'); plot(l,ns)
xlabel ('Wavelength (um)','FontSize',20);
ylabel ('Spacer refractive index (um)','FontSize',20); title(GraphTitle,'FontSize',20);
saveas (gcf,GraphTitle,'fig');
    if datac==1
        data=table(l,ns');
        writetable(data)% save data
    end
end
if test==5
clear all
nSiO2=1/1.42;
nITO=1/1.8;
thickness_SiO2=0.015;
Thickness_ITO=nITO.*thickness_SiO2/nSiO2;
plot_range=0:0.00625:0.45;
plot_range2=plot_range+Thickness_ITO;
end