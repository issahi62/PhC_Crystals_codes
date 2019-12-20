clc
clear all
%input
variable=-3.5:0.5:6;
wavelength=0.6:0.005:0.9;
toggle_r=2; % for 0 r(0.0), for 1 total reflection, for 2 absorption 
cell_number=4; %how naby cells in file
data_01='Voltage_';
data_02='Voltage_';

if toggle_r==0
    tr=1;
elseif toggle_r==1
    tr=2;
else
tr=3;
end
%load
cd('data1')
fid=fopen([data_01 num2str(0) '.dat'],'r'); delete_header =fgetl(fid);
x=fscanf(fid,'%f');
fclose(fid); 
absorbance=zeros(max(length(x)/cell_number),1); 
wavelength01=zeros(max(length(x)/cell_number),1); 
j=1;
for i=1:cell_number:max(length(x))
    wavelength01(j,1)=x(i);
    absorbance(j,1)=x(i+tr);
    j=j+1;
end
for i=1:length(variable)-1
fid=fopen([data_01 num2str(i) '.dat'],'r');
delete_header =fgetl(fid); x=fscanf(fid,'%f');
fclose(fid);
    j=1;
    for k=1:cell_number:max(length(x))
        absorbance(j,i+1)=x(k+tr);
j=j+1;
    end
end