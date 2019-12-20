function [kx, ky, vek, lvek, w, ws, wint, dw, D_w, n, dimk, dim2, dim3] = dispersion_relation(dim1,ord)
c = 1/13; 
dim_temp = (ord-1)*ord/2 + (ord+1)*(ord+2)/2; 
dim = 2*dim_temp-1; 
a = 500e-09; 
A_cell = a^2*sqrt(3)/2; 
R = 0.48*a; 
A_hole = R^2*pi; 
A = zeros(dim, dim); 
B = zeros(dim, dim); C = zeros(dim, 1);
dim2 = floor(dim1/sqrt(3)); dim3 = floor(dim1*2/sqrt(3)); 
dimk = dim1+dim2+dim3;
ewp = zeros(dim, dim); 
gamma = zeros(dim, dim); beta=zeros(dim,dim); 
[kx, ky] = meshgrid(linspace(0,1, dim2), linspace(0, 1, dim1)); 
k = []; 
wv1 =[]; 
u =0;
vek3 = NaN(dimk, dim); 
vek4 = vek3; 

for ii = 0:ord
    for jj = 0:ii 
        if(ii>0)
            t = t((ii)*ii)/2+(1+jj);
        else
            t = 1;
        end 
        C(t) = jj;
        D(t) = ii-jj;
        if (jj>0 && t<ii)
            C(dim_temp-u) = -C(t);
            D(dim_temp-u) = D(t); 
            u = u+1; 
        end 
    end 
end 
C = [flipud(C(2:dim_temp)); -C(1:dim_temp)]; 
D = [flipud(D(2:dim_temp)); -D(1:dim_temp)]; 

wb = waitbar(0, ['working ... ', num2str(ceil(0/dim1*100)), '% done']);

for l = 1:dim1
    mend = round(l/sqrt(3)); 
    if mend==0
        mend =1; 
    end 
    if mend>dim2
       mend=dim2; 
    end 
    
    for m =1:mend 
        for ii = 1:dim 
            for jj = 1:dim 
                A(ii, jj) = C(ii)-C(jj); 
                B(ii, jj) = D(ii) - D(jj); 
                if((abs(A(ii,jj)) + abs(B(ii,jj))) <= ord) 
                    if (A(ii,jj) ==0 && B(ii,jj)== 0)
                        gamma(ii,jj) = c+ (1-c)*A_hole/A_cell; 
                    else 
                        G = 2*pi/a*sqrt(A(ii,jj)^2+(2*B(ii,jj)-A(ii,jj))^2/3); 
                        gamma(ii,jj) = (1-c)/A_cell*3*pi*R/G*besselj(1, G*R);
                    end
                end 
                beta(ii,jj) = ((C(ii)-A(ii,jj)+kx(l,m)/3)^2+((2*(D(ii)-B(ii,jj)-(C(ii)-A(ii,jj)))/sqrt(3)+ky(l,m)/sqrt(3))^2));
            end 
        end
        ewp = gamma.*beta; 
        vek1 = sort(sqrt(eig(ewp(:, :))));
        
        if l == 1 && m == 1 
            lvek = length(vek1); 
            vek2 = cell(1,lvek); 
        end 
        for n =1:lvek 
            vek2{n}(l,m) = vek1(n);
        end
        if kx(l,m) ==0
            vek3(l,:) = vek1; 
        elseif ky(l,m) ==1
            vek3(dim1-1+m, :) =vek1; 
        elseif m == mend
                vek2(dim1+dim2+dim3-2-l, :) = vek1;
        end
        [o,p] = find(isnan(vek3(1:dim1+dim2+dim3-100,:)));
        vek4 = [vek3(1:o(1)-1,:); vek3(o(end)+1:end, :)];
        waitbar(l/dim1, wb, ['working ... ', num2str(ceil(l/dim1*100)), '% done']); 
    end 
end 
close(wb); 
vek = vek2; 
w = vek4; 

for n = 1:lvek
    vek{n}(vek{n}==0) = NaN; 
    vek{n}(vek{n}>1) = NaN; 
    wv1 = [wv1, vek2{n}]; 
end 

wv2 = wv1(:); 
N_w = numel(wv2); 
wmax = ma(wv2); 
n_av = 500; 
n = ceil(N_w/n_av); 

wint = linspace(0, wmax, n); 
dw = wint(2); 

D_w = zeros(1,n); 

for l=1:N_w
    m = ceil(wv2(l)*(n/wmax)); 
    if m>0 && m<= n 
        D_w(m) = D_w(m)+1; 
    elseif m>n
        m = n; 
        D_w(m) = D_w(m)+1; 
    end 
end 

ws = wint*10^15; 
dw = dw*10^15; 


 
                
    
             
         