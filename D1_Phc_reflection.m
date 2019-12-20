%MATLAB Program Reflectance Spectrum Computation of 1D PhC
%The program computes the reflectance spectrum
%of the user-defined layered structure by the method
%described in this chapter. The solution of the
%linear equations system is carried out by Cramer?s
%method. For this reason, two matrices are composed.
%The first one contains coefficients at the unknowns
%of the system while the second one contains free
%terms.
%Input parameters: layers thicknesses and refractive
%indices, the range of the wavelengths for the
%spectrum computation.

%Output data: The plot of the layered structure
   %reflectance spectrum within required range
   %-------------------------------------------------
   clear all;
   %The array lambda defines the set of the wavelength
   %where the reflectance computation will be carried
%out.
lambda=1:0.01:2;
%The array layerWidth defines the thickness of each
%layer. The number of elements in this array defines
%the number of layers in the structure.
layerWidth= [0.15,0.25, 0.15]; 
layerWidth=[0.15 0.25 0.15 0.25 0.15...
                              0.25 0.15 0.25 0.15];
%The array layerRI defines the refractive indices of
%layers. I should have the same length as layerWidth.
%Here, it is assumed that surrounding material
%refractive index equala 1 i.e. the structure is
%surrounded by air.
layerRI =[1, 3.5, 1]; 
layerRI=[3.5 1 3.5 1 3.5 1 3.5 1 3.5];

%The cycle over wavelength. The variable wavelength
%is counter and it takes on values from 1 to length
%of the array lambda
for wavelength=1:length(lambda)
%Here we compute the wave number
  k=2*pi/lambda(wavelength);
%Clearing the arrays from the values written at
%previous cycle step.
  clear coefMatrix;
  clear freeMemberMatrix;
%Two-dimensional array coefMatrix contains the
%coefficients at unknowns of the
%system of linear equations. Because the coefficient
%A0 is user-defined (as we determined from initial
%conditions) first two rows of the matrix contains
%coefficients B0, A1, B1 only. Consequently, they
%should be defined outside the cycle. Moreover, the
%coordinate of the first interface is assumed to be
%equal to zero so all the exponents here vanish.
coefMatrix(1,1:3)=[1 -1 -1];
coefMatrix(2,1:3)=[-1i*k*1 -1i*k*layerRI(1) ...
                               1i*k*layerRI(1)];
%The variable xCurrent contains the coordinate of
%current interface i.e. the interface we currently
%writing boundary conditions for
  xCurrent=layerWidth(1);
%The following cycle scans over all interfaces. The
%scanning begins from the interface between two first
%layers
  for countBoundary=1:length(layerWidth)-1
%Defining the function equality condition at the
%boundary.
    coefMatrix(countBoundary*2+1,2*countBoundary:...
2*countBoundary+3)=...
     [exp(1i*k*layerRI(countBoundary)*xCurrent)...
      exp(-1i*k*layerRI(countBoundary)*xCurrent)...
     -exp(1i*k*layerRI(countBoundary+1)*xCurrent)...
     -exp(-1i*k*layerRI(countBoundary+1)*xCurrent)];
%Defining the function derivatives equality condition
%at the boundary.
  coefMatrix(countBoundary*2+2,2*countBoundary:...
2*countBoundary+3)=...
        [1i*k*layerRI(countBoundary)*...
exp(1i*k*layerRI(countBoundary)*xCurrent)...
         -1i*k*layerRI(countBoundary)*...
exp(-1i*k*layerRI(countBoundary)*xCurrent)...
        -1i*k*layerRI(countBoundary+1)*...
exp(1i*k*layerRI(countBoundary+1)*xCurrent)...
        1i*k*layerRI(countBoundary+1)*...
exp(-1i*k*layerRI(countBoundary+1)*xCurrent)];
%The coordinate of current interface is increased by
%the thickness of current layer. After that it
%contains the coordinate of the next interface
    xCurrent=xCurrent+layerWidth(countBoundary+1);
  end

%The boundary conditions at the last interface are
%defined manually as it was done for the first one
%because we constituted the coefficient BN+1 to be
%equal to zero.
    coefMatrix(length(layerWidth)*2+1,2*...
length(layerWidth):2*length(layerWidth)+2)=...
   [exp(1i*k*layerRI(length(layerRI))*xCurrent)...
    exp(-1i*k*layerRI(length(layerRI))*xCurrent)...
   -exp(1i*k*1*xCurrent)];
coefMatrix(length(layerWidth)*2+2,2*...
   length(layerWidth):2*length(layerWidth)+2)=...
        [1i*k*layerRI(length(layerRI))*...
   exp(1i*k*layerRI(length(layerRI))*xCurrent)...
        -1i*k*layerRI(length(layerRI))*...
   exp(-1i*k*layerRI(length(layerRI))*xCurrent)...
        -1i*k*1*exp(1i*k*1*xCurrent)];
   %Defining column matrix containing the array of free
   %terms.
   %Nonzero elements of this array are only two first
   %elements because we assumed A0 to be equal to 1. For
   %all other equations free terms equal to zero
     freeMemberMatrix=zeros(length(layerWidth)*2+2,1);
     freeMemberMatrix(1)=-1;
     freeMemberMatrix(2)=-1*1i*k*1;
     %Cycle for Cramer?s method solution
   %Cashing the matrix to some temporary variable.
       tempMatrix=coefMatrix;
   %Substituting free-terms column matrix instead of
   %column with number corresponding the number of
   %solution we are searching for.
       tempMatrix(:,1)=freeMemberMatrix;
   %Searching for the solution by dividing the
   %determinant of temporary matrix by the
   %determinant of the matrix composed of coefficients
   %at the unknowns
       solution(wavelength)=...
                    1/det(coefMatrix)*det(tempMatrix);
end
   %Plotting the solution
   plot(lambda,abs(solution).^2,'LineWidth',2);
   xlabel('\lambda, \mum','FontSize',16);
   ylabel('Refrectance, r.u.','FontSize',16);