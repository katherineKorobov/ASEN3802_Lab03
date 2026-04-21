function C_D=Total_Drag(C_di,c_d)
%Inputs: sectional drag coef., induced drag coef.
%outputs: Total drag coef
%Author: John Heflin and Kiana Watson

mat = load('NACA_0012_cd.mat');
cd_0012 = mat.sorted_data; 

mat = load('NACA_2412_cd.mat');
cd_2412 = mat.sorted_data; 

C_Dm= C_di + c_d



end
