function f = SPL_Reg_Soft_linear( v,lambda1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% -----------------Reg_Soft_linear-------------------


sum11 =0;
 for i= 1:5400 
     xx= (v(i).^2)-(2*v(i));
     sum11= sum11 + xx;
 end 
    f = 0.5 * lambda1 * sum11;
end

