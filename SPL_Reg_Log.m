function f = SPL_Reg_Log( v,lambda1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%-----------------SPL_Reg_Log--------------------


 sum1 = 0;
  xx = 0;
  c = 1- lambda1;
 %logc = log(c);
    for i=1: 5400
    xx = (c * v(i)) - (c^v(i)/log(c));
    sum1= sum1+xx;
    end 
f= sum1;

end

