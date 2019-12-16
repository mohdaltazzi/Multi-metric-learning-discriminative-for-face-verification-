function f = SPL_Req_Hard(v,lambda1)


   % -----------------> hard Reqularzation (binary)
     sum1 = 0;  
      sum1 = sum(v);
    
     f = (-lambda1) * sum1;
end
%-----------------------------------------------------------------------

 %-----------------------------------------------> %soft Reqularaztion
 %sum11 =0;
 %for i= 1:5400 
  %   xx= (v(i).^2)-(2*v(i));
   %  sum11= sum11 + xx;
% end 
 %   f = 0.5 * lambda1 * sum11;
%end
 %-----------------------------------------------> % Log Reqularaztion 
  %sum1 = 0;
 % xx = 0;
 % c = 1- lambda;
 %logc = log(c);
   % for i=1: 5400
   % xx = (c * v(i)) - (c^v(i)/log(c));
   % sum1= sum1+xx;
   % end 
%f= sum1;
     
% end 
 
%----------------------------------------------------------------------- 
 %------------------------------------------------> Mixture Reqularaztion 
 
    % lambda1 = lambda;                 % lambda= 1e-3  and lambda1> lambda2>0 
    %lambda2 =  1e-5;                  %lambda2 = 1e-5
     
    % c= (lambda1 *lambda2)/(lambda1-lambda2);
     %  for i=1:5400
      %   cur_log = log(v(i)+((1/lambda1)*c)); 
      % end 
      %  f = -c * cur_log;
%end