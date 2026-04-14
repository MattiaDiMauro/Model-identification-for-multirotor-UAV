function qConj = quatConj( quat )
%	Calculates the conjugate of a quaternion
%
%   qConj = quatConj( quat )
%
%   This function calculates the conjugate of a quaternion written as a
%   4-by-1 vector with the scalar component as the fourth element of the
%   vector.
%
%   References:
%   [1] Markley, F. Landis, and John L. Crassidis. Fundamentals of spacecraft 
%       attitude determination and control. Vol. 33. New York: Springer, 2014.

qv = quat(1:3);
qs = quat(4);

qConj = [ -qv  ;
           qs ];
       
end

