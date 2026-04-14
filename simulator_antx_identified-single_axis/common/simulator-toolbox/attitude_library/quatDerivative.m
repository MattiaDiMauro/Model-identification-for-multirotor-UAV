function qDot = quatDerivative( ome, quat )
%   Calculates the quaternion derivative given the angular velocity
%
%   qDot = quatDerivative( ome, quat )
%
%   The angular velocity (ome) expressed in body frame must be written in 
%   radiants per seconds as a 3-by-1 vector.
%   The quaternion (quat) must be written as a 4-by-1 vector with the scalar
%   element as the fourth component.
%
%   References:
%	[1] Markley, F. Landis. "Attitude error representations for Kalman filtering." 
%       Journal of guidance control and dynamics 26.2 (2003): 311-317.

qDot = 0.5 * quatProd( [ome, 0]', quat);

end

