function euler = attToEuler( A )
%	Converts the attitude matrix A to the ZYX Euler angles (RPY)
%
%   euler = attToEuler( A )
%
%   The attitude matrix A represents the rotational matrix from the inertial
%   reference frame to the body frame of the vehicle.
%   The returned Euler angles vector will be written in radiants as a 3-by-1 vector.
%
%   References:
%   [1] Dreier, Mark E. Introduction to helicopter and tiltrotor simulation. 
%       AIAA (American Institute of Aeronautics & Ast, 2007.

phi = atan2( A(2,3), A(3,3) );
theta = - asin( A(1,3) );
psi = atan2( A(1,2), A(1,1) );

euler = [ phi, theta, psi ]';

end