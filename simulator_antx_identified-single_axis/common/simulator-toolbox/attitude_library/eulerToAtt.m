function A = eulerToAtt( euler )
%	Converts the ZYX Euler angles (RPY) to the attitude matrix A
%
%   A = eulerToAtt( euler )
%
%   The Euler angles vector must be written in radiants as a 3-by-1 vector.
%   The attitude matrix A represents the rotational matrix from the inertial
%   reference frame to the body frame of the vehicle.
%
%   References:
%   [1] Dreier, Mark E. Introduction to helicopter and tiltrotor simulation. 
%       AIAA (American Institute of Aeronautics & Ast, 2007.

phi = euler(1);
theta = euler(2);
psi = euler(3);

Ax = getAttitudeMatrix( 'x', phi );
Ay = getAttitudeMatrix( 'y', theta );
Az = getAttitudeMatrix( 'z', psi );

A = Ax * Ay * Az;

end