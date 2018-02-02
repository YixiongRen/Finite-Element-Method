function [ H,det,xbar,rhos ] = IPM( XX,R,S,nel,rho,itype)
%program
%   to evaluate the interpolation matrix hij 
%   at point(R,S) for a quadrilateral element

%---input variables--------------------------------------------------------
%   nel = number of element
%   rho = density at XX
%   itype = element type
%           eq.0 = axisymmetric
%           eq.1 = plane strain
%           eq.2 = plane stress
%   XX(2,3) = element node coordinates
%---output variables-------------------------------------------------------
%   H(2,6) = interpolation matrix H
%   det = det(J)
%   xbar = thickness in case axisymmetric
H = zeros(2,6);         % shape function hi in local coordinates (r,s)
P = zeros(2,3);         % dh/dr,ds
XJ = zeros(2,2);        % Jacobian matrix dx/dr,ds




RP = 1-R-S;
SP = R;
RM = S;


% shape function hij in local coordinates(R,S)
H(1,1) = RP;      %h11 = 1-r-s
H(1,2) = 0;               %h12 = 0
H(1,3) = SP;      %h13 = r
H(1,4) = 0;               %h14 = 0
H(1,5) = RM;      %h15 = s
H(1,6) = 0;               %h16 = 0

H(2,1) = 0;               %h21 = 0
H(2,2) = RP;      %h22 = 1-r-s
H(2,3) = 0;               %h23 = 0
H(2,4) = SP;      %h24 = r
H(2,5) = 0;               %h25 = 0
H(2,6) = RM;      %h26 = s



% natural coordinate derivatives of the shape function

%       1.with respect to R

P(1,1) = -1;
P(1,2) = 1;
P(1,3) = 0;

%       1.with respect to S

P(2,1) = -1;
P(2,2) = 0;
P(2,3) = 1;

% evaluate the Jacobian matrix at point(R,S)

for i = 1:2
    for j = 1:2
        dum = 0;
        for k = 1:3
            dum = dum + P(i,k)*XX(j,k);
        end
        XJ(i,j) = dum;
    end
end

% compute the determinant of the Jacobian matrix det(J) at point(R,S)

det = XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2);
if(le(det,0.00000001))
    error('Zero or negative Jacobian determinant for element %d\n',nel);
end


% compute the radius at point(R,S)
h = zeros(1,3);
h(1) = RP;      %h1 = 1-r-s
h(2) = SP;      %h2 = r
h(3) = RM;      %h3 = s

% density at point(R,S) from density at point XX calculation
rhos = 0;
for k =1:3
    rhos = rhos + h(k)*rho(1,k);
end



% in case of plane strain or plane stress analysis do not include the
% normal strain component
if (itype>0)
    xbar = 0;
    return;
end





xbar = 0;
for k =1:3
    xbar = xbar + h(k)*XX(1,k);
end




end


