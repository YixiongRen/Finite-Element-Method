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
%   XX(2,4) = element node coordinates
%---output variables-------------------------------------------------------
%   H(2,8) = interpolation matrix H
%   det = det(J)
%   xbar = thickness in case axisymmetric

RP = 1+R;
SP = 1+S;
RM = 1-R;
SM = 1-S;

% shape function hij in local coordinates(R,S)
H(1,1) = 0.25*RP*SP;      %h11 = 0.25(1+R)(1+S)
H(1,2) = 0;               %h12 = 0
H(1,3) = 0.25*RM*SP;      %h13 = 0.25(1-R)(1+S)
H(1,4) = 0;               %h14 = 0
H(1,5) = 0.25*RM*SM;      %h15 = 0.25(1+R)(1+S)
H(1,6) = 0;               %h16 = 0
H(1,7) = 0.25*RP*SM;      %h17 = 0.25(1+R)(1+S)
H(1,8) = 0;               %h18 = 0
H(2,1) = 0;               %h21 = 0
H(2,2) = 0.25*RP*SP;      %h22 = 0.25(1+R)(1+S)
H(2,3) = 0;               %h23 = 0
H(2,4) = 0.25*RM*SP;      %h24 = 0.25(1+R)(1+S)
H(2,5) = 0;               %h25 = 0
H(2,6) = 0.25*RM*SM;      %h26 = 0.25(1+R)(1+S)
H(2,7) = 0;               %h27 = 0
H(2,8) = 0.25*RP*SM;      %h28 = 0.25(1+R)(1+S)


% natural coordinate derivatives of the shape function

%       1.with respect to R

P(1,1) = 0.25*SP;
P(1,2) = -1*P(1,1);
P(1,3) = -0.25*SM;
P(1,4) = -P(1,3);

%       1.with respect to S

P(2,1) = 0.25*RP;
P(2,2) = 0.25*RM;
P(2,3) = -P(2,2);
P(2,4) = -P(2,1);

% evaluate the Jacobian matrix at point(R,S)

for i = 1:2
    for j = 1:2
        dum = 0;
        for k = 1:4
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
h = zeros(1,4);
h(1) = 0.25*RP*SP;      %h1 = 0.25(1+R)(1+S)
h(2) = 0.25*RM*SP;      %h2 = 0.25(1-R)(1+S)
h(3) = 0.25*RM*SM;      %h3 = 0.25(1-R)(1-S)
h(4) = 0.25*RP*SM;      %h4 = 0.25(1+R)(1-S)
% density at point(R,S) from density at point XX calculation
rhos = 0;
for k =1:4
    rhos = rhos + h(k)*rho(1,k);
end



% in case of plane strain or plane stress analysis do not include the
% normal strain component
if (itype>0)
    xbar = 0;
    return;
end





xbar = 0;
for k =1:4
    xbar = xbar + h(k)*XX(1,k);
end




end

