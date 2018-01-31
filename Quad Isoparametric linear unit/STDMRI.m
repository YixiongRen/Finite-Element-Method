function [ B,det,xbar,Tau ] = STDMRI( XX,R,S,nel,tau,itype)
%program
%   to evaluate the strain-displacement transformation matrix B
%   at point(R,S) for a quadrilateral element

%---input variables--------------------------------------------------------
%   nel = number of element
%   itype = element type
%           eq.0 = axisymmetric
%           eq.1 = plane strain
%           eq.2 = plane stress
%   tau = (3,4)or (4,4)if it's axisymmetric itype==0
%   in case of 4 strain at each XX tau(i,j) = strain i at point j
%   XX(2,4) = element node coordinates
%   S(8,8) = storage for stiffness matrix
%---output variables-------------------------------------------------------
%   B(4,8) = derivative operator B
%   det = det(J)
%   Tau = (4,1) 4 strain at gauss point(R,S)
%   xbar = thickness in case axisymmetric

H = zeros(1,4);         % shape function hi in local coordinates (r,s)
P = zeros(2,4);         % dh/dr,ds
XJ = zeros(2,2);        % Jacobian matrix dx/dr,ds
XJI = zeros(2,2);       % inverse Jacobian matrix J^-1
B = zeros(4,8);         % global derivative operator B



RP = 1+R;
SP = 1+S;
RM = 1-R;
SM = 1-S;

% shape function hi in local coordinates(r,s)
H(1) = 0.25*RP*SP;      %h1 = 0.25(1+R)(1+S)
H(2) = 0.25*RM*SP;      %h2 = 0.25(1-R)(1+S)
H(3) = 0.25*RM*SM;      %h3 = 0.25(1-R)(1-S)
H(4) = 0.25*RP*SM;      %h4 = 0.25(1+R)(1-S)


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

% compute inverse of the Jacobian matrix

dum = 1./det;
XJI(1,1) = XJ(2,2)*dum;
XJI(1,2) = -XJ(1,2)*dum;
XJI(2,1) = -XJ(2,1)*dum;
XJI(2,2) = XJ(1,1)*dum;

% evaluate global derivative operator B

k2 = 0;
for k = 1:4
    k2 = k2 + 2;
    B(1,k2-1) = 0;
    B(1,k2) = 0;
    B(2,k2-1) = 0;
    B(2,k2) = 0;    
    for i = 1:2
        B(1,k2-1) = B(1,k2-1) + XJI(1,i)*P(i,k);
        B(2,k2) = B(2,k2) + XJI(2,i)*P(i,k);
    end
    B(3,k2) = B(1,k2-1);
    B(3,k2-1) = B(2,k2);
end


% in case of plane strain or plane stress analysis do not include the
% normal strain component
if (itype>0)
    xbar = 0;
    Tau = zeros(4,1);
    for i = 1:3
        for j = 1:4
            Tau(i,1) = Tau(i,1) + H(j)*tau(i,j);
        end
    end
    return;
end

% compute the radius at point(R,S)

% Initial strain vector at point(R,S)

Tau = zeros(4,1);
for i = 1:4
    for j = 1:4
        Tau(i,1) = Tau(i,1) + H(j)*tau(i,j);
    end
end



xbar = 0;
for k =1:4
    xbar = xbar + H(k)*XX(1,k);
end

% evaluate the hoop strain-displacement relation

if(le(xbar,0.00000001))
    % for the case of zero radius equate radial to hoop strain
    for k = 1:8
        B(4,k) = B(1,k);
    end
    return;
else
    % non-zero radius
    dum = 1./xbar;
    k2=0;
    for k = 1:4
        k2 = k2+2;
        B(4,k2) = 0;
        B(4,k2-1) = H(k)*dum;
    end
    return;
end

end



