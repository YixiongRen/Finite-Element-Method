function [ S ] = tris( nel, itype, nint, thic, ym, pr, XX)
%program
%   to calculate isoparametric quadrilateral element stiffness matrix for
%   axisymmetric, plane stress, and plane strain conditions

%---input variables--------------------------------------------------------
%   nel = number of element
%   itype = element type
%           eq.0 = axisymmetric
%           eq.1 = plane strain
%           eq.2 = plane stress
%   nint = gauss numerical intergration order
%   thic = thickness of element
%   ym = young's modulus
%   pr = poisson's ratio
%   XX(2,3) = element node coordinates
%---output variables-------------------------------------------------------
%   S(6,6) = calculated element stiffness matrix
%
D = zeros(4,4);             % stiffness

                            % XX = zeros(2,4);
S = zeros(6,6);             % output stiffness matrix
         

XG1 = [0.1666666666667,0.1666666666667;...
        0.6666666666667,0.1666666666667;...
        0.1666666666667,0.6666666666667];                    % gauss-legendre sampling points when nint = 3
    
    
XG2 = [0.1012865073235,0.1012865073235;...
        0.7974269853531,0.1012865073235;...
        0.1012865073235,0.7974269853531;...
        0.4701420641051,0.0597158717898;...
        0.4701420641051,0.4701420641051;...
        0.0597158717898,0.4701420641051;...
        0.3333333333333,0.3333333333333];                       % gauss-legendre sampling points when nint = 7
    
    
XG3 = [0.0651301029022,0.0651301029022;...
        0.8697397941956,0.0651301029022;...
        0.0651301029022,0.8697397941956;...
        0.3128654960049,0.0486903154253;...
        0.6384441885698,0.3128654960049;...
        0.0486903154253,0.6384441885698;...
        0.6384441885698,0.0486903154253;...
        0.3128654960049,0.6384441885698;...
        0.0486903154253,0.3128654960049;...
        0.2603459660790,0.2603459660790;...
        0.4793080678419,0.2603459660790;...
        0.2603459660790,0.4793080678419;...
        0.3333333333333,0.3333333333333];                       % gauss-legendre sampling points when nint = 13

WGT1 = [0.3333333333333;0.3333333333333;0.3333333333333];       % gauss-legendre sampling weights when nint = 3

WGT2 = [0.1259391805448;...
        0.1259391805448;...
        0.1259391805448;...
        0.1323941527885;...
        0.1323941527885;...
        0.1323941527885;...
        0.225];                                                 % gauss-legendre sampling weights when nint = 7

WGT3 = [0.0533472356088;...
        0.0533472356088;...
        0.0533472356088;...
        0.0771137608903;...
        0.0771137608903;...
        0.0771137608903;...
        0.0771137608903;...
        0.0771137608903;...
        0.0771137608903;...
        0.1756152574332;...
        0.1756152574332;...
        0.1756152574332;...
        -0.1495700444677];                                      % gauss-legendre sampling weights when nint = 13

%  Stress-strain law
F = ym/(1. +pr);
G = F*pr/(1.-2.*pr);
H = F+G;
%  plane strain analysis
D(1,1) = H;
D(1,2) = G;
D(1,3) = 0;
D(2,1) = G;
D(2,2) = H;
D(2,3) = 0;
D(3,1) = 0;
D(3,2) = 0;
D(3,3) = F/2;

if(itype==1)
    thic = 1;
elseif(itype==0)
    %   axisymmetric analysis
    D(1,4) = G;
    D(2,4) = G;
    D(3,4) = 0;
    D(4,1) = G;
    D(4,2) = G;
    D(4,3) = 0;
    D(4,4) = H;
else
    %   for plane stress analysis condense stress-strain matrix
    for i = 1:3
        A = D(i,4)/D(4,4);
        for j = i:3
            D(i,j) = D(i,j) - D(4,j)*A;
            D(j,i) = D(i,j);
        end
    end
end
% calculate element stiffness
if (nint == 1)
    for i = 1:length(XG1(:,1))
        RI = XG1(i,1);
        SI = XG1(i,2);
        [B,det,xbar] = STDM(XX,RI,SI,nel,itype);
        if (itype > 0)
            xbar = thic;
        end
        WT = 0.5*WGT1(i)*xbar*det;
        S = S + B'*D*B*WT;
    end
elseif (nint == 2)
    for i = 1:length(XG2(:,1))
        RI = XG2(i,1);
        SI = XG2(i,2);
        [B,det,xbar] = STDM(XX,RI,SI,nel,itype);
        if (itype > 0)
            xbar = thic;
        end
        WT = 0.5*WGT2(i)*xbar*det;
        S = S + B'*D*B*WT;
    end    
else
    for i = 1:length(XG3(:,1))
        RI = XG3(i,1);
        SI = XG3(i,2);
        [B,det,xbar] = STDM(XX,RI,SI,nel,itype);
        if (itype > 0)
            xbar = thic;
        end
        WT = 0.5*WGT3(i)*xbar*det;
        S = S + B'*D*B*WT;
    end    
end







end



