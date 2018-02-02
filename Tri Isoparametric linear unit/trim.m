function [ M ] = trim( nel, itype, nint, thic, rho, XX)
%program
%   to calculate isoparametric quadrilateral element mass matrix for
%   axisymmetric, plane stress, and plane strain conditions

%---input variables--------------------------------------------------------
%   nel = number of element
%   rho = density at XX
%   itype = element type
%           eq.0 = axisymmetric
%           eq.1 = plane strain
%           eq.2 = plane stress
%   nint = gauss numerical intergration order
%   thic = thickness of element
%   XX(2,3) = element node coordinates
%---output variables-------------------------------------------------------
%   M(6,6) = calculated element mass matrix
%
M = zeros(6,6);
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
 

if (nint == 1)
    for i = 1:length(XG1(:,1))
        RI = XG1(i,1);
        SI = XG1(i,2);
        [H,det,xbar,rhos] = IPM(XX,RI,SI,nel,rho,itype);
        if (itype > 0)
            xbar = thic;
        end
        
        WT = 0.5*WGT(i)*xbar*det*rhos;
        M = M + H'*H*WT;
    end
elseif (nint == 2)
    for i = 1:length(XG2(:,1))
        RI = XG2(i,1);
        SI = XG2(i,2);
        [H,det,xbar,rhos] = IPM(XX,RI,SI,nel,rho,itype);
        if (itype > 0)
            xbar = thic;
        end
        
        WT = 0.5*WGT(i)*xbar*det*rhos;
        M = M + H'*H*WT;
    end    
else
    for i = 1:length(XG3(:,1))
        RI = XG3(i,1);
        SI = XG3(i,2);
        [H,det,xbar,rhos] = IPM(XX,RI,SI,nel,rho,itype);
        if (itype > 0)
            xbar = thic;
        end
        
        WT = 0.5*WGT(i)*xbar*det*rhos;
        M = M + H'*H*WT;
    end    
end








end

















