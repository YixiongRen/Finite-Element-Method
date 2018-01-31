function [ S ] = quads( nel, itype, nint, thic, ym, pr, XX)
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
%   XX(2,4) = element node coordinates
%---output variables-------------------------------------------------------
%   S(8,8) = calculated element stiffness matrix
%
D = zeros(4,4);             % stiffness

                            % XX = zeros(2,4);
S = zeros(8,8);             % output stiffness matrix
         
DB = zeros(1,4);            % D*B, buffer for BDB

XG = [0,0,0,0;...           % gauss-legendre sampling points
        -0.5773502691896,0.5773502691896,0,0;...
        -0.7745966692415,0,0.7745966692415,0;...
        -0.8611363115941,-0.3399810435849,0.3399810435849,0.8611363115941];
XG = XG';
WGT = [2,0,0,0;...          % gauss-legendre sampling weights
        1,1,0,0;...
        0.5555555555556,0.8888888888889,0.5555555555556,0;...
        0.3478548451375,0.6521451548625,0.6521451548625,0.3478548451375];
WGT = WGT';

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
ist = 3;
if(itype == 0)
    ist = 4;
end
for lx = 1:nint
    RI = XG(lx,nint);
    for ly = 1:nint
        SI = XG(ly,nint);
        %   evaluate derivative operator B and Jacobian determinant det(J)
        [B,det,xbar] = STDM(XX,RI,SI,nel,itype);
        if (itype > 0)
            xbar = thic;
        end
        WT = WGT(lx,nint)*WGT(ly,nint)*xbar*det;
        for j = 1:8
            for k = 1:ist
                for l = 1:ist
                    DB(k) = DB(k)+ D(k,l)*B(l,j);
                end
            end
            for k = j:8
                stiff = 0;
                for l = 1:ist
                    stiff = stiff + B(l,k)*DB(l);
                end
                S(k,j) = S(k,j) + stiff*WT;
            end
        end
    end
end
for i = 1:8
    for j = i:8
        S(i,j) = S(j,i);
    end
end

 

end








