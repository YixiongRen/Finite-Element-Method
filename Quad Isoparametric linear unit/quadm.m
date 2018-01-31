function [ M ] = quadm( nel, itype, nint, thic, rho, XX)
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
%   XX(2,4) = element node coordinates
%---output variables-------------------------------------------------------
%   M(8,8) = calculated element mass matrix
%
M = zeros(8,8);
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



for lx = 1:nint
    RI = XG(lx,nint);
    for ly = 1:nint
        SI = XG(ly,nint);
        [H,det,xbar,rhos] = IPM(XX,RI,SI,nel,rho,itype);
        if (itype > 0)
            xbar = thic;
        end
        
        WT = WGT(lx,nint)*WGT(ly,nint)*xbar*det*rhos;
        M = M + H'*H*WT;

    end
end



end





