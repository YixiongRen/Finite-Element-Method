function [ RB ] = RBody( nel, itype, nint, thic, fb, XX)
%program
%   to calculate isoparametric quadrilateral element mass matrix for
%   axisymmetric, plane stress, and plane strain conditions

%---input variables--------------------------------------------------------
%   nel = number of element
%   fb = element body force (2,4) fbx,fby at 4 corner points
%   itype = element type
%           eq.0 = axisymmetric
%           eq.1 = plane strain
%           eq.2 = plane stress
%   nint = gauss numerical intergration order
%   thic = thickness of element
%   XX(2,4) = element node coordinates
%---output variables-------------------------------------------------------
%   RB(8,1) = calculated element body force at u1v1, u2v2, u3v3, u4v4 
%
RB = zeros(8,1);            % body force for u1v1 u2v2 u3v3 u4v4... 
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
        [B,det,xbar] = STDMRI(XX,RI,SI,nel,itype);
        if (itype > 0)
            xbar = thic;
        end
        
        WT = WGT(lx,nint)*WGT(ly,nint)*xbar*det*fB;
        RB = RB + H'*WT;

    end
end



end


