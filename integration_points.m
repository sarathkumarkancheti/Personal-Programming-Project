%% Integration points 
%__________________________________________________________________________
%INPUT: Integration point, number of coordinates(2 for 2D, 3 for 3D)
%__________________________________________________________________________
%OUTPUT: Local coordinate values(ξ,η,ζ) and weights(w)
%__________________________________________________________________________
function [xi,eta,zeta,w] = integration_points(i,ncoord)

   if ncoord == 2
     Xi = [-1/sqrt(3),1/sqrt(3),1/sqrt(3),-1/sqrt(3)]; %first coordinate
     Eta = [-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3)];%second coordinate
     Ci = [0,0,0,0];
     W = [1,1,1,1];  
   elseif ncoord == 3
       Xi = [-1/sqrt(3),1/sqrt(3),-1/sqrt(3),1/sqrt(3),-1/sqrt(3),1/sqrt(3),...
            -1/sqrt(3),1/sqrt(3)];
       Eta = [-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3),-1/sqrt(3),-1/sqrt(3),...
             1/sqrt(3),1/sqrt(3)];
       Ci = [-1/sqrt(3),-1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3),...
            1/sqrt(3),1/sqrt(3)];
       W = [1,1,1,1,1,1,1,1];
   end
   xi = Xi(i);
   eta = Eta(i);
   zeta = Ci(i);
   w = W(i);
end


  



  
 












