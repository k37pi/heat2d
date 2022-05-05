function Nbc = neumann(phi)
[M,N] = size(phi);
phi(1,:) = phi(2,:); phi(M,:) = phi(M-1,:); %y boundary
phi(:,1) = phi(:,2); phi(:,N) = phi(:,N-1); %x boundary

BW = im2bw(-phi);
BW_filled = imfill(BW,'holes');
boundaries = bwboundaries(BW_filled);
for b = 1:length(boundaries) % how many boundaries
    for s = 1:length(boundaries{b}) % go row-wise in b-th boundary
        r = boundaries{b}(s,1); c = boundaries{b}(s,2);
        if r~=M && r~=1 && c~=1 && c~=N
            [R,C,d] = location(phi,r,c);
            if d(1) == 1 %up
                phi(R(1),C(1)) = phi(R(1)+1,C(1));
            elseif d(1) == 2 %down
                phi(R(1),C(1)) = phi(R(1)-1,C(1));
            elseif d(1) == 3 %left
                phi(R(1),C(1)) = phi(R(1),C(1)-1);    
            else %right
                phi(R(1),C(1)) = phi(R(1),C(1)+1);
            end
        end                      
    end    
end    

Nbc = phi; 
end   
