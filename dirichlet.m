function Dbc = dirichlet(phi)
[M,N] = size(phi);
phi(1,:) = 0; phi(M,:) = 0; %y boundary
phi(:,1) = 0; phi(:,N) = 0; %x boundary

BW = im2bw(-phi);
BW_filled = imfill(BW,'holes');
boundaries = bwboundaries(BW_filled);
for b = 1:length(boundaries) % how many boundaries
    for s = 1:length(boundaries{b}) % go row-wise in b-th boundary
        r = boundaries{b}(s,1); c = boundaries{b}(s,2);
        if r~=M && r~=1 && c~=1 && c~=N
            [R,C,~] = location(phi,r,c);
            phi([R,C]) = 0;
        end                      
    end    
end    

Dbc = phi; 
end   
