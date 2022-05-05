function [R,C,d] = location(phi,r,c) 
   
    if phi(r+1,c) >= 0
        d = 1; % up = 1
        R = r+1; C = c; 
    end    
    if phi(r-1,c) >= 0
        d = 2; % down = 2
        if ~exist("R",'var')
            R = r-1;
        else 
            R = [R;r-1];
        end   
        
        if ~exist("C",'var')
            C = c;
        else 
            C = [C;c];
        end
    end
    
    if phi(r,c+1) >= 0
        d = 4; % right = 4
        if ~exist("R",'var')
            R = r;
        else 
            R = [R;r];
        end   
        
        if ~exist("C",'var')
            C = c+1;
        else 
            C = [C;c+1];
        end 
    end
    
    if phi(r,c-1) >= 0
        d = 3; % left = 3
        if ~exist("R",'var')
            R = r;
        else 
            R = [R;r];
        end   
        
        if ~exist("C",'var')
            C = c-1;
        else 
            C = [C;c-1];
        end
    end    
        
end   
