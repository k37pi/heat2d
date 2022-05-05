M = 51; N = 101;
D = 500; dx = 0.5; dy = dx; dt = 0.05; its = 500; 
bc = "D"; % "D" for Dirichlet, "N" for Neumann

% IC
T0 = zeros(M,N);  
x1 = 10; x2 = 30; y1 = 4; y2 = 24; Temp0 = 50;
T0((y1+1):(y2+1),(x1+1):(x2+1)) = Temp0; 

x1 = 70; x2 = 95; y1 = 30; y2 = 45; Temp0 = 70;
T0((y1+1):(y2+1),(x1+1):(x2+1)) = Temp0; 

imagesc((T0)); %imagesc(flipud(T0)); 
%xticks(0:10:(N-1)); yticks(0:10:(M-1)); 
xlabel("x"); ylabel("y");
set(gca,'YDir','normal'); axis equal;

% Scratch somewhere
s1x = 40:60; s1y = 20:30;
T0(s1y,s1x) = -1;
[mg1x,mg1y] = meshgrid(s1x,s1y);
mg1 = [reshape(mg1x,size(mg1x,2)*size(mg1x,1),1),reshape(mg1y,size(mg1x,2)*size(mg1x,1),1)]; 

%{
s2x = 20:30; s2y = 30:35;
T0(s2y,s2x) = -1;
[mg2x,mg2y] = meshgrid(s2x,s2y);
mg2 = [reshape(mg2x,size(mg2x,2)*size(mg2x,1),1),reshape(mg2y,size(mg2x,2)*size(mg2x,1),1)]; 
%}

cx = 25; cy = 32; R = 5; mg2x = []; mg2y = [];
for r = (cy-r):(cy+r)
    for c = (cx-r):(cx+r)
        if (r-cy)^2+ (c-cx)^2 <= R^2
            T0(r,c) = -1;
            mg2x = [mg2x,c]; mg2y = [mg2y,r];
        end    
    end
end    
mg2 = transpose([mg2x;mg2y]);

mg = [mg1;mg2];

imagesc((T0));
xlabel("x"); ylabel("y");
set(gca,'YDir','normal'); axis equal;

% EDIT TO INCLUDE WEIRD GRIDS

T_old = T0; 
T_new = T_old;

%Semi-implicit time 
alpha = dt*D/dx^2; beta = dt*D/dy^2; gamma = 1+2*(alpha+beta); % coefficients
idx = 0; im = {1:its}; nImages = length(im);
for n = 1:its
    %{
    T(y_range,x_range) = T(y_range,x_range) + beta*(T(y_range+1,x_range)+ T(y_range-1,x_range)) +... 
                            alpha*(T(y_range,x_range+1)+ T(y_range,x_range-1));
    T(y_range,x_range) = T(y_range,x_range)/(gamma);      
    
    T(2:(M-1),2:(N-1)) = T(2:(M-1),2:(N-1)) + alpha*(T(3:M,2:(N-1))+ T(1:(M-2),2:(N-1))) +... 
                            beta*(T(2:(M-1),3:N)+ T(2:(M-1),1:(N-2)));
    T(2:(M-1),2:(N-1)) = T(2:(M-1),2:(N-1))/(gamma);       
    %}
    
    for i = 2:(M-1)
        for j = 2:(N-1)
            if ~ismember([j,i],mg,'rows') %(ismember(i,sy) && ismember(j,sx)) 
                T_new(i,j) = T_old(i,j) + alpha*(T_old(i,j+1)+ T_old(i,j-1)) +... 
                            beta*(T_old(i+1,j)+ T_old(i-1,j));
            end
        end
    end    
    T_new = T_new/(gamma);
    
    if lower(bc) == "d"
        T_new = dirichlet(T_new); % EDIT TO INCLUDE WEIRD SHAPES
    elseif lower(bc) == "n"
        T_new = neumann(T_new); % EDIT TO INCLUDE WEIRD SHAPES
    end 
    %imagesc(T_new); 
    %xlabel("x"); ylabel("y");
    %set(gca,'YDir','normal'); 
    %sgtitle([num2str(n) ' Iterations']); axis equal;
    %drawnow; 
    
    fig = figure(1);
    idx = idx + 1;
    imagesc(T_new)
    xlabel("x"); ylabel("y");
    set(gca,'YDir','normal'); 
    sgtitle([num2str(n) ' Iterations']); axis equal;
    drawnow;
    frame = getframe(fig);
    im{idx} = frame2im(frame);
    
    T_old = T_new;
end    

filename = "dx"+dx+"_dy"+dy+"_dt"+dt+"_D"+D+"_bc"+bc+".gif"; % Specify the output file name

for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
end

