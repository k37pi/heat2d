close all; clear vars-global; 

% https://hplgit.github.io/fdm-book/doc/pub/diffu/pdf/diffu-4print.pdf

alld = 1:4; % directions up(1), down(2), left(3), right(4)


M = 51; N = 101;
D = 500000; dx = 0.5; dy = dx; dt = 0.00005; its = 1000; 
idx = 1; im = {1:its}; nImages = length(1:its); % for gif
bc = "d"; % "D" for Dirichlet, "N" for Neumann

% IC
T0 = zeros(M,N);   
x1 = 10; x2 = 30; y1 = 4; y2 = 24; Temp0 = 50;
T0((y1+1):(y2+1),(x1+1):(x2+1)) = Temp0; % source 1

x1 = 70; x2 = 95; y1 = 30; y2 = 45; Temp0 = 70;
T0((y1+1):(y2+1),(x1+1):(x2+1)) = Temp0; % source 2

x1 = 55; x2 = 65; y1 = 1; y2 = 35; Temp0 = 100;
T0((y1+1):(y2+1),(x1+1):(x2+1)) = Temp0; % source 3

nsources = 3;

imagesc((T0)); %imagesc(flipud(T0)); 
%xticks(0:10:(N-1)); yticks(0:10:(M-1)); 
xlabel("x"); ylabel("y");
set(gca,'YDir','normal'); %axis equal;

% Scratch somewhere
s1x = 40:60; s1y = 20:30; % rectangular scratch 1
T0(s1y,s1x) = -1;
[mg1x,mg1y] = meshgrid(s1x,s1y);
mg1 = [reshape(mg1x,size(mg1x,2)*size(mg1x,1),1),reshape(mg1y,size(mg1x,2)*size(mg1x,1),1)]; 

%{
s2x = 20:30; s2y = 30:35; % rectangular scratch 2
T0(s2y,s2x) = -1;
[mg2x,mg2y] = meshgrid(s2x,s2y);
mg2 = [reshape(mg2x,size(mg2x,2)*size(mg2x,1),1),reshape(mg2y,size(mg2x,2)*size(mg2x,1),1)]; 
%}

cx = 25; cy = 32; R = 5; mg2x = []; mg2y = []; % circular scratch 1
for r = (cy-R):(cy+R)
    for c = (cx-R):(cx+R)
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

points_row = [10;19;38;35;20]; 
points_col = [18;59;83;14;80]; 
legend_text = strings(length(points_row)+1,1);
for i = 1:length(points_row)
    legend_text(i) = "p"+i;
end    
legend_text(length(points_row)+1) = "Mean T";

temp_pt = zeros(its+1,length(points_row));
mean_temp = zeros(its+1,1);

temp_pt(1,:) = diag(T0(points_row,points_col));
mean_temp(1) = mean(mean(T0(T0>=0)));

%co = colororder();
%nco = co(2:7,:);

fig = figure(1);
subplot(2,1,1);imagesc(T0); hold on
for i = 1:length(points_row)
    plot(points_col(i),points_row(i),'*','MarkerSize',7,'linewidth',1.6); hold on
end    
hold off; %colororder(nco);
xlabel("x"); ylabel("y");
set(gca,'YDir','normal'); 
title('Initial'); colorbar; %axis equal; 
drawnow;
subplot(2,1,2);
for i = 1:length(points_row)
    plot(temp_pt(:,i)); hold on
end
plot(mean_temp,'r'); hold off; xlabel("time"); ylabel("T"); xlim([0 its])    
lh=legend(legend_text); %lh = findobj(gcf, 'Type', 'legend');
neworder = [2:(length(points_row)+1),1];
lh.PlotChildren = lh.PlotChildren(neworder); 
%lh.Position(1) = 0.8;
%lh.Position(2) = 0.4;
frame = getframe(fig); % collect frames
im{idx} = frame2im(frame);

%% Semi-implicit time --
alpha = dt*D/dx^2; beta = dt*D/dy^2; gamma = 1+2*(alpha+beta); % coefficients

for n = 1:its
    for i = 2:(M-1)
        for j = 2:(N-1)
            if ~ismember([j,i],mg,'rows') 
                [nr,nc,d,nn] = location(T_old,i,j); 
                if nn == 4
                   %(ismember(i,sy) && ismember(j,sx)) 
                        T_new(i,j) = T_old(i,j) + alpha*(T_old(i,j+1)+ T_old(i,j-1)) +... 
                                    beta*(T_old(i+1,j)+ T_old(i-1,j));
                else
                    missing_dir = alld(~ismember(alld,d));
                    T_new(i,j) = T_old(i,j) + alpha*(T_old(i,j+1)+ T_old(i,j-1)) +... 
                                    beta*(T_old(i+1,j)+ T_old(i-1,j));
                    if ismember(1,missing_dir) %up
                        T_new(i,j) =  T_new(i,j)-beta*T_old(i+1,j);
                    end
                    if ismember(2,missing_dir) %down
                        T_new(i,j) =  T_new(i,j)-beta*T_old(i-1,j);
                    end
                    if ismember(3,missing_dir) %left
                        T_new(i,j) =  T_new(i,j)-alpha*T_old(i,j-1);   
                    end
                    if ismember(4,missing_dir) %right
                        T_new(i,j) =  T_new(i,j)-alpha*T_old(i,j+1);   
                    end          
                end
            end
        end
    end    
    T_new = T_new/gamma;
    
    if lower(bc) == "d"
        T_new = dirichlet(T_new); % EDIT TO INCLUDE WEIRD SHAPES
    elseif lower(bc) == "n"
        T_new = neumann(T_new); % EDIT TO INCLUDE WEIRD SHAPES
    end
    
    temp_pt(n+1,:) = diag(T_new(points_row,points_col));
    mean_temp(n+1) = mean(mean(T_new(T_new>=0)));

    fig = figure(1);
    idx = idx + 1;
    subplot(2,1,1);imagesc(T_new); hold on;
    for i = 1:length(points_row)
        plot(points_col(i),points_row(i),'*','MarkerSize',7,'linewidth',1.6); hold on
    end    
    hold off; %colororder(nco);
    xlabel("x"); ylabel("y");
    set(gca,'YDir','normal'); 
    title([num2str(n) ' Iterations']); colorbar; %axis equal;
    drawnow;
    subplot(2,1,2);
    for i = 1:length(points_row)
        plot(temp_pt(:,i)); hold on
    end     
    plot(mean_temp,"r"); hold on; xlabel("time"); ylabel("T"); xlim([0 its]);
    hold off;
    lh=legend(legend_text); %lh = findobj(gcf, 'Type', 'legend');
    neworder = [2:(length(points_row)+1),1];
    lh.PlotChildren = lh.PlotChildren(neworder);
    frame = getframe(fig); % collect frames
    im{idx} = frame2im(frame);
     
    T_old = T_new;
end    
%% End --

%% Saving gif --
%{
filename = "dx"+dx+"_dy"+dy+"_dt"+dt+"_D"+D+"_bc"+bc+"_nsrc"+nsources+".gif"; % Specify the output file name

for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
    elseif mod(idx,10)==0
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
end
%}
%% End ---

