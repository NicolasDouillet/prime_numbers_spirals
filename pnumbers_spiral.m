function S = pnumbers_spiral(N, option_display)
%% pnumbers_spiral : function to display the spiral of prime numbers laying
% on 6n +/- 1 circles, and check its sum.
%
% Author : nicolas.douillet (at) free.fr, 2023-2025.
%
%
% Inputs
%
% - N : integer scalar, N > 3. The first integer greater than the maximum prime wished number. May be or not be a prime number.
%
% - option_display : either numeric *1/0 or logical *true/false, the display option. 
%
%
% Output
%
% - S : integer scalar, equals to 1 + the sum of primes less or equal to N.


% TODO :
%
% + 2x2 squares  in rectangle triangles


%% Input parsing
if nargin < 2
   
    option_display = true;
    
    if nargin  < 1
       
        N = 67;
        
    end
    
end


%% Body
v = cat(2,1,primes(N)); % prime vector

M = cat(1,v(1:end-1),v(2:end),zeros(1,length(v)-1));
M = cat(2,[0;0;0],M);
Q = cat(1,v(1:end-1),zeros(1,length(v)-1),zeros(1,length(v)-1));

% Angle between two consecutive primes
beta = acos(M(1,1:end)./M(2,1:end)); % to build succesive rectangle triangles
beta(1,1) = 0;
theta = cumsum(beta);

R = @(alpha)[cos(alpha) -sin(alpha) 0;
    sin(alpha)  cos(alpha) 0;
    0           0          1];

P = Q(:,1);
k = 2;

% Compute prime spiral point coordinates
while k < 1 + size(M,2)
    
    P(:,k-1) = R(theta(1,k-1))*Q(:,k-1);
    k = k+1;
    
end

P(3,:) = v(1:end-1);


% Check expression of primes sum
S = sum(v);
disp(['check : mod(primes sum,6) < ',num2str(length(v)-1)]);
isequal(mod(S,6) < length(v)-1,true) % length(v)-1 since 1 not prime doesn't count (in prime indices)


%% Display
if option_display
    
    figure(1);
    set(gcf,'Color',[0 0 0]);
    
    plot3(P(1,:),P(2,:),P(3,:),'o','Color',[1 1 1],'Linewidth',3), hold on;
    line(P(1,:),P(2,:),P(3,:),'Color',[0.875 0 0],'Linewidth',2), hold on;
    
    % Rectangle triangles and Z axis display option
    for k = 1:size(P,2)
    
       line([0,P(1,k)],[0,P(2,k)],[P(3,k),P(3,k)],'Color',[1 1 1],'Linewidth',2), hold on; % [0 0 0.875]
    
    end
   
    line([0 0],[0 0],[v(1,1) v(1,end)],'Color',[1 1 1],'Linewidth',2), hold on;
    
    axis square, axis equal;
    grid on;
    view(2);
    % gca.Cipping = 'off';
    
    
    % Plot polar centered circles of prime radii    
    cmap = colormap(jet);    
    alpha = linspace(0,2*pi,361);
    base_circle_X = cos(alpha);
    base_circle_Y = sin(alpha);
    
    for k = v(4:end)
        
        line(k*base_circle_X,k*base_circle_Y,k*ones(1,numel(alpha)),'Color',cmap(round(k*size(jet,1)/N),:)), hold on;
        
    end
    
    xlabel('Signed integer unit','FontSize',16)
    ylabel('Signed integer unit','FontSize',16)
    title({['Spiral of prime numbers from 5 to ',num2str(N)],'laying on 6n \pm 1 radii / circles'},'FontSize',16);
    legend('Position of primes $$p = 6n \pm 1$$ on the circles','Distance $$l = \sqrt{p''^2-p^2}$$ between consecutives primes',...
           'Radii / values of primes','Interpreter','latex','FontSize',12,'Color',[0 0 0],'EdgeColor',[1 1 1],'TextColor',[1 1 1],'FontSize',12);

%     xlabel('Rayon (unit� enti�re)','FontSize',16);
%     ylabel('Rayon (unit� enti�re)','FontSize',16);
%     title({['Spirale des nombres premiers de 5 � ',num2str(N)],'bas�e sur les cercles de rayons 6n \pm 1'},'FontSize',16); % + vue polaire   
%     legend('Position des premiers $$p = 6n \pm 1$$ sur les cercles','Distance $$l = \sqrt{p''^2-p^2}$$ entre les premiers consecutifs',...
%            'Rayons / valeurs des premiers','Interpreter','latex','FontSize',12,'Color',[0 0 0],'EdgeColor',[1 1 1],'TextColor',[1 1 1]);
    
    box on;
    set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
    
end


end % pnumbers_spiral