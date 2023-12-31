function L = pnumber_squares_spiral(N, option_display)
%% pnumber_squares_spiral : function to display the spiral of prime number
% squares, compute and check its total length from p_n = 5 equals L = 24k + 1.  
%
% Author, copyright and support : nicolas (dot) douillet (at) free (dot) fr, 2023.
%
%
% Inputs
%
% - N : integer scalar, N > 3. The first integer greater than the maximum wished prime number. May be or not be a prime number.
%
% - option_display : either numeric *1/0 or logical *true/false, the display option. 
%
%
% Output
%
% - L : integer scalar, the total length of the prime number squares spiral.


%% Input parsing
if nargin < 2
   
    option_display = true;
    
    if nargin  < 1
       
        N = 67;
        
    end
    
end


%% Body
u = 1:N;
v = cat(2,1,u(isprime(u))); % prime vector

M = cat(1,v(1:end-1),v(2:end),zeros(1,length(v)-1));   
M = cat(2,[0;0;0],M);
Q = cat(1,v.^2,zeros(1,length(v)),zeros(1,length(v)));

% Angle between two consecutive primes
beta = acos(M(1,1:end)./M(2,1:end)); % to build succesive rectangle triangles
beta(1,1) = 0;
theta = cumsum(beta);
 
R = @(alpha)[cos(alpha) -sin(alpha) 0;
             sin(alpha)  cos(alpha) 0;
             0           0          1];

P = Q(:,1);
k = 2;

% Compute prime squares spiral point coordinates
while k < 1 + size(M,2)
    
    P(:,k-1) = R(theta(1,k-1))*Q(:,k-1);                
    k = k+1;
    
end

P(3,:) = zeros(1,size(P,2));


% Spiral length check
L = 25 + cumsum(diff(sqrt(sum(P(:,4:end).^2,1)),1,2),2);
disp('check : mod(L-1,24) == 0');
isequal(mod(round(L)-1,24) == 0,ones(1,size(L,2)))


%% Display
if option_display
    
    figure(1);
    set(gcf,'Color',[1 1 1]);
    line(P(1,:),P(2,:),P(3,:),'Color',[0 0 1],'Linewidth',2), hold on;
    axis square, axis equal;
    grid on;
    view(2);
    gca.Cipping = 'off';
    
end


end % pnumber_squares_spiral