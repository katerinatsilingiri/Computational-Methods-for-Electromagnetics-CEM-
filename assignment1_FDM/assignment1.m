% Katerina Tsilingiri, 2806 %
% Program that implements the Finite Difference Method with 2 ways  %
% and the analytical method for comparison.                         %
% The limits in the x and y axis are predefined as 4 and the user   %
% can choose the step h.                                            %

disp("Current limits for the x-axis and y-axis: 4");
disp("Step h should be be in range (0, 1]");
max_x = 4;
max_y = 4;

% Define the step h %
prompt = "Step h:";
h = input(prompt);

% Discretize x & y based on x_max, y_max & h %
xx = h:h:max_x-h;
yy = h:h:max_y-h;

% Define the dimension of sparse A %
dimension = (max_x/h) + 1;
A = zeros(dimension,dimension); 

% Calculate the number of points in the x & y axis %
pointsx = max_x/h;
pointsy = max_y/h;

% Set the boundary conditions %
A(1, :) = 0; A(:, 1) = 0; A(pointsx+1, :) = 0;
A(:, pointsy+1) = 100;

%           ------------ Computational Method: FDM ------------        %
% 1st implementation                                                   %
% This method solves the given equations using the elements of the the %
% sparse array A. For the method to converge, it must be repeated      %
% because initially all of the values of the sparse array A are zero.  %
for k=1:100
    for i=(dimension-1):-1:2
        for j=(dimension-1):-1:2
            A(i,j) = (A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i, j-1))/4;
            %fprintf('A(%d, %d) = %d\n', i, j, A(i,j));
        end
    end
end

% Separate the boundary conditions in order to plot them %
A_1 = A(2:pointsx,2:pointsx); % A_1 is the Phi except the boundary conditions % 

% 2nd implementation %
% This method constructs the sparse array A using the diagonals of the %
% matrix. This method is parametrized based on the number of points in %
% the x-axis and the number of points in the y-axis.                   %
mainDiag = -4 * ones((pointsx-1)*(pointsy-1), 1); % create the main diagonal (with -4) %

onesDiag = ones((pointsx-1)*(pointsy-1), 1); % create the 3rd diagonal (with 1) %

oneZerosDiag = ones(((pointsx-1)*(pointsy-1))-1, 1); % create the 1st diagonal (with 1) %
                                                     % this diagonal should have 1's    %
                                                     % except the positions i*(number   %
                                                     % of points in x)                  %   

positionsofZeros = pointsx-1:pointsx-1:(pointsx-1)*(pointsy-1); % define the positions of the zeros %
                                                                % multiple of the number of points in x-axis %

% Create the diagonal vectors for the -1 and +1 diagons %
oneZerosDiag(positionsofZeros) = 0;
oneZerosDiag2 = ones(((pointsx-1)*(pointsy))-1, 1);
positionsofZeros = pointsx-1+1:pointsx-1:(pointsx-1)*(pointsy-1);
oneZerosDiag2(positionsofZeros) = 0;

A_2 = spdiags(mainDiag, 0, (pointsx-1)*(pointsy-1), (pointsx-1)*(pointsy-1));
A_2 = spdiags(onesDiag, pointsx-1, A_2);
A_2 = spdiags(onesDiag, -(pointsx - 1), A_2);
A_2 = spdiags(oneZerosDiag2, 1, A_2);
A_2 = spdiags(oneZerosDiag, -1, A_2);

% Create the b vector for the system Ax=b         %
% The positions of the values are chosen based on %
% the number of points in x and y axis.           %
b_init = zeros(pointsx-1,pointsy-1); 
b_init(1, :) = 0; b_init(:, 1) = 0; b_init(pointsx-1, :) = 0;
b_init(:, pointsy-1) = 100;
b_vec = reshape(-b_init, [(pointsx-1)*(pointsx-1), 1]);


% Solve the sparse system                               %
% x_sol_array is the Phi except the boundary conditions % 
x_sol = A_2\b_vec;
x_sol_array = reshape(x_sol, [(pointsx-1), (pointsx-1)]); 

%          ------------ Analytical Method ------------         %
% This method was mainly implemented to compare the results of %
% of the FDM Methods implemented above.                        %
% Define the maximum number of iterations                      % 
N = 5;
A_new = zeros(dimension,dimension);  % A_new is the matrix used in the analytical method %

% Define the boundary conditions %
A_new(1, :) = 0; A_new(:, 1) = 0; A_new(pointsx+1, :) = 0;
A_new(:, pointsy+1) = 100;

% x_point & y_point define the x,y values in the following equation %
x_point = h;
y_point = 0;

for i = 2:1:pointsx
    y_point = 0;
    for j = 2:1:pointsy
        y_point = y_point + h;
        for n = 1:2:N
            A_new(i,j) = A_new(i,j) + (( 4 * 100 * sin(n * pi * (x_point) / 4) * sinh(n * pi * (y_point) / 4) ) / (n * pi * sinh(n * pi * 4 / 4)));
        end 
    end
    x_point = x_point + h;
end


A_new2 = A_new(2:pointsx,2:pointsx); % A_new2 is the matrix used in the analytical method %

% Plot the 3 solutions for comparison %
figure(1)
[X,Y] = meshgrid(xx, yy);
tiledlayout(3,1)

% 1st plot: FDM (1st implementation) %
ax1 = nexttile;
surf(X,Y,A_1);
title(ax1,'Computational Solution using FDM (1st implementation)')
grid(ax1,'on')

% 2nd plot: FDM (2nd implementation) %
ax2 = nexttile;
surf(X,Y, x_sol_array);
title(ax2,'Computational Solution using FDM (2nd implementation)')
grid(ax2,'on')

% 3rd plot: Analytical implementation %
ax3 = nexttile;
surf(X,Y,A_new2);
title(ax3,'Analytical Solution')
grid(ax2,'on')

% Calculate the errors of each method %
errorFDM1_An = A_1 - A_new2;
errorFDM2_An = x_sol_array - A_new2;
errorFDM = A_1 - x_sol_array;

fprintf('Error (FDM1-Analytical):\n');
disp(errorFDM1_An);
fprintf('Error (FDM2-Analytical):\n');
disp(errorFDM2_An);
fprintf('Error (FDM1-FDM2):\n');
disp(errorFDM);
