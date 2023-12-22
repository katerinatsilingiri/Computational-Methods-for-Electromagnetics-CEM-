% Katerina Tsilingiri, 2806 %
% Program that implements the Method of Moments for the given input %
% file. The format of the input file is the same as FastCap:        %
% Ca                                                                %
% Qy1 x1 y1 x2 y2 x3 y3 x4 y4                                       %
% Qy2 ... ... ..                                                    %
%                                                                   %    
% Cb                                                                %
% Qy3 x1 y1 x2 y2 x3 y3 x4 y4                                       %
% where (x1,y1): lower left corner of quarter                       %
%       (x2, y2): upper left of corner of quarter                   %
%       (x3, y3): upper right of corner of quarter                  %
%       (x4, y4): lower right of corner of quarter                  %

fID = fopen('input.txt', 'r');

conductor_num = 0;     % indicates the number of conductor %

conductor_array = [];  % array that contains the coordinates of all %
                       % quarters of all the present conductors     %

new_row = zeros(1, 9); % new empty entry for the array                    %
                       % format of each entry :                           %
                       % (index of current conductor) (x1) (y1) (x2) (y2) %
                       %                              (x3) (y3) (x4) (y4) %

% Parse the input file & create the conductor array %                      
while ~feof(fID)
    line = fgetl(fID);
    data = strsplit(line,' ');
    
    % locate the conductor in the input file and increase the counter %
    if contains(data, 'C')
        conductor_num = conductor_num +1;
        new_row(:) = 0;
    end
    
    % locate the quarters of each conductor, parse them & insert them %
    % into the conductor_array                                        %
    if contains(data(1), 'Q')
        new_row(1) = conductor_num;
        new_row(2) = str2double(data(2));
        new_row(3) = str2double(data(3));
        new_row(4) = str2double(data(4));
        new_row(5) = str2double(data(5));
        new_row(6) = str2double(data(6));
        new_row(7) = str2double(data(7));
        new_row(8) = str2double(data(8));
        new_row(9) = str2double(data(9));
        conductor_array = [conductor_array; new_row];
    end
end

% Display the (initial) conductor array %
fprintf("Conductor array\n");
fprintf("    id    x1    y1    x2    y2    x3    y3    x4    y4\n");
disp(conductor_array);

quarters = height(conductor_array); % total number of quarters %

A1 = zeros(quarters, quarters); % array of the system Ax=b %

curr_center_x = 0;
curr_center_y = 0;
relative_center_x = 0;
relative_center_y = 0;
hx = 0;
hy = 0;

% For each quarter of all conductors in the conductor array : %
%       - calculate the center coordinates of the current     %
%         quarter                                             %
%       - for each quarter in the conductor array :           %
%           - get the height & width                          %
%           - calculate A[i][j]                               %

% Formulas to use in order to calculate the (x,y) of the center %
% and the surface of each quarter                               %    
%         center_x = (x4 - x1)/2 + x1                           %
%         center_x = (row(8) - row(2))/2 + row(2);              %
%         center_y = (y2 - y1)/2 + y1                           %  
%         center_y = (row(5) - row(3))/2 + row(3);              %
%         surface = (x4 - x1) * (y2 - y1)                       %
%         surface = ((row(8) - row(2))) * ((row(5) - row(3)));  %
disp("Iteration : 0");
fprintf("  Construction of A (%d x %d)\n", quarters, quarters);
for i = 1:quarters
    curr_center_x = (conductor_array(i, 8) - conductor_array(i, 2))/2 + conductor_array(i, 2);
    curr_center_y = (conductor_array(i, 5) - conductor_array(i, 3))/2 + conductor_array(i, 3);
    for j = 1:quarters
        hx = (conductor_array(j, 8) - conductor_array(j, 2));
        hy = (conductor_array(j, 5) - conductor_array(j, 3));
        if (i == j)
            A1(i, j) = (2*(hx * log((hy + sqrt(hx^2 + hy^2))/hx) + hy * log((hx + sqrt(hx^2 + hy^2))/hy))) / (4*pi*8.85*10^(-12));
        else
            relative_center_x = (conductor_array(j, 8) - conductor_array(j, 2))/2 + conductor_array(j, 2);
            relative_center_y = (conductor_array(j, 5) - conductor_array(j, 3))/2 + conductor_array(j, 3);
            A1(i, j) = ((hx * hy)/sqrt(((curr_center_x - relative_center_x)^2 + (curr_center_y - relative_center_y)^2))) / (4*pi*8.85*10^(-12));
        end 
    end
end

% Construct the b-parts for the systems Ax=b %
disp("  Construction of b");
b_cells = cell(1, 3);
for k = 1:conductor_num
    b_temp = zeros(quarters, 1); % b_temp is a temporary vector for each case %
                                 % this will be inserted into the cell array  %
    for i = 1:quarters
        if conductor_array(i,1) == k
            b_temp(i) = 1;
        end 
    end 
    b_cells{k} = b_temp; 
end 
    
% solve the systems Ax=b according to the number of conductors from the input file %
x_cells = cell(1,3);
disp("  Solving the systems");
for k=1:conductor_num
    x_cells{k} = A1\b_cells{k};
end

disp("  Calculation of Q arrays");
Q_array_temp = zeros(conductor_num, 1); % this is a temporary vector for the Q %
                                        % for each different system that has   %
                                        % been solved                          %    

Q_cells = cell(1, 3);
for i = 1:conductor_num
   for j = 1:quarters
        conductor_index = conductor_array(j,1);
        temp_x_vector = (x_cells{i});
        % Q = Q + x * surface %
        Q_array_temp(conductor_index) = Q_array_temp(conductor_index) + ( temp_x_vector(j)*(conductor_array(j, 8) - conductor_array(j, 2)) * (conductor_array(j, 5) - conductor_array(j, 3)));
   end 
   % Insert the computed values of Q of each conductor in the cell array and then %
   % initialize it with zeros for the next solution set                           %  
   Q_cells{i} = Q_array_temp; 
   Q_array_temp = zeros(conductor_num, 1);
end 
temp_x_vector = [];

% Initially compute the Cxx elements                    %
% In the case with 3 conductors, this loop will compute %
% the elements C11, C22 and C33                         %
 C_vector = zeros(conductor_num*2, 1);
 c_index = 0;
 for i = 1:conductor_num
     c_index = c_index + 1;
     for j = 1:conductor_num
        C_vector(c_index) = C_vector(c_index) + Q_cells{i}(j);
     end 
 end 
    
 % Then compute the Cxy elements                         %
 % In the case with 3 conductors, this loop will compute %
 % the elements C12, C13 and C23                         %
 division = floor(conductor_num/2);
 for i=1:conductor_num - division
    for j=1:conductor_num - i
        c_index = c_index + 1;
        C_vector(c_index) = - Q_cells{i+j}(i);   
    end 
 end 

disp("  C vector");
disp(C_vector);
fprintf("\n\n");

iteration = 1;
% Loop where the quarters are divided and the new systems Ax = b are %
% constructed and solved. Then the Q's for each case are computed    %
% in order to compute the C vector. Then based on the comparison     %
% the precvious iteration, the program terminates (threshold = 10^-3 %
% for the norms of C vector).                                        %
while (1)
    terminateFlag = 0;
    fprintf("Iteration : %d\n", iteration);
    
    iteration = iteration + 1;
    new_conductor_array = []; % empty the old array in order to store the newly divided quarters % 
    new_row = zeros(1, 9); % entry for the conductor array %
    
    disp("  Division of quarters");
    % The quarters need to be divided into 4 different new quarters % 
    for c = 1:quarters
        % Format of the rows in the conductor array                 %
        % row(1) = id                                               %
        % row(2) = x1 row(3) = y1       row(4) = x2 row(5) = y2     %
        % row(6) = x3 row(7) = y3       row(8) = x4 row(9) = y4     %
        
        % ------ 1st quarter ------ %
        new_row(:) = 0;
        new_row(1) = conductor_array(c, 1); % conductor id %
    
        % lower left x, y %
        new_row(2) = conductor_array(c, 2); 
        new_row(3) = conductor_array(c, 3);
        % upper left x, y %
        new_row(4) = conductor_array(c, 2); 
        new_row(5) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
    
        % upper right x, y %
        new_row(6) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2); 
        new_row(7) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
        % lower right x, y %
        new_row(8) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2); 
        new_row(9) = conductor_array(c, 3);
       
        new_conductor_array = [new_conductor_array; new_row];
    
        % ------ 2nd quarter ------ %
        new_row(:) = 0;
        new_row(1) = conductor_array(c, 1); % conductor id %
        
        % lower left x, y %
        new_row(2) = conductor_array(c, 2); 
        new_row(3) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
        % upper left x, y %
        new_row(4) = conductor_array(c, 4); 
        new_row(5) = conductor_array(c, 5);
    
        % upper right x, y %
        new_row(6) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2);
        new_row(7) = conductor_array(c, 5);
        % lower right x, y %
        new_row(8) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2); 
        new_row(9) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
        
        new_conductor_array = [new_conductor_array; new_row];
        
        % ------ 3rd quarter ------ %
        new_row(:) = 0;
        new_row(1) = conductor_array(c, 1); % conductor id %
        
        % lower left x, y %
        new_row(2) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2);
        new_row(3) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
        % upper left x, y %
        new_row(4) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2);
        new_row(5) = conductor_array(c, 5);
    
        % upper right x, y %
        new_row(6) = conductor_array(c, 6); 
        new_row(7) = conductor_array(c, 7);
        % lower right x, y %
        new_row(8) = conductor_array(c, 8); 
        new_row(9) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
        
        new_conductor_array = [new_conductor_array; new_row];
    
        % ------ 4th quarter ------ %
        new_row(:) = 0;
        new_row(1) = conductor_array(c, 1); % conductor id %
        
        % lower left x, y %
        new_row(2) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2);
        new_row(3) = conductor_array(c, 3);
        % upper left x, y %
        new_row(4) = (conductor_array(c, 8) - conductor_array(c, 2))/2 + conductor_array(c, 2);
        new_row(5) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
    
        % upper right x, y %
        new_row(6) = conductor_array(c, 8); 
        new_row(7) = (conductor_array(c, 5) - conductor_array(c, 3))/2 + conductor_array(c, 3);
        % lower right x, y %
        new_row(8) = conductor_array(c, 8);
        new_row(9) = conductor_array(c, 9);

        new_conductor_array = [new_conductor_array; new_row];
    end

    quarters = height(new_conductor_array);
    
    disp("  Number of Quarters:");
    disp(quarters);
    
    curr_center_x = 0;
    curr_center_y = 0;
    relative_center_x = 0;
    relative_center_y = 0;
    hx = 0;
    hy = 0;

    fprintf("  Construction of A (%d x %d)\n", quarters, quarters);
    A1 = zeros(quarters, quarters);

    for i = 1:quarters
        curr_center_x = (new_conductor_array(i, 8) - new_conductor_array(i, 2))/2 + new_conductor_array(i, 2);
        curr_center_y = (new_conductor_array(i, 5) - new_conductor_array(i, 3))/2 + new_conductor_array(i, 3);
        for j = 1:quarters
            hx = (new_conductor_array(j, 8) - new_conductor_array(j, 2));
            hy = (new_conductor_array(j, 5) - new_conductor_array(j, 3));
            if (i == j)
                A1(i, j) = (2*(hx * log((hy + sqrt(hx^2 + hy^2))/hx) + hy * log((hx + sqrt(hx^2 + hy^2))/hy))) / (4*pi*8.85*10^(-12));
            else
                relative_center_x = (new_conductor_array(j, 8) - new_conductor_array(j, 2))/2 + new_conductor_array(j, 2);
                relative_center_y = (new_conductor_array(j, 5) - new_conductor_array(j, 3))/2 + new_conductor_array(j, 3);
                A1(i, j) = ((hx * hy)/sqrt(((curr_center_x - relative_center_x)^2 + (curr_center_y - relative_center_y)^2))) / (4*pi*8.85*10^(-12));
            end 
        end
    end
        
    % Construct the b-parts for the systems Ax=b %
    disp("  Construction of b");
    for k = 1:conductor_num
        b_temp = zeros(quarters, 1); 
        for i = 1:quarters
            if new_conductor_array(i,1) == k
                b_temp(i) = 1;
            end 
        end 
        b_cells{k} = b_temp; 
    end 
        
    % Solve the systems Ax=b according to the number of conductors from the input file %
    disp("  Solving the systems");
    for k=1:conductor_num
        x_cells{k} = A1\b_cells{k};
    end
    
    A1 = [];
    
    [size1, size2] = size(A1);

    % Store the C's of the previous iteration before computing the new values %
    % The old values will be used to calculate the norms                      %        
    C_vector_old = C_vector;
    C_vector = [];

    % Compute the new values of Q %
    disp("  Calculation of Q arrays");
    Q_array_temp = zeros(conductor_num, 1); 
    for i=1:conductor_num
       for j = 1:quarters
            conductor_index = new_conductor_array(j,1);
            temp_x_vector = (x_cells{i});
            Q_array_temp(conductor_index) = Q_array_temp(conductor_index) + ( temp_x_vector(j)*(new_conductor_array(j, 8) - new_conductor_array(j, 2)) * (new_conductor_array(j, 5) - new_conductor_array(j, 3)));
       end 
       Q_cells{i} = Q_array_temp;
       Q_array_temp = zeros(conductor_num, 1);
    end 
    
    % Initially compute the Cxx elements %
    % In the case with 3 conductors, this loop will compute %
    % the elements C11, C22 and C33                         %
     C_vector = zeros(conductor_num*2, 1);
     c_index = 0;
     for i = 1:conductor_num
         c_index = c_index + 1;
         for j = 1:conductor_num
            C_vector(c_index) = C_vector(c_index) + Q_cells{i}(j);
         end 
     end 
    
     % Then compute the Cxy elements %
     % In the case with 3 conductors, this loop will compute %
     % the elements C12, C13 and C23                         %
     division = floor(conductor_num/2);
     for i=1:conductor_num - division
        for j=1:conductor_num - i
            c_index = c_index + 1;
            C_vector(c_index) = - Q_cells{i+j}(i);
            
        end 
     end 
    
    disp("  C vector (current)");
    disp(C_vector);

    disp("  Check if the difference of the norms of C is less than 10^3:");
    disp((norm(C_vector - C_vector_old))/norm(C_vector_old));
    if ((iteration ~= 1) && ((norm(C_vector - C_vector_old))/norm(C_vector_old)) <= 10^(-3))
        terminateFlag = 1;
        break;
    end

    if (terminateFlag == 1)
        break;
    end 

    conductor_array = new_conductor_array;

    % empty the cell arrays and vectors to reduce memory consumption %
    x_cells = {};
    b_cells = {};
    Q_cells = {};
    A1= [];
    Q_array_temp = [];
    b_temp = [];
    new_conductor_array = [];

    fprintf("\n");
end