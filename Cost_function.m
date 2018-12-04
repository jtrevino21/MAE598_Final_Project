function TotalCost = Cost_function(x,u,z,A0,n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%USER DEFINED INPUT 
% The u matrix is user defined input describing the given exisitng machine 
% where u(:,1)represents the types of machine and u(:,2) is the machine 
% number for each type. u(i,j)is sorted from minimum to maximum
% u = [1 1 ;2 1; 2 2; 3 1; 3 2];

% z matrix is the user defined coordinates of each machine where z(:,1)
% represent the horizontal (x) coordinate and z(:,2) represents the
% vertical (y) coordinate. 
% z = [2 4 ;4 6 ;4 4 ;5 4 ;5 2];

% function [Throughput, TotalCost] = Design_opt_project(u,z)
% n matrix defines the total number of machines given
%n = length(u);

% Ps is the production score (unit/time) and it represents the production 
% capacity of each machine. Ps is a scalar value with an assummption that 
% the capacity of all the machine is same,  
%Ps = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step to computate the throughput based on all the feasible adjacent matrix 
% Throughput is a vector which give the total thoughput for each type of 
% machine 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The below code explains the maximum possible ones that can be present in
% the adjacent matrix for the given number and types of machines
% A0 = ones(n,n);
% for i = 1:n
%     for j = 1:n
%         if (j <= i)
%             A0(i,j) = 0;
%          elseif (u(j,1) - u(i,1) > 1) 
%             A0(i,j) = 0;
%          elseif u(j,1) - u(i,1) == 0
%              A0(i,j) = 0;
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Based on the number of ones in A0 matrix calculating all the possible
% % combinations 
% k = nnz(A0);

% x = [1 , 1 , 1 , 0,  1, 0];

% Based on all the possible combinations filtering out the non feasible
% combinations in the binary matrix
xx  = sum(u(:,1) == 2); % Defining total number of type 2 machines
yy = sum(u(:,1) > 1); % Defining total number of machine for type 2 and greater 

% Perform elimination on non feasible Adjacency matrix

if (sum(x(1,1:xx)) >= 1) && (sum(x(1,:)) <= yy)  
   x(1,:) = x(1,:); 
else 
    x(1,:) = zeros(1,k);
end 

[row, col] = find(A0 > 0);
b = [row, col];

A = zeros(n,n);
    for j = 1:k
        A(b(j,1),b(j,2)) = x(j);
    end
 
% Eliminate all the non feasible adjacent matrix 

   
    for j = 1:n
        if sum(A(:,j)) > 1 % Eliminate matrix that has sum of 
%             element column >1 
            A = zeros(n,n);
        end
    end
    for m = 1:yy
        if sum(A(:,1+m)) == 0 && sum(A(1+m,:)) > 0 % Eliminate matrix 
%           that contain discrete connection in the production line
            A = zeros(n,n);
        end
    end

  if sum(A)==0
    A(:,:)=inf;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Steps for calculating the setup cost for each possible combination of the
% adjacent matrix. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the user defined coordinates, using the Eucledian distance formula
% calculating the distance score between each machines. D matrix is the
% distance score.

D = zeros(n,n);
for i = 1:n
    for j = 1:n
        D(i,j) = norm(z(i,:) - z(j,:));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the distance of the connecting machine from each adjacent matrix
% It is the multiplication of each possible combination of adjacent matrix 
% with D


A_distance = A.*D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit cost is the costs per unit length of the conveyer required to setup
% the connection between machine 
Unitcost = 15;

% Calculating the cost required for each adjacent matrix 

Cost = A_distance * Unitcost;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Cost defines the cost required for each adjacent matrix 
TotalCost = sum(Cost(:));
end