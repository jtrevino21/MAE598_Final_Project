function Throughput = Throughput_function(x,u,Ps,A0,n,k)
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
% Ps = [8 4 4 2 2];

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
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Po defines the production output for each machine. The output depends on 
% the connection between the machines that is defined in the adjacent matrix 

P_out = zeros(n,1);
    for j = 1:n
    P_out(j) = 1./sum(A(j,:));
    if P_out(j) == Inf
        P_out(j) = 0;
    else
        P_out(j) = P_out(j);
    end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As defines all the possible connections from each of the machines based on 
% the defined adjacent matrix. 
As = zeros(n,1);
q = max(u(:,1)) - 1;

A_sum = zeros(n);
    for j = 1:q
    A_sum = A_sum + A^j;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S is a binary matrix in which the diagonal element of the matrix defines 
% whether there is an output exist for each machine. 1 defines there is an
% output exist on the machine, otherwise value of 0 indicates no output at machine 

    S = zeros(n);
    for j = 1:n
    if (A(j,:)*ones(n,1) == 0) && (A(:,j)'*ones(n,1) == 1)
        S(j,j) = 1;
    else
        S(j,j) = 0;
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T is a cell that store all the binary matrix that represent the output 
% of each machine for each production line. The production line is defined 
% if and only if there a non zero value exist for each column of Matrix Z. 
% The non zero value is defined as output on machine. 
T = A_sum * S;

% Pw is a cell that stores matrix that represent the production weight 
% resulted from the product of T matrix and Po that gives the output weight 
% of each machine at each production line

Pw = T.* P_out;

% Temporary convert the zero value of Pw to one
Pw1 = Pw;
    for j = 1:n
        for m = 1:n
            if Pw(j,m) == 0
               Pw1(j,m) = 1;
            end
        end
    end
    
% Pww is a [1 x n] vector represent the output of each machine in unit/time
Pww = zeros(1,n);
    for j = 1:n
        if sum(Pw(:,j)) == 0
           Pww(:,j) = 0;
        else 
            Pww(:,j) = prod(Pw1(:,j));
         end
    end
 Pww = Pww.*Ps;   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the throughput for each possible combination of the adjacent
% matrix. Throughput is in the vector form [1 x mtype] which give the total
% thoughput for each type of machine 
mtype = max(u(:,1)); % Total number of machine type
Throughput = zeros(1,mtype);
     
    for i = 1:mtype
        a_out = find(u(:,1) == i);
        Throughput(:,i) = sum(Pww(a_out));
    end
end