clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%USER DEFINED INPUT 
% The u matrix is user defined input describing the given exisitng machine 
% where u(:,1)represents the types of machine and u(:,2) is the machine 
% number for each type. u(i,j)is sorted from minimum to maximum
u = [1 1;
    2 1; 2 2; 
    3 1; 3 2; 3 3;
    4 1;4 2;4 3;4 4]; 
   

z = [2 4 ;
    4 14 ;4 12 ;  
    7 12 ;7 10; 7 8;
    9 4;9 7;9 11;9 14];
   

Ps = [16 8 8 4 4 4 2 2 2 2];

% z matrix is the user defined coordinates of each machine where z(:,1)
% represent the horizontal (x) coordinate and z(:,2) represents the
% vertical (y) coordinate. 
% z = [2 4 ;4 6 ;4 4 ;5 4 ;5 2; 6 3];

% function [Throughput, TotalCost] = Design_opt_project(u,z)
% n matrix defines the total number of machines given
n = length(u);

% Ps is the production score (unit/time) and it represents the production 
% capacity of each machine. Ps is a scalar value with an assummption that 
% the capacity of all the machine is same,  
% Ps = [8 4 4 2 2 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step to computate the throughput based on all the feasible adjacent matrix 
% Throughput is a vector which give the total thoughput for each type of 
% machine 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The below code explains the maximum possible ones that can be present in
% the adjacent matrix for the given number and types of machines
A0 = ones(n,n);
for i = 1:n
    for j = 1:n
        if (j <= i)
            A0(i,j) = 0;
         elseif (u(j,1) - u(i,1) > 1) 
            A0(i,j) = 0;
         elseif u(j,1) - u(i,1) == 0
             A0(i,j) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the number of ones in A0 matrix calculating all the possible
% combinations 
k = nnz(A0);
I = dec2bin(0:2^k-1) - '0'; %All the possible combination of one [2^k x k] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on all the possible combinations filtering out the non feasible
% combinations in the binary matrix
xx = sum(u(:,1) == 2); % Defining total number of type 2 machines
yy = sum(u(:,1) > 1); % Defining total number of machine for type 2 and greater  
At = zeros(2^k,k);
% Perform elimination on non feasible Adjacency matrix
for i = 1:2^k
    if (sum(I(i,1:xx)) >= 1) && (sum(I(i,:)) <= yy)  
       At(i,:) = I(i,:); 
    end
end
At(all(~At,2),:) = [];

% At = @(k)A_total(k);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input the combination of binary matrix into the adjacent matrix and 
% perform further filtering to obtain the possible feasible Adjacent matrix.

s = length(At);
[row, col] = find(A0 > 0);
b = [row, col];
Ac = cell(s,1); %Input the combination of ones into the Adjacent Matrix 
for i = 1:s
    A = zeros(n,n);
    a = At(i,:);
    for j = 1:k
        A(b(j,1),b(j,2)) = a(j);
    end
    Ac{i} = A;
end

% Eliminate all the non feasible adjacent matrix 
for i = 1:s
    aa = Ac{i};
    for j = 1:n
        if sum(aa(:,j)) > 1 % Eliminate matrix that has sum of 
%             element column >1 
            Ac{i} = [];
        end
    end
    for m = 1:yy
        if sum(aa(:,1+m)) == 0 && sum(aa(1+m,:)) > 0 % Eliminate matrix 
%           that contain discrete connection in the production line
            Ac{i} = [];
        end
    end
end
Ac = Ac(~cellfun('isempty',Ac));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Po defines the production output for each machine. The output depends on 
% the connection between the machines that is defined in the adjacent matrix 
r = size(Ac,1);
Po = cell(r,1);
for i = 1:r
    aa = Ac{i};
    P_out = zeros(n,1);
    for j = 1:n
    P_out(j) = 1./sum(aa(j,:));
    if P_out(j) == Inf
        P_out(j) = 0;
    else
        P_out(j) = P_out(j);
    end
    end
    Po{i} = P_out;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As defines all the possible connections from each of the machines based on 
% the defined adjacent matrix. 
As = cell(r,1);
q = max(u(:,1)) - 1;
for i = 1:r
    aa = Ac{i};
    A_sum = 0;
    for j = 1:q
    A_sum = A_sum + aa^j;
    end
    As{i} = A_sum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S is a binary matrix in which the diagonal element of the matrix defines 
% whether there is an output exist for each machine. 1 defines there is an
% output exist on the machine, otherwise value of 0 indicates no output at machine 

S = cell(r,1);
for i = 1:r
    aa = Ac{i};
    So = zeros(n);
    for j = 1:n
    if (aa(j,:)*ones(n,1) == 0) && (aa(:,j)'*ones(n,1) == 1)
        So(j,j) = 1;
    else
        So(j,j) = 0;
    end
    end
    S{i} = So;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T is a cell that store all the binary matrix that represent the output 
% of each machine for each production line. The production line is defined 
% if and only if there a non zero value exist for each column of Matrix Z. 
% The non zero value is defined as output on machine. 
T = cell(r,1);
for i = 1:r
    T{i} = As{i} * S{i};
end

% Pw is a cell that stores matrix that represent the production weight 
% resulted from the product of T matrix and Po that gives the output weight 
% of each machine at each production line
Pw = cell(r,1);
for i = 1:r
    Pw{i} = T{i} .* Po{i};
end

% Temporary convert the zero value of Pw to one
Pw1 = cell(r,1);
for i = 1:r
    aa = Pw{i};
    for j = 1:n
        for m = 1:n
            if aa(j,m) == 0
               aa(j,m) = 1;
            end
        end
    end
        Pw1{i} = aa;
end

% Pww is a [1 x n] vector represent the output of each machine in unit/time
Pww = cell(r,1);
for i = 1:r
    aa = Pw1{i};
    bb = Pw{i};
    Pww{i} = prod(aa);
    cc = Pww{i};
    for j = 1:n
        if sum(bb(:,j)) == 0
           cc(j) = 0;
         end
    end
    Pww{i} = cc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the throughput for each possible combination of the adjacent
% matrix. Throughput is in the vector form [1 x mtype] which give the total
% thoughput for each type of machine 
mtype = max(u(:,1)); % Total number of machine type
Throughput = cell(r,1);
for j = 1:r
    aa = Pww{j}.*Ps;
    bb = Throughput{j};
    for i = 1:mtype
        a_out = find(u(:,1) == i);
        bb(i) = sum(aa(a_out));
    end
    Throughput{j} = bb;
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
A_distance = cell(r,1);
for i = 1:r
    aaa = Ac{i};
    A_distance{i} = aaa.*D;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit cost is the costs per unit length of the conveyer required to setup
% the connection between machine 
Unitcost = 15;

% Calculating the cost required for each adjacent matrix 
Cost = cell(r,1);
for i = 1:r
    aaaa = A_distance{i};
    Cost{i} = aaaa * Unitcost;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Cost defines the cost required for each adjacent matrix 
TotalCost = zeros(r,1);
for i = 1:r
   aaaaa = Cost{i};
    TotalCost(i) = sum(aaaaa(:));
end
% end

Throughput = cell2mat(Throughput);
% Result = Throughput;
% Result(:,1) = TotalCost(:,1);
% scatter3(Result(:,2), Result(:,3), Result(:,1))
% xlabel('A');
% ylabel('B');
% zlabel('Cost');
% labels = num2str((1:r)','%d');
% text(Result(:,2), Result(:,3), Result(:,1), labels, 'horizontal','left', 'vertical','bottom')
% figure
% scatter(Result(:,2), Result(:,1))
% labels = num2str((1:r)','%d');
% text(Result(:,2), Result(:,1), labels, 'horizontal','left', 'vertical','bottom')
% xlabel('A')
% ylabel('Cost')
% figure
% scatter(Result(:,3), Result(:,1))
% labels = num2str((1:r)','%d');
% text(Result(:,3), Result(:,1), labels, 'horizontal','left', 'vertical','bottom')
% xlabel('B')
% ylabel('Cost')