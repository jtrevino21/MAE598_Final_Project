u = [1 1 ;
    2 1; 2 2; 
    3 1; 3 2; 3 3;
    4 1; 4 2; 4 3;4 4] ; 
   

z = [2 4 ;
    4 14; 4 12 ;  
    7 12; 7 10; 7 8;
    9 4;9 7;9 11;9 14] ;   

Ps = [16 8 8 4 4 4 2 2 2 2];
n = length(u);
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

TP = @(x) Throughput_function(x,u,Ps,A0,n,k);
C = @(x) Cost_function(x,u,z,A0,n,k);
fitnessfcn = @(x)[TP(x), C(x)];
% rng default % For reproducibility
options = gaoptimset('PopulationType','bitstring','PopulationSize',2000);%'Generations',1000,'PopulationSize',300);
[x,fval] = gamultiobj(fitnessfcn,k,[],[],[],[],[],[],options);

[m,n] = size(fval);
for i = 1:m
    for j = 1:n
        if isnan(fval(i,j))
           fval(i,j) = 0;
        end
    end
    if sum(fval(i,:)) == 0
       x(i,:) = 0;
    end
end
fval(all(~fval,2),:) = [];
fval(:,1) = [];
x(all(~x,2),:) = [];

        
