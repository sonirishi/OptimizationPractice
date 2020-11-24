A = [3 -1; 6 3; -1 0; 0 -1; 4 1];
B = [0;40;0;0;16];
C=[30 10];
cvx_begin
    variable x(2)
    maximize C*x
    subject to
        A*x <= B;
cvx_end
    
        