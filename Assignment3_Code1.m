A = [50 30; 24 22; -1 0; 0 -1]
B = [2400;2100;-45;-5]
cvx_begin
    variable x(2)
    maximize sum(x)
    subject to
        A*x <= B;
cvx_end
    
        