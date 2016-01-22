function test_fun(a,varargin)
 
b = 2; % set defaults for b and c
c = 3;
 
extract_varargin; % override b and c if specified
 
fprintf('a is %d, b is %d, c is %d\n',a,b,c);