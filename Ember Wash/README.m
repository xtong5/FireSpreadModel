%changes
%
% 
%% fileLine.m
% line 71: "spot", ignite spots at 
% (41,151) (81,151) (121,151) (161,151)  
% (51,101) (101,101) (151,101)
% (67,51) (134,51)
% the initial state is defined as 'spot' in geom.m, line 58
% 
% line 94-108: the first step 
% the (velx,vely) are zeros so the first cos and sin are set to be zeros. 
%
% line 128-132: calculate the hypotenuse then cos and sin using previous
% step's (velx,vely). Then compute new (velx,vely)
%
%% solvers.m
% line 207-211: compute the velocities. 
% rdbino = binornd(1,0.01,o.N,o.N);
% velx = psix + etay + cos.*exprnd(o.s,o.N,o.N)*rdbino;
% vely = psiy - etax + o.s + sin.*exprnd(o.s,o.N,o.N)*rdbino;