java c
COSC2500/COSC7500—Semester 2, 2024
MODULE 3. ODES
Exercises—due 3:00 pm Friday 11th October
Required
The grade of 1–5 awarded for these exercises is based on the required exercises. If all of the required exercises are completed correctly, a grade of 5 will be obtained.
R3.1 From Sauer computer problems 6.1:
Using Euler’s method and a 4th-order Runge–Kutta method, solve the fol-lowing IVPs on the interval [0, 1] where y(0) = 1 and
a) y′ = t
b) y′ = 2(t + 1)y
c) y′ = 1/y2
for step sizes h = 0.1 × 2−k, for 0 ≤ k ≤ 5 (i.e., h = 0.1, 0.05, 0.025, ...). Plot the solutions for h = 0.1, 0.05, 0.025, comparing with the exact solution. Plot a log–log plot of the error at t = 1 as a function of the step size h.
Also compare with the solution obtained using Matlab’s ode45 or another adaptive-step Runge–Kutta code.
A fixed-step Runge–Kutta code, RK4.m, is available on Blackboard.
R3.2 Consider the systems of linear DEs from Sauer computer problems 7.1–7.3, questions 1–2:
a) i. y′′ = y + 3/2et, y(0) = 0, y(1) = 3/1e.
ii. y′′ = (2 + 4t2)y, y(0) = 1, y(1) = e.
b) i. y′′ = 3y − 2y′, y(0) = e3, y(1) = 1.
a) Solve these DEs using the shooting method.
b) Solve these DEs using the finite difference method. Use the finite dif-ference substitutions for the derivatives to convert the DEs to linear systems, and solve the linear systems using the backslash operator or otherwise.
R3.3 Discuss and critique the following Matlab and C codes. Briefly compare them with the code for Matlab’s ode45 (you can display it in the command window using type ode45 or open it in the editor with open ode45).
function [t ,y , last ] = RK4 (f , tspan , yi , dt )
% Standard RK4 algorithm .
% This code is a standard RK4 algorithm with a
% fixed step - size .
% It can be used to solve explicit first - order ODEs .
% The local truncation error is O( dt ^5).
%" Time " is the name for the x - axis variable .
%
% Input :
% f is a function handle for the first - order ODE
% - dy / dt = f(t ,y)
% - f is assumed to be a row vector .
% tspan is a vector that contains ti and tf ,
% e.g. tspan = [ti , tf ]
% ti is the initial time
% tf is the final time
% yi is the initial values for the ODE .
% - yi is assumed to be a vector .
% dt is the step size
%
% Output :
% t is a vector of time values .
% y is the approximate solution curve for the
% initial value problem of the ODE .
% last is the last y values - useful if you are
% using the last values for something
%
% Hint - call the function like this :
% [T ,Y] = RK4 (f , tspan ,yi , dt )
ti = tspan (1); tf = tspan (2);
num_steps = ceil (( tf - ti )/ dt );
% ceil forces it to be an integer
t = linspace ( ti , tf , num_steps +1). ’;
% creates time vector , then transposes
if size ( yi ,2) > 1
yi = yi . ’;
end
% Initialising y
y = [ yi , zeros ( size (yi ,1) , num_steps )];
% Application of the RK4 algorithm .
for n = 1: num_steps
k1 = dt * f ( t ( n ) , y (: , n ));
k2 = dt * f ( t ( n ) + 0.5* dt , y (: , n ) + 0.5* k1 );
k3 = dt * f ( t ( n ) + 0.5* dt , y (: , n ) + 0.5* k2 );
k4 = dt * f ( t ( n ) + dt , y (: , n ) + k3 );
y (: , n +1) = y (: , n ) + (1/6)*( k1 + 2* k2 + 2* k3 + k4 );
end
y = y . ’;
last = y ( end ,:);
end
void runge4 ( double x , double y [] , double step )
{
double h = step /2.0 , /* the midpoint */
t1 [ N] , t2 [N ] , t3 [N ] ,
k1 [ N] , k2 [N ] , k3 [N ] , k4 [ N ];
int i ;
for ( i =0; i < N ; i ++)
{ t1 [ i ]= y [ i ]+0.5*( k1 [ i ]= step * f (x , y , i ));}
for ( i =0; i < N ; i ++)
{ t2 [ i ]= y [ i ]+0.5*( k2 [ i ]= step * f ( x +h , t1 , i ));}
for ( i =0; i < N ; i ++)
{ t3 [ i ]= y [ i ]+ ( k3 [ i ]= step * f ( x +h , t2 , i ));}
for ( i =0; i < N ; i ++)
{ k4 [ i ]= step * f ( x + step , t3 , i );}
for ( i =0; i < N ; i ++)
{ y [ i ]+=( k1 [ i ]+2* k2 [ i ]+2* k3 [i ]+ k4 [ i ])/6.0;}
}
double f ( double x , double y [] , int i )
{
/* derivative of first equation */
if ( i ==0) return ( exp ( -2* x ) - 3* y [0]);
if ( i !=0) exit (2);
}
R3.4 Discuss, with critical analysis, a paper from the research literature.
A list of papers will be available on Blackboard if you have difficulty finding a suitable paper to discuss.
Additional
Attempts at these exercises can earn additional marks, but will not count towards the grade of 1–5 for the exercises. Completing all of these exercises does not mean that 6 marks will be obtained—the marks depend on the quality of the answers. It is possible to earn all 6 marks without completing all of these additional exercises.
A3.1 Estimate the error in a numerical solution obtained using Euler’s method from the convergence of the solution as the step size is changed. Then com-pare your estimate of the error with actual error determined from compar-ison with the analytical solution.
A3.2 Below is a model that is a variation of the Lotka–Volterra equations (i.e., predator–prey equations):
x′ = 1.3x − 0.5xy
y′ = 0.2xy − 0.5y
Choose an ODE solver and solve this system for both x0 and y0 equal to 12 between t = 0 and t = 500. Plot x and y vs t and then a x-y plot. Briefly comment on the behaviour of the system and whether this model seems realistic (e.g., to model a system consisting of predators and prey).
One possible way to solve a system in Matlab using ode45 is to create an anonymous function that returns a column vector. To make it work with most of the Matlab ODE solvers, the inputs should be t and an arbitrary value (say) u. Every incidence of x is replaced with the indexed variable u(1) and every代 写COSC2500/COSC7500—Semester 2, 2024 MODULE 3. ODESR
代做程序编程语言 incidence of y is replaced with u(2).
For reference on how you would represent the system in Matlab, the fol-lowing example is provided. For the system
x′ = xy + x
y′ = y
the corresponding anonymous function code is
func = @ (t , u ) [u (1)* u (2) + u (1); u (2)];
% You should not be solving this equation
% for the module exercise ,
% this is an example ...
A3.3 Investigate Sauer computer simulation 6.3.2 (the pendulum). For some spe-cific tasks, see computer problems 6.3.
A3.4 Investigate Sauer computer simulation 6.3.3 (orbital mechanics). For some specific tasks, see computer problems 6.3.
A3.5 Sauer reality check 6—The Tacoma Narrows bridge (Sauer pg 322 in 2E, 328 in 1E). As well as the suggested activities, try using Matlab’s ode45 solver.
A3.6 Simulate the motion of a paper airplane.
(Googling for "paper airplane" "differential equations" lift drag will provide a convenient starting point.)
A3.7 Consider the systems of non-linear DEs from Sauer computer problems 7.1–7.3, questions 3–4:
a) i. y′′ = 18y2, y(1) = 3/1, y(2) = 12/1.
ii. y′′ = 2e−2y(1 − t2), y(0) = 0, y(1) = ln2.
b) i. y′′ = ey, y(0) = 1, y(1) = 3.
ii. y′′ = sin y′, y(0) = 1, y(1) = −1.
Solve these using the finite difference method, using the non-linear code from Sauer or otherwise. You might like to compare your solutions with those obtained using the shooting method.
A3.3 Solve the DEs from the required exercises and/or A3.7 using the finite ele-ment method, collocation method, or similar method.
A3.4 Estimate the error in a numerical solution obtained using the finite dif-ference method from the convergence of the solution as the step size is changed. Then compare your estimate of the error with actual error de-termined from comparison with the analytical solution.
A3.5 Solve the system of equations in Gramotnev and Nieminen 1999. The ap-pendix in the paper describes a finite difference solution of the system. Can you solve it using the shooting method?
A3.6 Model a multi-species predator–prey relationship using a homogeneous linear matrix DE. If we begin with known populations of each species, we have an initial value problem; if If we know a mixture of final and initial populations, we have a BVP. (You might like to try N = 5.)
Programming hints
R3.2(a) For the shooting method, you’re combining and initial value problem solver and a root-finder. You can use whichever IVP solver (from this mod-ule or elsewhere) and root-finder (from Module 2 or elsewhere) you prefer. The only tricks are that
1. Your IVP solver will want a system of 1st order ODEs. You are starting with a second order ODE. So, the first step is to convert the 2nd order ODE into a system of 2 1st order ODEs, as described at the start of Module 4.
2. You need to write a function that the root-finder will find the roots of. For the shooting method, we guess the missing initial condition, and solve the system of ODEs as an IVP. We can then compare our solution to the other boundary condition (which we haven’t used yet). We want the difference between our solutions and this boundary condition to be zero. That’s our root-finding problem. So we need a function like
function final_error = bc_mismatch (x )
% Find the roots of this for R5 .1 shooting method
% Solve the IVP
[t , y ]= ode45 ( @r51_1a ,[0 1] ,[0 x ]);
% Compare with final BC
bc = (1/3) * exp (1);
final_error = y ( end ,1) - bc ;
return
The parameter passed to this function, x, is your guess for the missing initial condition. Since your root-finder will only want to pass one parameter to the function it’s finding roots of, it’s easiest to hard-code things like the ODE system function name, the boundary conditions, etc. in the file as in this example. Depending on which order you write your 2 1st order ODEs in, you might need to swap the order of the 2 boundary conditions passed to your IVP solver, and to compare the final value y(end,2).
As far as your root-finder cares, this function is just another function to find roots of. You can also calculate values for different x and plot a graph; this will be easiest to do using a loop to cycle through the different values of x.
R3.2(b) If you want to solve these as linear systems, you’ll need to write your own code. Once you have your linear system, solving it is easy. There are two steps you will need to do first:
1. Substitute the finite difference approximations into your ODE. Here is an example, with a complicated ODE:
Our ODE: y′′ = y′ + t sin t
Our finite difference approximation for y′′:
y′′n = (yn+1 − 2yn + yn−1)/h2
Note that yn = y(tn) and yn+1 = y(tn+1) = y(tn + h).
Let’s use the forward difference for y′: y′n = (yn+1 − yn)/h.
Finally, we’ll replace t by tn.
So, our ODE becomes:

We want to write this in our usual form. for a linear equation in y1, y2, etc., as An yn−1 + Bn yn + Cn yn+1 = Dn. Re-arranging, we get

2. Generate the linear system. The easy way is to put our initial and final boundary conditions as the 1st and last rows.
% How many points do we want ?
N = 10;
% Pre - allocate memory for matrix and vector
A = zeros (N , N );
c = zeros (N ,1);
tn = linspace (0 ,1 , N ); % Need to change for 2( a),
% since it goes to 1.5
h = t (2) - t (1);
% Boundary conditions
A (1 ,1) = 1;
c (1) = 0; % First BC
A (N , N ) = 1;
c ( N ) = exp (1)/3; % End BC
% Make the matrix
for n = 2:( N -1)
A (n ,n -1) = 1/ h ^2;
A (n , n ) = 1/ h - 2/ h ^2
A (n , n +1) = 2/ h ^2 - 1/ h ;
c ( n ) = tn ( n ) * sin ( tn ( n ));
end
Finally, solve and plot!







         
加QQ：99515681  WX：codinghelp  Email: 99515681@qq.com
