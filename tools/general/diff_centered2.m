function du2 = diff_centered2(u,x)
% DIFF CENTERED - Centered difference for computing 2nd derivatives 
% du2 = diff_centered2(u,x)
%
% Computes the 2nd derivative of a vector using the  finite
% difference approximation. Second order accurate (c.f. Laveque 2007 
% "Finite Difference Methods for Ordinary and Partial Differential 
% Equations"). Edge derivatives are computed using the left and right sided
% approximations with second order accuracy.
% 
% The method assumes equal spacing of grid x, but this implementation of it
% does not require equal spacing. The centered difference is order (deltax^2)
% accurate, thus accuracy goes as the square of the spacing. An easy way to
% increase accuracy (in most practical applications) is to use finer grid.
%
% INPUTS
% u - the function of x. Can be a matrix of independent rows, columns
%     associated with x (NOTE DIFFERENCE WITH DIFF.M, which is the matrix
%     of row diffs). 
% x - the grid (usually equal spacing)
%
% OUTPUTS
% du2/dx2 - same size as inputs
%
% EXAMPLE
% % Need to be careful using this with degrees vs radians:
% x = 1:180; % degrees
% % Must convert all inputs to the same grid
% d2u = diff_centered2(sind(x),x*pi/180);
% plot(x,cosd(x),'-bo'), hold on
% plot(x,du,'r*')
%
% SEE ALSO 
% del2.m


% 10 Feb 2017

% Note: Might look into CFD notes for alternative methods if more accuracy
% is needed, cf compact_finite_diff.m which uses Lele 92 method but needs
% to be generalized


% check for test case
if strcmp('--t',u), test_case, return, end

% Preallocate
du2 = NaN(size(u));
x = repmat(x(:).',size(du2,1),1);

% Get first point from formula found on internet (fd.pdf)
% http://shahroodut.ac.ir/fa/download.php?id=1111113740
% Right sided second order scheme
du2(:,1) = ( 2*u(:,1) - 5*u(:,2) + 4*u(:,3) - u(:,4) ) ./ ((x(:,2)-x(:,1)).^2) ;

% Centered Approximation, Eqn 1.13
i = 2:(size(u,2)-1);
du2(:,i) = ( u(:,i+1) - 2*u(:,i) + u(:,i-1) )./ (( x(:,i)-x(:,i-1)).^2) ;

% Get last point
% Left sided second order scheme
i = size(u,2);
du2(:,i) = ( 2*u(:,i) - 5*u(:,i-1) + 4*u(:,i-2) - u(:,i-3) ) ./ ((x(:,i)-x(:,i-1)).^2) ;



end


function test_case
% FROM MATHWORKS:
% Use the diff function to approximate partial derivatives with the syntax 
% Y = diff(f)/h, where f is a vector of function values evaluated over 
% some domain, X, and h is an appropriate step size.
% 
% For example, the first derivative of sin(x) with respect to x is cos(x),
% and the second derivative with respect to x is -sin(x). You can use diff 
% to approximate these derivatives.

% Mathworks test, modified for matrix inputs
% NOTE:    diff(X), for a matrix X, is the matrix of row differences,
% [X(2:n,:) - X(1:n-1,:)].

h = 0.001;       % step size
X = -pi:h:pi;    % domain
f = sin(X);      % range
f = [f; f; f;];


x = 1:180; % degrees

% Must convert all inputs to the same grid
d2u = diff_centered2(sin(X),X);
plot(X,-sin(X),'r*'), hold on
plot(X,d2u,'-b.')

% error is order 10e-8
figure
plot(X,abs(d2u+sin(X)))



keyboard

% this function 
d2f = diff_centered2(f,X);

% looks good?
h3 = plot(X,d2f,'-g.');

keyboard

legend([h2(1) h1(1)  h3(1)],'u(x)', 'du/dx','FEM version')


% % COMPARE ERRORS?
% figure
% plot(df(1,:) - cos(X)), hold on
% plot(Y(1,:) - cos(X(1:size(Y,2))),'-r.')





% % TEST X NOTE
% % sind converts x into radians, so the inputs need to be consistent?
% 
% x = 1:0.01:180;
% du = diff_centered(sind(x),x*pi/180);
% plot(x,cosd(x),'-bo'), hold on
% plot(x,du,'r*')
% 
% keyboard

% TEST 2
% Laveque 2007, pg 4
% Example 1.1. Let u(x) = sin(x) and x = 1; 
% thus we are trying to approximate u'(x 0:5403023. 
x = 1;
u = sin(x);

h = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3];

disp('h D3')

for i = 1:numel(h)
    
    t = [-2 -1 0 1 2]*h(i);
    
    du = diff_centered(sin(x + t), t);
    
    err(i) = 0.5403023 - du(3); 
    
    % disp([num2str(h(i)) ' ' num2str(err)])
end

% Compare with figure 1.2
h1 = loglog(h,err,'-bo'); hold on

h0 = loglog(h, [9.0005e-4 2.251e-4 9.005e-6 2.2513e-6 9.005e-8],'r*');

axis([1e-4 1 1e-10 1e-2])

legend('diff_centered.m','Laveque Fig 1.2')

keyboard


end