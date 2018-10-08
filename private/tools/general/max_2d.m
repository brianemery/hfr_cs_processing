function  [mx,r,c] = max_2d(x)
% MAX 2D - find global max of a matrix
% [mx,r,c] = max_2d(x)


[cx,j] = max(x);

[mx,c] = max(cx);

r = j(c);

end