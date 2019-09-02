function [ y ] = LinearInterpolate( x,x1,y1,x2,y2 )
%UNTITLED12 Summary of this function goes here
%  The following function is a standard linear interpolation routine. W/e
%  Units
%1/13/16. Johnny Jaffee

%Inputs
    %x = the desired independent var. pt you want the corresponding y val
    %for. Falls between x1 and x2
    %x1 = the lower bound x value given.
    %y1 corresponding y value to given x1 lower bound
    %x2 = upper given x val
    %y2 = corresponding y value to x2
%Outputs
    %y = found value corresponding to x
    
%checks that the given x fall between x1 and x2. If it is below or equal to
%x1 sets y = y1, and V.V. Willmprint error message
if x < x1
    disp('Error: x is below x1. y set equal to y1.');
    y = y1;
else
if x > x2
    disp('Error: x is above x2. y set equal to y2.');
    y = y2;
else
if x == x1
    y = y1;
else
if x == x2
    y = y2
else
    %calculates the slope
    m = (y2-y1)/(x2-x1);
    %calculates y intercept
    b = y1-m*x1;
    %calculates y, using y = mx+b
    y = m*x+b;
end

end
end
end
end

