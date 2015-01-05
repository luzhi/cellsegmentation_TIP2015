function [InterpolationPointsXYOnLine] = GetInterpolationPointsXY(CentroidXY, BoundaryXY, numPointsOnLine)
% Get (x,y) of Interpolation Points on Each Line.

xCentroid = CentroidXY(1);
yCentroid = CentroidXY(2);

xBoundary = BoundaryXY(1);
yBoundary = BoundaryXY(2);


LineLength = sqrt(double((xCentroid-xBoundary)^2 + (yCentroid-yBoundary)^2));
LineLength = round(LineLength);

LengthStep = LineLength / numPointsOnLine;

InterpolationPointsXYOnLine = zeros(numPointsOnLine + 1, 2);

InterpolationPointsXYOnLine(1,1) = xCentroid;
InterpolationPointsXYOnLine(1,2) = yCentroid;
InterpolationPointsXYOnLine(numPointsOnLine + 1,1) = xBoundary;
InterpolationPointsXYOnLine(numPointsOnLine + 1,2) = yBoundary;

for i = 2 : numPointsOnLine
    x_i = (xCentroid * (LineLength + 1 - i * LengthStep) + xBoundary * ( i * LengthStep - 1)) / LineLength;
    y_i = (yCentroid * (LineLength + 1 - i * LengthStep) + yBoundary * ( i * LengthStep - 1)) / LineLength;
    
    InterpolationPointsXYOnLine(i,1) = x_i;
    InterpolationPointsXYOnLine(i,2) = y_i;        
end
end