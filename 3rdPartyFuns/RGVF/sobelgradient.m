function B = sobelgradient(A)
%calculate the gradient map of A, sqrt(Ax^2+Ay^2)
[m,n]=size(A);
C = BoundMirrorExpand(A);

h = fspecial('sobel');
v=h';

ResX = filter2(h,C,'valid');
ResY = filter2(v,C,'valid');

D = sqrt(ResX.^2 + ResY.^2);

MaxGra = max(max(D));
MinGra = min(min(D));

B = 255*(D(:,:)-MinGra) / (MaxGra - MinGra);
