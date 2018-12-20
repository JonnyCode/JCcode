set(0, 'DefaultAxesFontName','Palatino')
set(0, 'DefaultAxesFontSize', 16)

M = [3 -2 1
	 2  1 3
	-1  2 1];

figure(1)
Plot3DVector([3 2 -1], 'b');
hold on
Plot3DVector([-2 1 2], 'g');
Plot3DVector([1 3 1], 'r');
hold off
xlabel('x')
ylabel('y')
zlabel('z')
rotate3d on

