clf
#load('func1.txt')
#load('func2.txt')
#load('K_d_arr1.txt')
#load('K_d_arr2.txt')
#load('K_d.txt')
#load('func.txt')
X = load('X.txt');
Z = load('Z.T.txt');
H = load('H.txt');
hold on
# figure('name','banana')
hFig = figure(1);
hSurf = surface(X, Z, H);
view(30, 50)
#while(ishandle(hFig)) sleep(1); end

#plot(K_d_arr1, func1, 'k')
#plot(K_d_arr2, func2, 'k')
#plot(K_d, func, 'diamond', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
#grid('on')
#xlabel('Kd')
#ylabel('Function Value')
#axis([K_d_arr1(1), K_d_arr2(length(K_d_arr2)), func2(1), func2(length(func2))])