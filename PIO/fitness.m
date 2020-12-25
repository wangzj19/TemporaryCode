function R = fitness(X,in)
s=0;
S=0;
% ���Ժ���һ��Sphere[-100,100]
% for i = 1:in
%     s = s + X(i).^2;
% end
% R = s;
% ���Ժ�������Rosenbrock Function[-10,10]
% for ii = 1:(in-1)
% 	xi = X(ii);
% 	xnext = X(ii+1);
% 	new = 100*(xnext-xi^2)^2 + (xi-1)^2;
% 	s = s + new;
% end
% R = s;
%���Ժ�������Zakharov Function****[-5, 10]
% for i = 1:in
%     xi = X(i);
%     s = s + xi^2;
%     S = S + 0.5*i*xi;
% end
% R = s + S^2 + S^4;
%���Ժ����ģ�McCormick Function[-3,4]
% x1 = X(1);
% x2 = X(2);
% term1 = sin(x1 + x2);
% term2 = (x1 - x2)^2;
% term3 = -1.5*x1;
% term4 = 2.5*x2;
% R = term1 + term2 + term3 + term4 + 1;
%���Ժ����壺Griewank Function[-600,600]
% for ii = 1:in
% 	xi = X(ii);
% 	s = s + xi^2/4000;
% 	S = S * cos(xi/sqrt(ii));
% end
% R = s -S + 1;
% ���Ժ�������Ackley Function��-40 40��
for i = 1:in
    s = s + X(i).^2;
    S = S + cos(2*pi*X(i));
end
R = -20*exp(-0.2*sqrt((1/in)*s))-exp((1/in)*S)+20+exp(1);
%���Ժ����ߣ�Holder Table Fuction**********ע�⣬�������������޸��������䣬
%���������Ա���x������Ϊ[-10,10],������������Сֵ
% x1 = X(1);
% x2 = X(2);
% fact1 = sin(x1)*cos(x2);
% fact2 = exp(abs(1 - sqrt(x1.^2+x2.^2)/pi));
% R = -abs(fact1*fact2);

%���Ժ����ˣ�Schaffer Function[-100,100]
% R =((sin(sqrt(X(:,1).^2+X(:,2).^2))).^2-0.5)./((1+0.001*(X(:,1).^2+X(:,2).^2)).^2)+0.5;
% %���Ժ����ţ�Booth Function[-10,10]
% x1 = X(1);
% x2 = X(2);
% term1 = (x1 + 2*x2 - 7)^2;
% term2 = (2*x1 + x2 - 5)^2;
% R = term1 + term2;
%���Ժ���ʮ��Branin Function*********����[-5,10]
% x1 = X(1);
% x2 = X(2);
% a = 1;b = 5.1/(4*pi^2);c = 5/pi;r = 6;ss = 10;t = 1/(8*pi);
% term1 = a * (x2 - b*x1^2 + c*x1 - r)^2;
% term2 = ss*(1-t)*cos(x1);
% R = term1 + term2 + ss;

%���Ժ����ţ�Shubert Function***x[-10��10]
% x1 = X(1);
% x2 = X(2);
% for i = 1:5
% 	new1 = i * cos((i+1)*x1+i);
% 	new2 = i * cos((i+1)*x2+i);
% 	s = s + new1;
% 	S = S + new2;
% end
% R = s * S;
%���Ժ���ʮ����Rastrigin
% for i = 1:in
%     s = s + X(i)^2 - 10*cos(2*pi*X(i));
% end
% R = s+300;
end