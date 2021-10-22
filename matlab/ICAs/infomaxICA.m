function [w,wz]=infomaxICA(x)

%runica 				% calculates wz, w and uu. The matrix x gets
				% overwritten by a sphered version of x. 
%%%%%%%%%%%%%%%Start RunICA
% Script to run ICA on a matrix of images. Original by Tony Bell. 
% Modified by Marian Stewart Bartlett.

%Assumes image gravalues are in rows of x. Note x gets overwritten.
%Will find N independent components, where N is the number of images.

%There must be at least 5 times as many examples (cols of x) as the
%dimension of the data (rows of x). 

N=size(x,1); P=size(x,2); M=N;	 %M is dimension of the ICA output
%spherex;                          % remove first and second order stats from x
% SPHEREX - spheres the training vector x.
%    Requires x, P, to be predefined, and defines mx, c, wz.

mx=mean(x');%subtracting mean\n');
x=x-(ones(P,1)*mx)';
c=cov(x');%'calculating whitening filter\n');
wz=2*inv(sqrtm(c));
x=wz*x;%whitening\n');
%%%%%End spherex
xx=inv(wz)*x;                     % xx thus holds orig. data, w. mean extracted.

%******** setup various variables
w=eye(N); count=0; perm=randperm(P); sweep=0; Id=eye(M);
oldw=w; olddelta=ones(1,N*M); angle=1000; change=1000;

%******** Train. outputs a report every F presentations.
% Watch "change" get small as it converges. Try annealing learning 
% rate, L, downwards to 0.0001 towards end.
% For large numbers of rows in x (e.g. 200), you need to use a low 
% learning rate (I used 0.0005). Reduce if the output blows 
% up and becomes NAN. If you have fewer rows, use 0.001 or larger.

B=50; L=0.001; F=5000; 
for I=1:200%I means iteration
    %sep96;
    x=x(:,perm);
    sweep=sweep+1; t=1;
    noblocks=fix(P/B);
    BI=B*Id;
    for t=t:B:t-1+noblocks*B,
        count=count+B;
        u=w*x(:,t:t+B-1);
        w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
    end;
end
% B=50; L=0.0003; F=5000; 
% for I=1:200
%     %sep96;
%     x=x(:,perm);
%     sweep=sweep+1; t=1;
%     noblocks=fix(P/B);
%     BI=B*Id;
%     for t=t:B:t-1+noblocks*B,
%         count=count+B;
%         u=w*x(:,t:t+B-1);
%         w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
%     end;
% end
% B=50; L=0.0002; F=5000; 
% for I=1:200
%     %sep96;
%     x=x(:,perm);
%     sweep=sweep+1; t=1;
%     noblocks=fix(P/B);
%     BI=B*Id;
%     for t=t:B:t-1+noblocks*B,
%         count=count+B;
%         u=w*x(:,t:t+B-1);
%         w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
%     end;
% end
% B=50; L=0.0001; F=5000; 
% for I=1:200
%     %sep96;
%     x=x(:,perm);
%     sweep=sweep+1; t=1;
%     noblocks=fix(P/B);
%     BI=B*Id;
%     for t=t:B:t-1+noblocks*B,
%         count=count+B;
%         u=w*x(:,t:t+B-1);
%         w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
%     end;
% end
