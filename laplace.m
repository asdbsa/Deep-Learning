clear all;
%数据读取
d=fopen('a.mat');
data=fread(d,inf, 'int16');
fclose(d);
x=zeros(1,10480);
for j=1:10480
    x(j)=data(j+100);
end
x=x-mean(x);
%定义laplace小波；
fs = 10240;
for f=1:50;
   for ks=1:40;
      for nn=1:50;
         lp(f,ks,nn)=exp(-((ks/200)/sqrt(1-(ks/200)^2))*((f*10+1500)*2*pi)*nn/fs)*exp(-i*(f*10+1500)*2*pi*nn/fs);
         clp(f,ks,nn)=conj(lp(f,ks,nn));
         alp(f,ks,nn)=abs(lp(f,ks,nn));
      end
   end
end

%做相关滤波
y=zeros(1,2048);
z=zeros(1,2048);
k=zeros(1,2048);
y(1)=0; %记录最大的kr;
z(1)=0; %记录最大的f;
k(1)=0; %记录最大的ks;
for tao=1:2048;
    for m=tao:(tao+49);
        a(:,:,m)=clp(:,:,m-tao+1).*x(m);  %x(t)与lp小波共轭相乘；
        ax(m)=abs(x(m));                  %|x(m)|;
        ax2(m)=(ax(m))^2;                 %|x(m)|^2;
        alp2(:,:,m)=(alp(:,:,m-tao+1)).^2;%|lp|^2;
    end
    x2=sqrt(sum(ax2));   
    lp2(:,:)=sqrt(sum(alp2,3)); 
    b(:,:)=sum(a,3);
    for f=1:50;
        for ks=1:40;
            c=lp2(f,ks);
            kr=sqrt(2)*abs(b(f,ks))/(x2*c);
        
            if kr>y(tao);
               y(tao)=kr;
               z(tao)=f*10+1500;
               k(tao)=ks/200;
            end
        end
     end
end

%直线拟合剔除直流，
t=1:2048;
A=polyfit(t,y,2);
aa=A(1);bb=A(2);cc=A(3);
for i=1:2048
    y(i)=y(i)-(aa*i^2+bb*i+cc);
end
for i=1:10;
    y=y-mean(y);
end
figure;
t=1:2048;
subplot(211);
plot(t/10240,x(t));hold on; ylabel('x(t)');xlabel('时间t/s');
subplot(212);
plot(t/10240,y(t));hold on; ylabel('kr(i)');xlabel('时间t/s');
%直线拟合剔除直流后，求包络谱；
nfft=10240;
hkr=hilbert(y);
p=abs(fft(hkr,nfft));
figure;
plot((0:nfft/2-1)/nfft*fs,p(1:nfft/2));hold on; ylabel('幅值');xlabel('f /Hz');

    
