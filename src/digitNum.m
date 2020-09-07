function digitNum = digitNum(x)
a=mod(x,1); %the fraction part
div=0.1;
digitNum=0;
while abs(a-0)>1e-12
    a=mod(a,div);
    digitNum=digitNum+1;
end