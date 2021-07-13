function Q=cat2bin(Tem,type)

% This function transfrom a categorical pattern to binary pattern

% Input:
%       Tem: A matirx or image with categorical varibale. 
%       The categorical variable should be  define by numerical value like 1, 2, 3, etc

%       type: this is the category for which transfomarion need to be
%       performed. Suppose type=2, All category '2' in Tem will be replaced
%       by 1 and non-'2' will be replaced by 0;

% OUTPUT:
%       Q: The is a binary matrix or image generated from Tem

    k=find(Tem(:)==type);
    k1=find(Tem(:)~=type);
    Q=Tem;
    Q(k)=1;
    Q(k1)=0;
    
    
    
    