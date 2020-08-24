function Xhat=ppxa(Y,Sd,St, test_ind) %Xhat

global R

M = ones(size(Y));          % weight matrix or binary mask
M(test_ind) = 0;
 
R=opRestriction(numel(Y),find(M(:)==1));%opRestriction(numel(Y), test_ind);


global theta
theta=5;

noi=20;

global pp mu1 mu2 lamda 


%lambda=0.1;

Sd = preprocess_PNN(Sd,pp);
St = preprocess_PNN(St,pp);  

% Laplacian Matrices    
Dd = diag(sum(Sd)); Lr = Dd - Sd;  
if(det(Dd)==0)
    Dd=0.1*eye(size(Dd))+Dd;
end
Lr = (Dd^(-0.5))*Lr*(Dd^(-0.5));

Dt = diag(sum(St)); Lc = Dt - St;
if(det(Dt)==0)
    Dt=0.1*eye(size(Dt))+Dt;
end
Lc = (Dt^(-0.5))*Lc*(Dt^(-0.5));

X1=randn(size(Y));
X2=randn(size(Y));
X3=randn(size(Y));
X4=randn(size(Y));
X5=randn(size(Y));

Xhat2=randn(size(Y));

X_prev = (1/theta)*(X1+X2+X3+X4+X5);

%binary mask
%%Rop=opRestriction(numel(Y),find(M(:)==1));
%%R=opToMatrix(Rop); %(dim of R is #elements2keep X  #totalEelementsInX)

 %%A = eye(size(R,2)) + theta*(R'*R);

 
   for k=1:noi
       
       %Xhat1vec=lsqr( [sqrt(theta/2)*R; sqrt(1/2)*I] ,[sqrt(theta/2)*Yvec;sqrt(1/2)*X1vec] ); Xhat1=reshape(Xha1tvec, size(Y))
        %Xhat1=mldivide(sqrt(theta/2)*R; sqrt(1/2)*eye(size(R))] , [sqrt(theta/2)*Y;sqrt(1/2)*X1]);
       
        % you need to solve (I + theta*R?*R)x = theta*R?*y + x_1
       
        %%b = theta*R'*Rop(Y(:),1) + X1(:);
        b= theta*R(R(Y(:),1),2) + X1(:);
        %%xhat1 = conjgrad(A,b);%lsqr(A,b,1e-100,30);%lsqr(A,b);
       
        xhat1=lsqr(@fun_handle,b ,[],20);
        Xhat1 = reshape(xhat1, size (M));  
       
        %%%% made Xhat1 upate iterative
       
        [U, S, V] = svd(X2);
        sigma = sign(S).*max(0,abs(S)-(theta*lamda/2) );  
        Xhat2=U*(sigma)*V';
       
       
        Xhat3=min(max(X3,0),1);
        %T=min(max(X3,0),1); Xhat3= (T>=0.5);
       
        
        Xhat4 = pinv( eye(size(Lr,1))+2*theta*mu1*Lr)*X4;
        %Xhat4 = sylvester(mu1*eye(size(X4,1)),theta*Lc,(1/2)*X4); %
       
        Xhat5 = X5*pinv(2*theta*mu2*Lc + eye(size(Lc,1)));
        %Xhat5=(sylvester(mu2*eye(size(X5,2)),theta*Lr,(1/2)*X5'))';
       
        Xhat=(1/theta)*(Xhat1+Xhat2+Xhat3+Xhat4+Xhat5); 
        
       
        X1=X1+(2*Xhat-X_prev-Xhat1); %   mid var
        X2=X2+(2*Xhat-X_prev-Xhat2);
        X3=X3+(2*Xhat-X_prev-Xhat3);
        X4=X4+(2*Xhat-X_prev-Xhat4);
        X5=X5+(2*Xhat-X_prev-Xhat5);
       
        X_prev = Xhat ; %Xhat1 0.5087
       
        %keep this for debug
       %crit(k) = 1/2*norm(Rop(Y(:),1)-R*Xhat(:))^2 + ...
           %lamda.*sum(abs(svd(Xhat(:)))) + mu1*trace(Xhat'*Lr*Xhat)+mu2*trace(Xhat*Lc*Xhat');
           crit(k) = 1/2*norm(R(Y(:),1)-R(Xhat(:),1))^2 + ...
           lamda.*sum(abs(svd(Xhat(:)))) + mu1*trace(Xhat'*Lr*Xhat)+mu2*trace(Xhat*Lc*Xhat');
    end
   
  Xhat= X_prev;
 Xhat(Xhat<0)=0.01;
 
 
end


function y =fun_handle(x,direction)

global R theta
Rt=opTranspose(R);

I=opDiag(numel(x), 1);

    if strcmp(direction,'notransp') %y=A*x

        y= (theta*R( (R(x,1)) ,2) ) + (I(x,1) ); 
   %or y= (theta*Rt( (R(x,1)) ,1) ) + (I(x,1) ); 
       
       
    elseif strcmp(direction,'transp') %y=A'*x

        y= (theta*R( (R(x,1)) ,2) ) + (I(x,1) );
    end
   
end

function [S,p]=preprocess_PNN(S,p)
%preprocess_PNN sparsifies S by keeping, for each drug/target, the "p"
% nearest neighbors (NNs) and discarding the rest. 

    NN_mat = zeros(size(S));
    for j=1:length(NN_mat)
        [~,indx] = sort(S(j,:),'descend');
        indx = indx(1:p+1);     % keep drug/target j and its "p" NNs
        NN_mat(j,indx) = 1;
    end
    NN_mat = (NN_mat+NN_mat')/2;
    S = NN_mat .* S;

end
