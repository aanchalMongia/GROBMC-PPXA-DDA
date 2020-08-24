function getParameters(classifier,ds)
%
% This function retrieves optimal parameter values based on the supplied
% algorithm, the CV scenario used and the dataset on which prediction will
% be performed.
%
% INPUT:
%  classifier:  algorithm to be used for DDA prediction
%  ds:          dataset (4:NR, 3:GPCR, 2:IC, 1:E)
%
% OUTPUT:
%  params:   optimal parameter values for given inputs
%
global pp mu1 mu2 lamda
    
    switch classifier

        case 'ppxa'
           
                    switch(ds)
                        
                        case 1
                            pp=5;
                            lamda=0.1;%0.5; auc0.95
                            mu1=0.05;%0.001;%0.1;
                            mu2=0.1;%0.001;%0.1;
                        case 2
                            pp=5;
                            lamda=0.1;%0.5;%1;
                            mu1=0.05;%0.001;%0.1;
                            mu2=0.1;%0.001;%0.1;
                   
                    end
                    
            end
            

end