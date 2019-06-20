function [error] = getRMSE(Res,Test)
%RMSE returns the root mean squared error for the non-zeros values in Test.

        nnz = Test > 0;
        error = sqrt(sum((Res(nnz) - Test(nnz)).^2) / sum(nnz(:) > 0));                 
           
end
