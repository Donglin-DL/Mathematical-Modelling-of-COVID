function [out] = logLikeNB(data, mean, p, m, n)
%The function "logLikeNB" returns the log combined likelihood of a
%negative binomial probability model with m data sets (n observations in 
%each data set). The negative binomial probability model is parameterized
%in terms of the "mean" number of counts and the probability "p", where the
%variance is equal to the mean divided by p, 0 < p < 1.

sum1 = 0;

%m is the number of data sets
for j = 1:m

    %n is the number of observations (assumes the same number of 
    %observations in each data set)
    for i = 1:n

        r_temp = (p(j)*mean(j, i))/(1 - p(j));

        sum1 = sum1 + gammaln(data(j, i) + r_temp) - gammaln(data(j, i) + 1) - gammaln(r_temp) + r_temp*log(p(j)) + data(j, i)*log(1 - p(j));

    end

end

out = sum1;

end