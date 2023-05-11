
function [p95]=cumpdf(data_x)
    p_level=0.95;
    [pdf,bin]=hist(data_x,1:max(data_x));
    pdf=pdf./sum(pdf);
    pdf_cumsum=cumsum(pdf);
    temp_var=find(pdf_cumsum<=p_level);
    p95 = interp1(pdf_cumsum(temp_var(end):(temp_var(end)+1)),bin(temp_var(end):(temp_var(end)+1)),p_level,'spline');
end