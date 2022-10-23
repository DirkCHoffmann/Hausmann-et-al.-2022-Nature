function randpkloc = rand_matrix_CircularShift(pklocs)

randpkloc = zeros(size(pklocs,1),size(pklocs,2));

for roi=1 : size(pklocs,2)
    shift = round(rand * size(pklocs,1));
    randpkloc(:,roi) = circshift(pklocs(:,roi),shift);
end

end