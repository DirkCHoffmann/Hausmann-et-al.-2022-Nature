function randpkloc = rand_matrix_Scrambled(pklocs)

randpkloc = zeros(size(pklocs,1),size(pklocs,2));

for roi=1 : size(pklocs,2)
    randframe = randperm(size(pklocs,1));
    for frame=1 : size(pklocs,1)
       randpkloc(frame, roi) = pklocs(randframe(frame), roi);
    end
    clearvars randframe
end

end