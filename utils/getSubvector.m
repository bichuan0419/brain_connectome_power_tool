function subvector = getSubvector(Clist, CID, i)
% split the vector Clist into subvectors using the values of CID as the lengths of these subvectors. 
% Then, when given an index i, the corresponding subvector from Clist is selected.
    if i == 1
        start_idx = 1;
    else
        start_idx = sum(CID(1:i-1)) + 1;
    end
    end_idx = sum(CID(1:i));
    
    subvector = Clist(start_idx:end_idx);
end
