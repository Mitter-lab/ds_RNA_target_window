"""
Fasta file format import module.
Reads file --> returns dictionary {header:sequence}
Sequence is capitalised.  Not checking for correct """

def get_ref_f_strand(ref_filename):
    """returns the provided strand only
    """

    ref_dict={} #the stored dictionary
    key='' #Sequence key
    full_len_seq = '' #Sequence strand

    with open(ref_filename, 'rU') as loaded_ref:

        for line in loaded_ref:
            
            clean_line=line.strip()
            
            if line[0] == '>' and full_len_seq == '':
                key = clean_line
            elif line[0] == '>' and full_len_seq != '':
                ref_dict.update({key:full_len_seq})
                key = clean_line
                full_len_seq = ''
            
            elif line[0] == '' and full_len_seq != '':
                ref_dict.update({key:full_len_seq})
                key=clean_line
                full_len_seq = ''   
            
            elif line[0] =='':
                pass
            
            else: 
                full_len_seq += clean_line.upper()
        
    	ref_dict.update({key:full_len_seq})



    return ref_dict