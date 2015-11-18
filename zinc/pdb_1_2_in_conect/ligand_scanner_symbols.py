from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp

pdbpath = '/Users/rafalpwiewiora/pdb/pdb/*/*.ent.gz'
cutoff = 0.3
metal_name = 'ZN'
ppn = 8


def ligand_scanner(file):
    
    #print(file)
    
    ok_file_count, error_file_count, metal_count, files_with_metal_count = 0, 0, 0, 0
    extra_in_cutoff_count_by_atom, extra_in_CONECT_count_by_atom, both_extra_CONECT_higher_count_by_atom = 0, 0, 0
    both_extra_cutoff_higher_count_by_atom, both_extra_both_equal_by_atom, extra_in_cutoff_count_by_file = 0, 0, 0
    extra_in_CONECT_count_by_file, both_extra_count_by_file, both_equal_count_by_atom, both_equal_count_by_file = 0, 0, 0, 0
    CONECT_ligand_numbers = []
    extra_in_cutoff_atoms_residues = []
    CONECT_only_ligand_numbers = []
    CONECT_and_extra_ligand_numbers = []
    
    
    metal_ligand_dict_by_cutoff = {}
    metal_ligand_dict_by_CONECT = {}
    metal_ligand_dict_in_cutoff_not_CONECT = {}
    metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT = {}
    metal_ligand_dict_in_CONECT_not_cutoff = {}
    extra_in_cutoff_for_file_analyzed = False
    extra_in_CONECT_for_file_analyzed = False
    both_extra_for_file_analyzed = False
    both_equal_for_file_analyzed = False
    dictionary_of_process_counts = {}
    
    # initialize the results dictionary (dictionary_of_process_counts)        
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['both_equal_count_by_atom'] = both_equal_count_by_atom
    dictionary_of_process_counts['extra_in_cutoff_count_by_atom'] = extra_in_cutoff_count_by_atom
    dictionary_of_process_counts['extra_in_CONECT_count_by_atom'] = extra_in_CONECT_count_by_atom
    dictionary_of_process_counts['both_extra_CONECT_higher_count_by_atom'] = both_extra_CONECT_higher_count_by_atom
    dictionary_of_process_counts['both_extra_cutoff_higher_count_by_atom'] = both_extra_cutoff_higher_count_by_atom
    dictionary_of_process_counts['both_extra_both_equal_by_atom'] = both_extra_both_equal_by_atom
    dictionary_of_process_counts['both_equal_count_by_file'] = both_equal_count_by_file
    dictionary_of_process_counts['extra_in_cutoff_count_by_file'] = extra_in_cutoff_count_by_file
    dictionary_of_process_counts['extra_in_CONECT_count_by_file'] = extra_in_CONECT_count_by_file
    dictionary_of_process_counts['both_extra_count_by_file'] = both_extra_count_by_file     
    dictionary_of_process_counts['CONECT_ligand_numbers'] = CONECT_ligand_numbers
    dictionary_of_process_counts['extra_in_cutoff_atoms_residues'] = extra_in_cutoff_atoms_residues
    dictionary_of_process_counts['CONECT_only_ligand_numbers'] = CONECT_only_ligand_numbers
    dictionary_of_process_counts['CONECT_and_extra_ligand_numbers'] = CONECT_and_extra_ligand_numbers
    
    
    # load file
    try:
        traj = md.load_pdb(file)
        ok_file_count += 1
        dictionary_of_process_counts['ok_file_count'] = ok_file_count
    except:
        error_file_count += 1
        dictionary_of_process_counts['error_file_count'] = error_file_count
        return None

    
    # Select atoms and pairs
    topo = traj.topology
    metal_select_name = 'name %s and resname %s' % (metal_name, metal_name)
    metal_atoms = topo.select(metal_select_name)
    metal_all_pairs = topo.select_pairs(metal_select_name, 'symbol O or symbol N or symbol S or symbol Cl')
    
    # No metal - skip, metal - add to count 
    if not metal_atoms.size:
        return dictionary_of_process_counts
    
    files_with_metal_count += 1
    metal_count += len(metal_atoms)
    
    # Prep dictionary with metal indices as keys
    for i in metal_atoms:
        metal_ligand_dict_by_cutoff[i] = []
        metal_ligand_dict_by_CONECT[i] = []
        metal_ligand_dict_in_cutoff_not_CONECT[i] = []
        metal_ligand_dict_in_CONECT_not_cutoff[i] = []
        metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i] = []
        
        
    # Compute distances
    try:
        metal_all_distances = md.compute_distances(traj, metal_all_pairs, periodic = False)
    except:
        error_file_count += 1
        dictionary_of_process_counts['error_file_count'] = error_file_count
        #return dictionary_of_process_counts   
        return None
        
    # Populate the metal_ligand_dict_by_cutoff with atoms closer than cutoff
    for i in range(len(metal_all_distances[0])):
        if metal_all_distances[0, i] < cutoff:
            if metal_all_pairs[i, 0] in metal_atoms:
                metal_ligand_dict_by_cutoff[metal_all_pairs[i, 0]].append(metal_all_pairs[i, 1])
            elif metal_all_pairs[i, 1] in metal_atoms:
                metal_ligand_dict_by_cutoff[metal_all_pairs[i, 1]].append(metal_all_pairs[i, 0])
                
    # Look at the CONECT record through topology.bonds
    for i in metal_atoms:
        for bond in topo.bonds:
            if bond[0].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[1].index)
            elif bond[1].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[0].index)
                
    for i in metal_atoms:
        CONECT_ligand_numbers.append(len(metal_ligand_dict_by_CONECT[i]))  
        if len(metal_ligand_dict_by_CONECT[i]) == 1 or len(metal_ligand_dict_by_CONECT[i]) == 2:
            print(file)
            return None
                  
                
    # Compare the contents of metal_ligand_dict_by_cutoff and metal_ligand_dict_by_CONECT, correct for cutoff recognizing extra atoms for residues already in CONECT
    for i in metal_atoms:
        metal_ligand_dict_in_cutoff_not_CONECT[i] = [x for x in metal_ligand_dict_by_cutoff[i] if x not in metal_ligand_dict_by_CONECT[i]]
        metal_ligand_dict_in_CONECT_not_cutoff[i] = [x for x in metal_ligand_dict_by_CONECT[i] if x not in metal_ligand_dict_by_cutoff[i]]    
    
    for i in metal_atoms:
        for j in metal_ligand_dict_in_cutoff_not_CONECT[i]:
            if not any(topo.atom(j).residue == topo.atom(k).residue for k in metal_ligand_dict_by_CONECT[i]):
                metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i].append(j)
    
    
    # compare the cutoff and CONECT data and make conclusions
    for i in metal_atoms:
        
        if not metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i] and not metal_ligand_dict_in_CONECT_not_cutoff[i]:
            both_equal_count_by_atom += 1
            CONECT_only_ligand_numbers.append(len(metal_ligand_dict_by_CONECT[i]))
            
            if not both_equal_for_file_analyzed:
                both_equal_count_by_file += 1
                both_equal_for_file_analyzed = True
        
        elif metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i] and not metal_ligand_dict_in_CONECT_not_cutoff[i]:
            extra_in_cutoff_count_by_atom += 1
            extra_in_cutoff_atoms_residues.append([])
            for j in metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]:
                extra_in_cutoff_atoms_residues[-1].append((topo.atom(j).residue.name, topo.atom(j).name))
                
            CONECT_and_extra_ligand_numbers.append(len(metal_ligand_dict_by_CONECT[i]) + len(metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]))    
                
            if not extra_in_cutoff_for_file_analyzed:
                extra_in_cutoff_count_by_file += 1
                extra_in_cutoff_for_file_analyzed = True
        
        elif metal_ligand_dict_in_CONECT_not_cutoff[i] and not metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]:
            extra_in_CONECT_count_by_atom += 1
            if not extra_in_CONECT_for_file_analyzed:
                extra_in_CONECT_count_by_file += 1
                extra_in_CONECT_for_file_analyzed = True
                 
        elif metal_ligand_dict_in_CONECT_not_cutoff[i] and metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]:
            
            if len(metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]) > len(metal_ligand_dict_in_CONECT_not_cutoff[i]):
        
                both_extra_cutoff_higher_count_by_atom += 1
                if not both_extra_for_file_analyzed:
                    both_extra_count_by_file += 1
                    both_extra_for_file_analyzed = True
            
            elif len(metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]) < len(metal_ligand_dict_in_CONECT_not_cutoff[i]):
                    
                both_extra_CONECT_higher_count_by_atom += 1
                if not both_extra_for_file_analyzed:
                    both_extra_count_by_file += 1
                    both_extra_for_file_analyzed = True
                     
            elif len(metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]) == len(metal_ligand_dict_in_CONECT_not_cutoff[i]):
                
                both_extra_both_equal_by_atom += 1
                if not both_extra_for_file_analyzed:
                    both_extra_count_by_file += 1
                    both_extra_for_file_analyzed = True
    
    # update the results dictionary (ok_file_count and error_file_count have already been done)

    dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['both_equal_count_by_atom'] = both_equal_count_by_atom
    dictionary_of_process_counts['extra_in_cutoff_count_by_atom'] = extra_in_cutoff_count_by_atom
    dictionary_of_process_counts['extra_in_CONECT_count_by_atom'] = extra_in_CONECT_count_by_atom
    dictionary_of_process_counts['both_extra_CONECT_higher_count_by_atom'] = both_extra_CONECT_higher_count_by_atom
    dictionary_of_process_counts['both_extra_cutoff_higher_count_by_atom'] = both_extra_cutoff_higher_count_by_atom
    dictionary_of_process_counts['both_extra_both_equal_by_atom'] = both_extra_both_equal_by_atom
    dictionary_of_process_counts['both_equal_count_by_file'] = both_equal_count_by_file
    dictionary_of_process_counts['extra_in_cutoff_count_by_file'] = extra_in_cutoff_count_by_file
    dictionary_of_process_counts['extra_in_CONECT_count_by_file'] = extra_in_CONECT_count_by_file
    dictionary_of_process_counts['both_extra_count_by_file'] = both_extra_count_by_file   
    dictionary_of_process_counts['CONECT_ligand_numbers'] = CONECT_ligand_numbers
    dictionary_of_process_counts['extra_in_cutoff_atoms_residues'] = extra_in_cutoff_atoms_residues
    dictionary_of_process_counts['CONECT_only_ligand_numbers'] = CONECT_only_ligand_numbers
    dictionary_of_process_counts['CONECT_and_extra_ligand_numbers'] = CONECT_and_extra_ligand_numbers
                
                                                          
    #print(dictionary_of_process_counts)
    #return dictionary_of_process_counts              
    return None                                        
                     
def ligand_scanner_all_database(pdbpath):
    
    dictionary_of_database_results = {}
    dictionary_of_database_results['ok_file_count'] = 0
    dictionary_of_database_results['error_file_count'] = 0
    dictionary_of_database_results['files_with_metal_count'] = 0
    dictionary_of_database_results['metal_count'] = 0
    dictionary_of_database_results['both_equal_count_by_atom'] = 0
    dictionary_of_database_results['extra_in_cutoff_count_by_atom'] = 0
    dictionary_of_database_results['extra_in_CONECT_count_by_atom'] = 0
    dictionary_of_database_results['both_extra_CONECT_higher_count_by_atom'] = 0
    dictionary_of_database_results['both_extra_cutoff_higher_count_by_atom'] = 0
    dictionary_of_database_results['both_extra_both_equal_by_atom'] = 0
    dictionary_of_database_results['both_equal_count_by_file'] = 0
    dictionary_of_database_results['extra_in_cutoff_count_by_file'] = 0
    dictionary_of_database_results['extra_in_CONECT_count_by_file'] = 0
    dictionary_of_database_results['both_extra_count_by_file'] = 0  
    dictionary_of_database_results['CONECT_ligand_numbers'] = []
    dictionary_of_database_results['extra_in_cutoff_atoms_residues'] = []
    dictionary_of_database_results['CONECT_and_extra_ligand_numbers'] = []
    dictionary_of_database_results['CONECT_only_ligand_numbers'] = []
    
    pool.map(ligand_scanner, glob.iglob(pdbpath))
    
    return None

    #for dictionary_of_process_counts in pool.map(ligand_scanner, glob.iglob(pdbpath)):
    #    for i in dictionary_of_process_counts:
#            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
            
#    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    dictionary_of_database_results = ligand_scanner_all_database(pdbpath)
    

# Write results
f = open('ligand_scanner_results.txt', 'w')
f.write(str(dictionary_of_database_results))
f.close()
