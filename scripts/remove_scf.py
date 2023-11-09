# -*- coding: utf-8 -*-
'''
# handel fasta files，keep some scaffolds or chromosomes
# need to provide the fasta file to process, the scaffold or chr list to keep, and the generated file name
# Usage : python remove_scf.py raw.fa scf.list fin.fa
'''

import sys
import time

fasta_file_name = sys.argv[1]
scf_file_name = sys.argv[2]
write_file_name = sys.argv[3]

# process_entry
def main_func():
    print(fasta_file_name, scf_file_name, write_file_name)

    now_s = time.time()
    group_list = []

    with open(fasta_file_name, "r") as fasta_file:
        read_list = fasta_file.readlines()
        fasta_file.close()
        for line in read_list:
            if line.startswith('>'):
                group_list.append([line])
            else:
                group_list[len(group_list) - 1].append(line)

    with open(scf_file_name, "r") as num_file:
        num_list = num_file.readlines()
        num_file.close()

    num_map = {}
    for line in num_list:
        num = line.replace('\n', '')
        num_map[num] = num

    for line in num_map:
        print(line)

    for item in group_list:
        print(item[0].replace('>', '').replace('\n', '')) 
    write_list = []
    with open(write_file_name, "w") as write_file:
        for item in group_list:
            num = item[0].replace('>', '').replace('\n', '')
            if num in num_map:
                write_list += item
        write_file.writelines(write_list)
        write_file.close()

    print('use time {}s, write {} lines'.format(str(time.time() - now_s), str(len(write_list))))


if __name__ == '__main__':
    print('aaa')
    main_func()
