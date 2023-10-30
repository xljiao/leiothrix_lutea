# -*- coding: utf-8 -*-
import time


def main_func():
    now_s = time.time()

    # gene list
    gene_lst = get_lines_from_file('gene_lst.txt')

    # list of raw fasta files
    file_list = get_lines_from_file('file.lst')

    global_map = {}

    # Adjust the original data format
    for file_name in file_list:
        with open(file_name+".fas", "r") as origin_file:
            lines = origin_file.readlines()
            new_lines = []
            for line in lines:
                if line.startswith('>'):
                    temp_list = [line]
                    new_lines.append(temp_list)
                else:
                    new_lines[len(new_lines) - 1].append(line)

            result_lines = []
            for new_line in new_lines:
                temp = new_line[0]
                for item in gene_lst:
                    if item in temp:
                        new_line.insert(0, exe_str(temp, item))
                        new_line.pop(1)
                        result_lines.append(new_line)
                        break

            global_map[file_name.split('.')[0]] = result_lines

    global_map = convert_map(global_map)
    for gene_name, item in global_map.items():
        write_result(item, gene_name + '.fas')

    print(u'done, total generate {} files, ues time {} seconds'.format(str(len(global_map)), str(time.time() - now_s)))


def exe_str(line, sub_str):
    # Remove redundant strings. For examples, '>NC_0255561_bb_bb_bb_bb; 9593-9772; +; nad3_1' -> '>nad3_1'
    return '{}{}'.format('>', line[line.find(sub_str):])


def get_lines_from_file(file_name):
    with open(file_name, "r") as target_file:
        lines = target_file.readlines()
        new_lines = []
        # Replaces special characters in the line
        for line in lines:
            if '\r\n' in line:
                new_lines.append(line.replace('\r\n', ''))
            elif '\n' in line:
                new_lines.append(line.replace('\n', ''))
            else:
                new_lines.append(line)
        target_file.close()
        return new_lines


def write_result(nested_list, file_name):
    result_file = open(file_name, 'w')
    for inner_list in nested_list:
        for line in inner_list:
            result_file.write(line)


def convert_map(global_list):
    '''
     {
        'name1': [list1],
        'name2': [list2],
        'name3': [list3]
     }
    '''
    gene_map = {}
    for name, item_list in global_list.items():
        for item in item_list:
            gene_names = get_gene_names(item[0])
            gene_name = gene_names[0]
            new_name = ''
            if len(gene_names) > 1:
                new_name = '{}_{}'.format(name, gene_names[1])
            item.insert(0, '{}{}{}'.format('>', name if new_name == '' else new_name, '\n'))
            item.pop(1)
            if gene_name in gene_map:
                gene_map[gene_name].append(item)
            else:
                gene_map[gene_name] = [item]
    return gene_map


def get_gene_names(line):
    # Get gene names, for examples '>nad4l' -> 'nad4l'
    # In particular, when there is an underscore, return a list of size 2, for examples '>nad3_1' -> ['nad3', '1']
    return line.replace('>', '').replace('\n', '').split('_')


if __name__ == '__main__':
    main_func()

