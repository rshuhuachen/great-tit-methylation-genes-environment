"""
Autor: Michelle Memelink
Date: 03-01-2022
Description: This script filters the results from the mQTL analysis.
The trans file is filtered by the distance between the SNP
site and the CpG site.
"""
try:
    import gzip
    import re
except ImportError:
    raise ImportError("An error has occurred when importing a module")


def ReadInfo(path_info):
    """
    the function reads in the info file and extracts only the information that
    applies to the python script in which this function is run (filter_matrixeqtl.py).
    This concerns input and output pathways and file names. This information
     is then stored in a dictionary. If the info-file does not contain
    information for this python script, the script will stop.
    :param path_info: contains the pathway of the file containing all input and
    output paths of the files used and created.
    :return dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    """
    dict_path_info = {}
    x = False
    try:
        with open(path_info, "r") as file:
            for line in file:
                if line.__contains__("filter_matrixeqtl.py"):
                    x = True
                elif line == "\n":
                    x = False
                elif x:
                    line = line.split("=")
                    dict_path_info[line[0]] = line[1].replace("\n", "")
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    if len(dict_path_info) == 0:
        print("The info-file does not contain any information associated with "
              "filter_matrixeqtl.py.")
        exit()

    return dict_path_info

def ReadTrans(dict_path_info):
    """
    This function creates the "output2_trans_filtered.txt" file. This file only
    contains trans SNP-CpG relationships where they are more than 5 million base
    pairs apart.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: -
    """
    index_snp, index_gene = "", ""

    try:
        file_filterd = open(dict_path_info['output'] + "output2_trans_filtered.txt", "w")
        with open(dict_path_info['trans_data'], "r") as file:
            for line in file:
                if line.__contains__("snps"):
                    file_filterd.write(line)
                    line = line.split(",")
                    for x in range(len(line)):
                        if line[x].__contains__("snps"):
                            index_snp = x
                        elif line[x].__contains__("gene"):
                            index_gene = x
                else:
                    CheckRelationship(line, index_snp, index_gene, file_filterd)
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


def CheckRelationship(line, index_snp, index_gene, file_filterd):
    """
    This function checks per SNP-CpG relationships whether the SNP and the CpG
    are at least 5 million base pairs apart. If this is the case, the line from
    the original trans file is saved in the new file.
    :param line: is a string and consisting of a line from the file from
    function "ReadTrans(..)"
    :param index_snp: Is a string and contains the index where SNP is given
    in the file
    :param index_gene: Is a string and contains the index where CpG site is
    given in the file
    :param file_filterd: this specifies the file to be created in the function
    "ReadTrans(..)"
    :return: -
    """
    try:
        split = line.split(",")
        snp = split[index_snp].split("_")
        gene = split[index_gene].split("_")
        if snp[0] + "_" + snp[1] == gene[0] + "_" + gene[1]:
            difference = abs(int(snp[2].replace('"', ''))-int(gene[2]
                                                              .replace('"', '')))
            if difference > 5000000:
                file_filterd.write(line)
        else:
            file_filterd.write(line)
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


if __name__ == '__main__':
    # Change the pathway below to the correct custom pathway:
    path_info = "/home/info.txt"
    dict_path_info = ReadInfo(path_info)

    ReadTrans(dict_path_info)
