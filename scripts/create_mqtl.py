"""
Autor: Michelle Memelink
Date: 03-01-2022
Description: This function reads the files containing the cis and trans
SNP-CpG relations. Subsequently, all relations that are located on the
chromosomes are shown in an mQTL graph which is also developed by this script.
"""
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("An error has occurred when importing a module")


def ReadInfo(path_info):
    """
    the function reads in the info file and extracts only the information that
    applies to the python script in which this function is run (create_mqtl.py).
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
                if line.__contains__("create_mqtl.py"):
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
              "create_mqtl.py.")
        exit()

    return dict_path_info


def ReadInfoFile(dict_path_info):
    """
    This function reads the file which contains the following information:
    chromosome name and length of the chromosome. Then this length is stored
    in a list.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return info_length: is a list containing lists containing (2D list) the
    length of the chromosome and also the total length of the beginning of the
    genome, as: [[400, 7060], [340, 7460], .....]
    """
    info_length = []
    index_refseq, index_size = "", ""
    count = 0

    try:
        with open(dict_path_info['info'], "r") as file:
            for line in file:
                line = line.split("\t")
                for x in range(len(line)):
                    if line[x].__contains__("RefSeq"):
                        index_refseq = x
                    if line[x].__contains__("Size"):
                        index_size = x
                if index_refseq != "" and index_size != "" \
                        and not line[index_refseq].__contains__("RefSeq"):
                    count += int(line[index_size])
                    info_length.append([line[index_refseq], count])
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return info_length


def ReadResults(list_pathways, info_length):
    """
    This function reads into the files containing the cis and the trans SNP-CpG
    relations. Subsequently, the SNP sites and CpG sites where they are
    located on the genome are stored in lists. However, both the SNP and the
    CpG must be on a chromosome for both to be stored.
    :param list_pathways: is a list containing the pathways and file names to be
    read, as: ["/home/output2_cis_filtered.txt", "/hom/output2_trans_filtered.txt"]
    :param info_length: is a list containing lists containing (2D list) the
    length of the chromosome and also the total length of the beginning of the
    genome, as: [[400, 7060], [340, 7460], .....]
    :return list_snps: is a list containing all SNP locations where they are
    located on the genome, as: [400, 565665, .....]
    :return list_meth: is a list containing all CpG sites where they are
    located on the genome, as: [400, 565665, .....]
    """
    list_snps = []
    list_meth = []
    index_snp, index_meth = "", ""

    try:
        for file in list_pathways:
            with open(file, "r") as file:
                for line in file:
                    if line.startswith('""'):
                        line = line.split(",")
                        for x in range(len(line)):
                            #dit kan misschien ook wel zodner for loop, dat je index ophaalt
                            if line[x].__contains__("snps"):
                                index_snp = x
                            if line[x].__contains__("gene"):
                                index_meth = x
                    elif index_snp != "" and index_meth != "" and not line.__contains__("NW_"):
                        line = line.split(",")
                        split_snp = line[index_snp].replace('"', '').split("_")
                        split_meth = line[index_meth].replace('"', '').split("_")
                        for list in range(len(info_length)):
                            if info_length[list][0] == split_snp[0] + "_" + split_snp[1]:
                                if list != 0:
                                    index = list - 1
                                    list_snps.append(
                                        int(info_length[index][1]) + int(
                                            split_snp[2]))
                                else:
                                    list_snps.append(int(split_snp[2]))
                            if info_length[list][0] == split_meth[0] + "_" + split_meth[1]:
                                if list != 0:
                                    index = list - 1
                                    list_meth.append(
                                        int(info_length[index][1]) + int(
                                            split_meth[2]))
                                else:
                                    list_meth.append(int(split_meth[2]))
        print("The number of SNP-CpG relationships found to show: " + str(len(list_snps)))
        print("The number of unique SNP-locations found: " + str(
            len(set(list_snps))))
        print("The number of unique CpG-sites found: " + str(len(set(list_meth))))
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return list_snps, list_meth


def CreatePlot(list_snps, list_meth, info_length, dict_path_info):
    """
    This function creates a graph that shows where all SNP-CpG relationships
    are located on the genome, by creating an mQTL graph.
    :param list_snps: is a list containing all SNP locations where they are
    located on the genome, as: [400, 565665, .....]
    :param list_meth: is a list containing all CpG sites where they are
    located on the genome, as: [400, 565665, .....]
    :param info_length: is a list containing lists containing (2D list) the
    length of the chromosome and also the total length of the beginning of the
    genome, as: [[400, 7060], [340, 7460], .....]
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: -
    """
    try:
        length = len(info_length) - 1
        fig, ax = plt.subplots()
        ax.set_xlim(0, info_length[length][1])
        ax.set_ylim(0, info_length[length][1])
        fig.set_size_inches(30, 30, forward=True)
        ax.scatter(list_snps, list_meth)
        plt.locator_params(axis="x", nbins=20)
        plt.yticks(fontsize=150)
        plt.xticks(fontsize=100)
        ax.tick_params(axis="x", colors="white")
        ax.tick_params(axis="y", colors="white")

        for list in info_length:
            plt.axvline(x=int(list[1]), color='r', linestyle='--')
            plt.axhline(y=int(list[1]), color='r', linestyle='--')

        ax.set_xlabel("SNP location",
                      fontsize=35)
        ax.set_ylabel('Methylation location (CpG-site)', fontsize=35)
        ax.set_title("SNP-CpG relations",
                     fontsize=50)
        fig.savefig(dict_path_info['output'] + 'snp_cpg_relations.png', dpi=100)
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


if __name__ == '__main__':
    # Change the pathway below to the correct custom pathway:
    path_info = "/home/info.txt"
    dict_path_info = ReadInfo(path_info)

    info_length = ReadInfoFile(dict_path_info)
    list_pathways = [dict_path_info['cis_data'], dict_path_info['trans_data']]
    list_snps, list_meth = ReadResults(list_pathways, info_length)
    CreatePlot(list_snps, list_meth, info_length, dict_path_info)
