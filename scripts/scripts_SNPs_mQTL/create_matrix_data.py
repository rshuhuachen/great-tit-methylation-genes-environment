"""
Autor: Michelle Memelink
Date: 03-01-2022
Description: By means of this script all input data required for the R-package
"matrixEQTL" can be generated. Five files are created: "Covariates.txt",
"SNP.txt", "GE.txt", "snpsloc.txt" and "geneloc.txt". Several input files are
required for this. The input files should contain general information
(gender and age), SNP data and methylation data.
"""

def ReadInfo(path_info):
    """
   the function reads in the info file and extracts only the information that
   applies to the python script in which this function is run
   (create_matrix_data.py). This concerns input and output pathways and file
   names. This information is then stored in a dictionary. If the info-file
   does not contain information for this python script, the script will stop.
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
                if line.__contains__("create_matrix_data.py"):
                    x = True
                elif line == "\n":
                    x = False
                elif x:
                    line = line.split("=")
                    dict_path_info[line[0]] = line[1].replace("\n", "")
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except Exception as error:
        print("An error has occurred: " + str(error))

    if len(dict_path_info) == 0:
        print("The info-file does not contain any information associated with "
              "create_matrix_data.py. ")
        exit()

    return dict_path_info


def ReadFam(dict_path_info):
    """
    This function reads the file containing information about the families and
    extracts, among other things, the gender, age. This is then saved in the
    file "Covariates.txt" which is created. Also, the headers for the file
    "SNP.txt" and "GE.txt" are created and stored in a list.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return list_ge: is a list containing the header consisting of all barcodes
    of the individuals, as: ["id", "BD_305306", "BD_305308", ...]
    :return list_snp: is a list containing the header consisting of all barcodes
    of the individuals, as: ["id", "BD_305306", "BD_305308", ...]
    """
    index_age, index_gender, index_barcode = "", "", ""
    list_age = []
    list_gender = []
    list_covariates = []
    list_ge = []
    list_snp = []

    try:
        list_ge.append("id")
        list_snp.append("id")
        file_covariates = open(dict_path_info['output'] + "Covariates.txt", "w")
        list_covariates.append("id")
        list_age.append("age")
        list_gender.append("gender")

        with open(dict_path_info['info_fam'], "r") as file:
            for line in file:
                line = line.replace("\n", "")
                split = line.split("\t")
                if split.__contains__("Role"):
                    for index in range(len(split)):
                        if split[index].__contains__("Role"):
                            index_age = index
                        if split[index].__contains__("Mol_sex"):
                            index_gender = index
                        if split[index].__contains__("RNR"):
                            index_barcode = index
                elif index_age != "" and index_gender != "" and index_barcode != "":
                    CreateCovariates(split, index_age, index_gender,
                                             list_gender, list_age)
                    list_covariates.append(split[index_barcode])
                    list_snp.append(split[index_barcode])
                    list_ge.append(split[index_barcode])

        file_covariates.write("\t".join(list_covariates))
        file_covariates.write("\n")
        file_covariates.write("\t".join(list_age))
        file_covariates.write("\n")
        file_covariates.write("\t".join(list_gender))
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return list_ge, list_snp


def CreateCovariates(split, index_age, index_gender, list_gender, list_age):
    """
    This function creates two lists containing header order, age and age for
    all individuals.
    :param split: is a list containing the information of a line from the file
    containing general information about the individuals, as:
    ["chick", "BD_305306", "female", ...]
    :param index_age: is a sting which contains the index number for age
    :param index_gender: is a sting which contains the index number for gender
    :param list_gender: is a list containing the gender for all individuals in
    header order, as:  ["0", "0", "1", ...]
    :param list_age: is a list containing age for all individuals in header
    order, as: ["0", "0", "1", ...]
    :return: -
    """
    try:
        if split[index_age] == "Chick":
            list_age.append("0")
        elif split[index_age] == "Parent":
            list_age.append("1")
        gender = split[index_gender].split(" ")
        if gender[2] == "Female":
            list_gender.append("0")
        elif gender[2] == "Male":
            list_gender.append("1")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


def ReadNames(dict_path_info):
    """
    This function retrieves the RefSeq name and the GenBank name for all
    chromosomes and scaffold and stores them in a dictionary.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return dict_names: is a dictionary containing as key the GenBank name and
    as value the RefSeq name for all chromosomes and scaffols, as:
    {'NC_031768.1': 'chr1', 'NC_031769.1': 'chr2', ......}
    """
    index_refseq, index_genbank = "", ""
    dict_names = {}
    try:
        with open(dict_path_info['info_chr_scaf'], "r") as file:
            for line in file:
                if line.__contains__("RefSeq"):
                    line = line.split("\t")
                    for x in range(len(line)):
                        if line[x].__contains__("GenBank"):
                            index_genbank = x
                        elif line[x].__contains__("RefSeq"):
                            index_refseq = x
                else:
                    line = line.split("\t")
                    dict_names[line[index_genbank]] = line[index_refseq]
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return dict_names


def ReadSNP(dict_path_info, list_snp):
    """
    This function reads the file containing the genotype per SNP and per
    individual. This information is then stored in a dictionary. A list is
    also created indicating the exact location of the SNP.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :param list_snp: is a list containing the header consisting of all barcodes
    of the individuals, as: ["id", "BD_305306", "BD_305308", ...]
    :return list_loc: is a list containing all SNP locations, as:
    ["NC_031768.1_2243", "NC_031768.1_56565", ...]
    :return dict_snp_info: is a dictionary containing as key the barcode of an
    individual and as value a list containing all genotypes for all SNPs, as:
    {'BD_305306': ['0/1', '0/0'], 'BD_305308': ['0/1', '0/0'], ...}
    """
    list_loc = []
    dict_snp_info = {}

    try:
        with open(dict_path_info['snp_data'], 'r') as file:
            first_line_info = file.readline()
            first_line_info = first_line_info.replace("\n", "")
            split_info = first_line_info.split("\t")
            split_info2 = []
            for i in split_info:
                if i.__contains__("_mergedsorted"):
                    split = i.split("_")
                    split_info2.append(split[0] + "_" + split[1])
                else:
                    split_info2.append(i)

            for line in file:
                line = line.split("\t")
                loc = line[0] + "_" + line[1]
                list_loc.append(loc)

                for i in range(1, len(list_snp)):
                    index = split_info2.index(list_snp[i])
                    if list_snp[i] in dict_snp_info.keys():
                        dict_snp_info[list_snp[i]] += [line[index]]
                    else:
                        dict_snp_info[list_snp[i]] = [line[index]]
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return list_loc, dict_snp_info


def AddSnps(list_loc, list_snp, dict_snp_info, dict_path_info):
    """
    This function creates the file "SNP.txt" and appends the header to it.
    Next, one SNP locations per line are listed, indicating the genotype
    (0, 1 or 2) for each individual.
    :param list_loc: is a list containing all SNP locations, as:
    ["NC_031768.1_2243", "NC_031768.1_56565", ...]
    :param list_snp: is a list containing the header consisting of all barcodes
    of the individuals, as: ["id", "BD_305306", "BD_305308", ...]
    :param dict_snp_info: is a dictionary containing as key the barcode of an
    individual and as value a list containing all genotypes for all SNPs, as:
    {'BD_305306': ['0/1', '0/0'], 'BD_305308': ['0/1', '0/0'], ...}
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: -
    """
    try:
        file = open(dict_path_info['output'] + "SNP.txt", "w")
        file.write("\t".join(list_snp))
        file.write("\n")
        length_bar = len(list_snp) - 1
        for i in range(len(list_loc)):
            file.write(list_loc[i] + "\t")
            count = 1
            for k in range(1, len(list_snp)): #split bevat de barcodes + id
                if count != length_bar:
                    addition = "\t"
                else:
                    addition = ""
                if dict_snp_info[list_snp[k]][i] == "0/0":
                    file.write("0" + addition)
                elif dict_snp_info[list_snp[k]][i] == "0/1" or \
                        dict_snp_info[list_snp[k]][i] \
                        == "1/0" or dict_snp_info[list_snp[k]][i] == "." or \
                        dict_snp_info[list_snp[k]][i] == "./." or \
                        dict_snp_info[list_snp[k]][i] == "./0" or \
                        dict_snp_info[list_snp[k]][i] == "./1" or \
                        dict_snp_info[list_snp[k]][i] == "0/." or \
                        dict_snp_info[list_snp[k]][i] == "1/.":
                    file.write("1" + addition)
                elif dict_snp_info[list_snp[k]][i] == "1/1":
                    file.write("2" + addition)
                count += 1
            file.write("\n")
        file.close()
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


def CreateSNPloc(list_loc, dict_names, dict_path_info):
    """
    This function creates the file "snpsloc.txt" in which an SNP is added per
    line and on which chromosome it is located and then the exact location.
    :param list_loc: is a list containing all SNP locations, as:
    ["NC_031768.1_2243", "NC_031768.1_56565", ...]
    :param dict_names: is a dictionary containing as key the GenBank name and
    as value the RefSeq name for all chromosomes and scaffols, as:
    {'NC_031768.1': 'chr1', 'NC_031769.1': 'chr2', ......}
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: -
    """
    try:
        file = open(dict_path_info['output'] + "snpsloc.txt", "w")
        file.write("snp" + "\t" + "chr" + "\t" + "pos" + "\n")

        for i in list_loc:
            file.write(i + "\t")
            split = i.split("_")
            chr_saf = split[0] + "_" + split[1]
            if chr_saf in dict_names.keys():
                file.write(dict_names[chr_saf] + "\t")
                file.write(split[2] + "\n")
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


def ReadPos(dict_path_info):
    """
    This function reads the file in which location contains information on
    the CpG sites. This information is then stored in a list.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return list_pos: is a list indicating the locations of the CpG sites,
    such as: ["NC_031768.1_2243", "NC_031768.1_56565", ...]
    """
    list_pos = []
    try:
        with open(dict_path_info['pos_data'], "r") as file:
            first_line = file.readline().replace('"', '').split(",")
            index_chr = first_line.index("chr")
            index_start = first_line.index("start")
            for line in file:
                if index_chr != "" and index_start != "":
                    line = line.replace('"', '').split(",")

                    list_pos.append(line[index_chr] + "_" + line[index_start])
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))
    return list_pos


def ReadMethFile(file_paths):
    """
    This function reads the files containing the methylation data. This data
    is then stored in a dictionary.
    :param file_paths: is a list containing the pathways and file names to be
    read, as: ["/home/dfchickCalcSig.txt", "/hom/dfparCalcSig.txt"]
    :return meth_dict: is a dictionary which contains as key the exact location
    of the CpG site and also which individual is involved. The value contains
    the methylation percentage, as:
    {'NC_031768.1_2243, BD_305306': '60.01', 'NC_031768.1_56565,
    BD_305308': '80.90', ...}
    """
    meth_dict = {}
    index_pos, index_meth, index_bar = "", "", ""

    try:
        for file in file_paths:
            with open(file, "r") as file:
                for line in file:
                    line = line.replace("\n", "")
                    if line.__contains__("CHR"):
                        line = line.split(",")
                        for x in range(len(line)):
                            if line[x].__contains__('"CHR_POS"'):
                                index_pos = x
                            if line[x].__contains__('"Methylation"'):
                                index_meth = x
                            if line[x].__contains__('"RNR"'):
                                index_bar = x
                    elif not line.__contains__("CHR"):
                        line = line.split(",")
                        if index_pos != "" and index_meth != "" and index_bar != "":
                            key = line[index_pos] + ", " + line[index_bar]
                            meth_dict[key] = line[index_meth]
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return meth_dict


def AddGe(list_pos, list_ge, meth_dict, dict_path_info):
    """
    This function creates the "GE.txt" file. Then the header is added to the
    file. Hereafter, the methylation percentage is indicated per CpG site and
    per individual.
    :param list_pos: is a list indicating the locations of the CpG sites,
    such as: ["NC_031768.1_2243", "NC_031768.1_56565", ...]
    :param list_ge: is a list containing the header consisting of all barcodes
    of the individuals, as: ["id", "BD_305306", "BD_305308", ...]
    :param meth_dict: is a dictionary which contains as key the exact location
    of the CpG site and also which individual is involved. The value contains
    the methylation percentage, as:
    {'NC_031768.1_2243, BD_305306': '60.01', 'NC_031768.1_56565,
    BD_305308': '80.90', ...}
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: -
    """
    try:
        file = open(dict_path_info['output'] + "GE.txt", "w")
        file.write("\t".join(list_ge))
        file.write("\n")
        length_bar = len(list_ge) - 1

        for meth in list_pos:
            file.write(meth + "\t")
            list_gene = []
            for bar in range(1, len(list_ge)):
                key = meth + ', ' + list_ge[bar]
                if key in meth_dict.keys():
                    if meth_dict[key] != "NA":
                        ge = meth_dict[key]
                        list_gene.append(ge)
                    else:
                        ge = meth_dict[key]
                        list_gene.append(ge)
                else:
                    list_gene.append("NA")
            for el in range(len(list_gene)):
                if el != length_bar - 1:
                    file.write(list_gene[el] + "\t")
                else:
                    file.write(list_gene[el])
            file.write("\n")
        file.close()
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))


def CreateGeneloc(list_pos, dict_names, dict_path_info):
    """
    This function creates the file "geneloc.txt" and per line it is indicated
    which CpG site it concerns and on which chromosome it is located. Then the
    exact start and end location is also indicated.
    :param list_pos: is a list indicating the locations of the CpG sites,
    such as: ["NC_031768.1_2243", "NC_031768.1_56565", ...]
    :param dict_names: is a dictionary containing as key the GenBank name and
    as value the RefSeq name for all chromosomes and scaffols, as:
    {'NC_031768.1': 'chr1', 'NC_031769.1': 'chr2', ......}
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
   {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: -
    """
    try:
        file = open(dict_path_info['output'] + "geneloc.txt", "w")
        file.write("geneid" + "\t" + "chr" + "\t" + "s1" + "\t" + "s2" + "\n")

        for i in list_pos:
            file.write(i + "\t")
            split = i.split("_")
            chr_saf = split[0] + "_" + split[1]
            if chr_saf in dict_names.keys():
                file.write(dict_names[chr_saf] + "\t")
                file.write(split[2] + "\t")
                file.write(split[2] + "\n")
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
    list_ge, list_snp = ReadFam(dict_path_info)
    dict_names = ReadNames(dict_path_info)

    #snp_files
    list_loc, dict_snp_info = ReadSNP(dict_path_info, list_snp)
    AddSnps(list_loc, list_snp, dict_snp_info, dict_path_info)
    CreateSNPloc(list_loc, dict_names, dict_path_info)

    #ge_files
    list_pos = ReadPos(dict_path_info)
    file_paths = [dict_path_info['meth_chick_data'], dict_path_info['meth_par_data']]
    meth_dict = ReadMethFile(file_paths)
    AddGe(list_pos, list_ge, meth_dict, dict_path_info)
    CreateGeneloc(list_pos, dict_names, dict_path_info)
