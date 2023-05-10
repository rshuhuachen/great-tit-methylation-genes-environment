"""
Author: Michelle Memelink
Date: 24-01-2022
Description: The file containing the barcodes is read in, provided the column
(Barcode_R2) is present, these barcodes are temporarily stored in a list.
Subsequently, the input file is searched for per barcode and an output file is
defined. The command line is called by means of a subprocess. The command
assumes that PICARD (from Conda) is installed and you are running the script in
this environment. The command causes each file to receive a unique RG-tag.
This makes it clear which data belongs to which individual.
"""
try:
    import subprocess
except ImportError:
    raise ImportError("An error has occurred when importing a module")


def ReadInfo(path_info):
    """
    the function reads in the info file and extracts only the information that
    applies to the python script in which this function is run (add_tags.py).
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
                if line.__contains__("add_tags.py"):
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
              "add_tags.py.")
        exit()

    return dict_path_info


def OpenFile(dict_path_info):
    """
    By means of this function, the barcodes of the organisms are obtained.
    These barcodes are placed in a list "list_barcodes", which is then returned.
    First, the column "Barcode_R2" is searched for, if it is present, the
    index of this column is stored in a string. then all barcodes are added to
    the list by means of the index. If the file does not contain the column
    containing 'Barcode_R2', the script will be stopped.
    :param dict_path_info: is a dictionary containing all input and output file
     names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return list_barcodes: is a list containing the barcodes of the organisms
    from which a file is to be obtained with RG-tag.
    """
    pathways = [dict_path_info['barcode1'], dict_path_info['barcode2']]
    list_barcodes = []
    index_barcode = ""
    try:
        for path in pathways:
            with open(path, "r") as file:
                first_line = file.readline().split("\t")
                for x in range(len(first_line)):
                    if first_line[x].__contains__("Barcode_R2"):
                        index_barcode = int(x) + 1
                for line in file:
                    if index_barcode != "":
                        line = line.split("\t")
                        list_barcodes.append(line[index_barcode])
                    else:
                        print("File is not in correct format")
                        exit()
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the correct "
              "information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return list_barcodes


def AddTags(list_barcodes, dict_path_info):
    """
    this function is able to define the input file location and name for each
    barcode from the input list. Also, the output location and file name are
    defined for each barcode. After defining an input and output files for a
    barcode, a subprocess is called. The command causes an RG tag to be added
    to the output file. This is done by means of the "PICARD" tool.
    :param list_barcodes: is a list containing the barcodes of the organisms
    from which a file is to be obtained with RG tag.
    :param dict_path_info: is a dictionary containing all input and output file
    names and paths described. with the format:
    {'barcode1': '/home/epiGBS2/data/barcode1_4.tsv', 'output': '/home/epiGBS2/'}
    :return: nothing
    """
    for code in list_barcodes:
        input = dict_path_info['input'] + code + \
                "_trimmed_filt_merged.1_bismark_bt2_pe.bam"
        output = dict_path_info['output'] + code + \
                 "_trimmed_filt_merged.1_bismark_bt2_pe_tag.bam"
        try:
            subprocess.run(
                ["picard", "AddOrReplaceReadGroups", "I=" + input,
                 "O=" + output,
                 "RGID=4", "RGLB=lib1", "RGPL=epiGBS2", "RGPU=" + code,
                 "RGSM=" + code])
        except FileNotFoundError:
            print("An error has occurred: wrong file of file path")
        except subprocess.CalledProcessError as error_pro:
            print("An error has occurred: " + error_pro.output)
        except Exception as error:
            print("An error has occurred: " + str(error))


if __name__ == '__main__':
    #Change the pathway below to the correct custom pathway:
    path_info = "/home/info.txt"

    dict_path_info = ReadInfo(path_info)
    list_barcodes = OpenFile(dict_path_info)
    AddTags(list_barcodes, dict_path_info)
