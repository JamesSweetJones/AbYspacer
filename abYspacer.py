#!/usr/bin/python3

import argparse
import sys
import re
import requests

class annotate():

    def __init__(self, aaseq, scheme):

        self.aaseq=aaseq
        self.scheme=scheme

    def __repr__(self):
        return "Annotation of VH or VL sequence using Kabat, Chothia, Contact, or IMGT scheme"
    def output(self, chain, lst, regionlst):
        """
        Prints the FR and CDR regions and their corresponding seq. It returns a `list` of 2 `dict`.
        :param chain: STRING, either "H" or "L" in uppercase
        :param lst:  LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :param regionlst: LIST, a list of peptides, each corresponds to a FR or CDR region
        :return: LIST, a list of 2 `dict`, The first dict consists of region: seq pairs. The second dict consists of number:residue pairs.
        """
        self.chain=chain
        self.lst=lst
        self.regionlst=regionlst

        self.regiondict, self.numberdict={}, {}

        for i in range (0, len(self.lst), 2):
            self.numberdict[self.lst[i]]=self.lst[i+1]


        if self.chain=="L":
            Region_list = [self.regionlst[0],self.regionlst[1],self.regionlst[2],self.regionlst[3],self.regionlst[4],self.regionlst[5],self.regionlst[6]]
            #print(Region_list)

            for region, seq in zip(["L-FR1", "L-CDR1", "L-FR2","L-CDR2", "L-FR3", "L-CDR3", "L-FR4"], self.regionlst):
                self.regiondict[region]=seq

            return [self.regionlst, self.numberdict]

        else:
            Region_list = [self.regionlst[0],self.regionlst[1],self.regionlst[2],self.regionlst[3],self.regionlst[4],self.regionlst[5],self.regionlst[6]]
            #print(Region_list)

            for region, seq in zip(["H-FR1", "H-CDR1", "H-FR2","H-CDR2", "H-FR3", "H-CDR3", "H-FR4"], self.regionlst):
                self.regiondict[region]=seq
            return [self.regionlst, self.numberdict]


    def analyze(self,chain, lst):
        """
        Define CDR and FR regions based on the numbered sequence returned from website
        :param chain: STRING, "H" or "L" in uppercase
        :param lst: LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :return: LIST, a list of strings, where each string is a peptide corresponding to the a region, in the order of: FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4
        :raises: `ValueError` if any of the FR or CDR region is missing
        """

        self.chain=chain
        self.lst=lst
        if self.chain=="L":
            self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4=["" for i in range (0, 7)]

            try:
                if self.scheme in ["kabat", "chothia"]:
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L24"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L24"), self.lst.index("L35"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L35"), self.lst.index("L50"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L50"), self.lst.index("L57"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L57"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L98"), len(self.lst), 2)])

                elif self.scheme =="contact":
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L30"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L30"), self.lst.index("L37"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L37"), self.lst.index("L46"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L46"), self.lst.index("L56"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L56"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L97"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L97"), len(self.lst), 2)])

                elif self.scheme =="martin":
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L24"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L24"), self.lst.index("L34"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L34"), self.lst.index("L50"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L50"), self.lst.index("L56"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L56"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L97"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L97"), len(self.lst), 2)])

                else: #IMGT scheme
                    self.L_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("L27"), 2)])
                    self.L_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("L27"), self.lst.index("L33"), 2)])
                    self.L_FR2="".join([self.lst[i+1] for i in range (self.lst.index("L33"), self.lst.index("L50"), 2)])
                    self.L_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("L50"), self.lst.index("L52"), 2)])
                    self.L_FR3="".join([self.lst[i+1] for i in range (self.lst.index("L52"), self.lst.index("L89"), 2)])
                    self.L_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4="".join([self.lst[i+1] for i in range (self.lst.index("L98"), len(self.lst), 2)])

                return [self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4]

            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured")
        else:
            self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4=["" for i in range (0, 7)]

            try:
                if self.scheme=="kabat":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H31"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H31"), self.lst.index("H36"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H36"), self.lst.index("H50"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H50"), self.lst.index("H66"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H66"), self.lst.index("H95"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H103"), len(self.lst), 2)])

                elif self.scheme=="chothia":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H26"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H26"), self.lst.index("H33"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H33"), self.lst.index("H52"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H52"), self.lst.index("H57"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H57"), self.lst.index("H95"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H103"), len(self.lst), 2)])

                elif self.scheme=="contact":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H30"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H30"), self.lst.index("H36"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H36"), self.lst.index("H47"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H47"), self.lst.index("H59"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H59"), self.lst.index("H93"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H93"), self.lst.index("H102"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H102"), len(self.lst), 2)])

                elif self.scheme=="martin":
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H26"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H31"), self.lst.index("H35"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H35"), self.lst.index("H50"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H50"), self.lst.index("H65"), 2)])
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H65"), self.lst.index("H95"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H95"), self.lst.index("H102"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H102"), len(self.lst), 2)])

                else: #IMGT scheme
                    self.H_FR1="".join([self.lst[i+1] for i in range (0, self.lst.index("H26"), 2)])
                    self.H_CDR1="".join([self.lst[i+1] for i in range (self.lst.index("H26"), self.lst.index("H34"), 2)])
                    self.H_FR2="".join([self.lst[i+1] for i in range (self.lst.index("H34"), self.lst.index("H51"), 2)])
                    self.H_CDR2="".join([self.lst[i+1] for i in range (self.lst.index("H51"), self.lst.index("H58"), 2)]) #51>57 (instead of 56)
                    self.H_FR3="".join([self.lst[i+1] for i in range (self.lst.index("H58"), self.lst.index("H93"), 2)])
                    self.H_CDR3="".join([self.lst[i+1] for i in range (self.lst.index("H93"), self.lst.index("H103"), 2)])
                    self.H_FR4="".join([self.lst[i+1] for i in range (self.lst.index("H103"), len(self.lst), 2)])

                return [self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4]

            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured in the `analyze()` method")

    def retrieve (self):

        self.url="http://www.bioinf.org.uk/abs/abnum/abnum.cgi"

        try:
            if self.scheme not in ["kabat", "chothia", "contact", "imgt", "martin"]:
                raise Exception

        except ValueError:
            print("Incorrect scheme mode. Must be one of the following (lowercase): kabat, chothia, contact, imgt")

        else:
            if self.scheme=="kabat":
                self.sche="-k"
            elif self.scheme=="chothia":
                self.sche="-c"
            elif self.scheme=="martin":
                self.sche="-m"

        try:
            self.d={"plain":1, "scheme":self.sche, "aaseq":self.aaseq}
            self.myPage=requests.get(self.url, params=self.d)
            self.text=self.myPage.text
            self.lst=self.text.split()

            if len(self.lst)>1:
                self.chain=self.lst[0][0]
                self.result=self.output(self.chain, self.lst, self.analyze(self.chain, self.lst))
                return self.result
            else:
                print("No annotation retrieved. Did you enter the complete VH or VL sequence?")
        except:
            print("An error occured in the `retrieve()` method")

def get_spaced_sequence(Heavy_seq, Light_seq, scheme):
    try:
        seq1=annotate(Heavy_seq,scheme)
        heavy_sequence_split = seq1.retrieve()[0]
        heavy_sequence_num   = seq1.retrieve()[1]
        heavy_sequence_num_list = list(heavy_sequence_num.keys())
        seq2=annotate(Light_seq,scheme)
        light_sequence_split = seq2.retrieve()[0]
        light_sequence_num   = seq2.retrieve()[1]
        light_sequence_num_list = list(light_sequence_num.keys())


        if scheme == 'kabat':
            Heavy_sorter = ['H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12','H13','H14','H15','H16','H17','H18','H19','H20','H21','H22','H23','H24','H25','H26','H27','H28','H29','H30','H31','H32','H33','H34','H35','H35A','H35B','H36','H37','H38','H39','H40','H41','H42','H43','H44','H45','H46','H47','H48','H49','H50','H51','H52','H52A','H52B','H52C','H53','H54','H55','H56','H57','H58','H59','H60','H61','H62','H63','H64','H65','H66','H67','H68','H69','H70','H71','H72','H73','H74','H75','H76','H77','H78','H79','H80','H81','H82','H82A','H82B','H82C','H83','H84','H85','H86','H87','H88','H89','H90','H91','H92','H93','H94','H95','H96','H97','H98','H99','H100','H100A','H100B','H100C','H100D','H100E','H100F','H100G','H100H','H100I','H100J','H100K','H101','H102','H103','H104','H105','H106','H107','H108','H109','H110','H111','H112','H113']
            Light_sorter = ['L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16','L17','L18','L19','L20','L21','L22','L23','L24','L25','L26','L27','L27A','L27B','L27C','L27D','L27E','L27F','L28','L29','L30','L31','L32','L33','L34','L35','L36','L37','L38','L39','L40','L41','L42','L43','L44','L45','L46','L47','L48','L49','L50','L51','L52','L53','L54','L55','L56','L57','L58','L59','L60','L61','L62','L63','L64','L65','L66','L67','L68','L69','L70','L71''L72','L73','L74','L75','L76','L77','L78','L79','L80','L81','L82','L83','L84','L85','L86','L87','L88','L89','L90','L91','L92','L93','L94','L95','L95A','L95B','L95C','L95D','L95E','L95F','L96','L97','L98','L99','L100','L101','L102','L103','L104','L105','L106','L106A','L107','L108','L109']
        elif scheme == 'chothia':
            Heavy_sorter = ['H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12','H13','H14','H15','H16','H17','H18','H19','H20','H21','H22','H23','H24','H25','H26','H27','H28','H29','H30','H31','H31A','H31B','H32','H33','H34','H35','H36','H37','H38','H39','H40','H41','H42','H43','H44','H45','H46','H47','H48','H49','H50','H51','H52','H52A','H52B','H52C','H53','H54','H55','H56','H57','H58','H59','H60','H61','H62','H63','H64','H65','H66','H67','H68','H69','H70','H71','H72','H73','H74','H75','H76','H77','H78','H79','H80','H81','H82','H82A','H82B','H82C','H83','H84','H85','H86','H87','H88','H89','H90','H91','H92','H93','H94','H95','H96','H97','H98','H99','H100','H100A','H100B','H100C','H100D','H100E','H100F','H100G','H100H','H100I','H100J','H100K','H101','H102','H103','H104','H105','H106','H107','H108','H109','H110','H111','H112','H113']
            Light_sorter = ['L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16','L17','L18','L19','L20','L21','L22','L23','L24','L25','L26','L27','L28','L29','L30','L30A','L30B','L30C','L30D','L30E','L30F','L31','L32','L33','L34','L35','L36','L37','L38','L39','L40','L41','L42','L43','L44','L45','L46','L47','L48','L49','L50','L51','L52','L53','L54','L55','L56','L57','L58','L59','L60','L61','L62','L63','L64','L65','L66','L67','L68','L69','L70','L71','L72','L73','L74','L75','L76','L77','L78','L79','L80','L81','L82','L83','L84','L85','L86','L87','L88','L89','L90','L91','L92','L93','L94','L95','L95A','L95B','L95C','L95D','L95E','L95F','L96','L97','L98','L99','L100','L101','L102','L103','L104','L105','L106','L106A','L107','L108','L109']
        elif scheme == 'martin':
            Heavy_sorter = ['H1','H2','H3','H4','H5','H6','H7','H8','H8A','H8B','H8C','H8D','H9','H10','H11','H12','H13','H14','H15','H16','H17','H18','H19','H20','H21','H22','H23','H24','H25','H26','H27','H28','H29','H30','H31','H31A','H31B','H32','H33','H34','H35','H36','H37','H38','H39','H40','H41','H42','H43','H44','H45','H46','H47','H48','H49','H50','H51','H52','H52A','H52B','H52C','H53','H54','H55','H56','H57','H58','H59','H60','H61','H62','H63','H64','H65','H66','H67','H68','H69','H70','H71','H72','H72A','H72B','H72C','H73','H74','H75','H76','H77','H78','H79','H80','H81','H82','H83','H84','H85','H86','H87','H88','H89','H90','H91','H92','H93','H94','H95','H96','H97','H98','H99','H100','H100A','H100B','H100C','H100D','H100E','H100F','H100G','H100H','H100I','H100J','H100K','H101','H102','H103','H104','H105','H106','H107','H108','H109','H110','H111','H112','H113']
            Light_sorter = ['L0','L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16','L17','L18','L19','L20','L21','L22','L23','L24','L25','L26','L27','L28','L29','L30','L30A','L30B','L30C','L30D','L30E','L30F','L31','L32','L33','L34','L35','L36','L37','L38','L39','L40','L40A','L41','L42','L43','L44','L45','L46','L47','L48','L49','L50','L51','L52','L52A','L52B','L52C','L52D','L52E','L53','L54','L55','L56','L57','L58','L59','L60','L61','L62','L63','L64','L65','L66','L67','L68','L68A','L68B','L68C','L68D','L68E','L68F','L68G','L68H','L69','L70','L71','L72','L73','L74','L75','L76','L77','L78','L79','L80','L81','L82','L83','L84','L85','L86','L87','L88','L89','L90','L91','L92','L93','L94','L95','L95A','L95B','L95C','L95D','L95E','L95F','L96','L97','L98','L99','L100','L101','L102','L103','L104','L105','L106','L107','L107A','L108','L109','L110']
        else:
            print("no scheme provided")
            quit()
        spaced_heavy_seq = ""
        for x in range(len(Heavy_sorter)):
            if Heavy_sorter[x] not in heavy_sequence_num_list:
                spaced_heavy_seq += "X"
            elif Heavy_sorter[x] in heavy_sequence_num_list:
                spaced_heavy_seq += str(heavy_sequence_num.get(Heavy_sorter[x]))

            #reordered_Heavy_dict = {x: Heavy_chain_ptm_counts[x] for x in Heavy_sorter}

        spaced_light_seq = ""
        for x in range(len(Light_sorter)):
            if Light_sorter[x] not in light_sequence_num:
                spaced_light_seq += "X"
            elif Light_sorter[x] in light_sequence_num_list:
                spaced_light_seq += str(light_sequence_num.get(Light_sorter[x]))
        return(spaced_heavy_seq,spaced_light_seq)
    except:
        pass

if __name__ == '__main__':

    my_parser = argparse.ArgumentParser(prog="abYspacer",
                                        usage='python %(prog)s [options]',
                                        description="Programme to generate numbered antibody sequences with 'X' to denote spaces in the numbering scheme")
    my_parser.add_argument('-i','--input', type=str,help='the path to paired fasta file')
    my_parser.add_argument('-o','--output', type=str,help='string of output fasta file')
    my_parser.add_argument('-s','--scheme', type=str,help='k/kabat, c/chothia, m/martin')

        #my_parser.print_help()
    args = my_parser.parse_args()
    input_fasta = args.input
    if input_fasta is None:
        print('No input was given. Exiting programme')
        sys.exit()


    scheme = args.scheme
    if scheme is None or scheme == "kabat" or scheme == "k":
        scheme = 'kabat'
    elif scheme == "chothia" or scheme == "c":
        scheme = 'chothia'
    elif scheme == 'martin' or scheme == "m":
        scheme = 'martin'
    else:
        print("no suitable numbering scheme provided")
        quit()
    output_name = args.output
    if output_name is None:
        output_name = str(input_fasta+"_"+str(scheme)+".faa")

    with open (input_fasta, "r") as f:
        output = open(output_name,"w+")
        for line in f:
            if line[0] == ">":
                if "_H|" in line or "_VH|" in line:
                    Heavy_seq_identifier = line.strip()
                    H_identifier = Heavy_seq_identifier.split("|")[1]
                    Heavy_seq = f.readline().strip()
                    Light_seq_identifier = f.readline().strip()
                    L_identifier = Light_seq_identifier.split("|")[1]
                    Light_seq = f.readline().strip()
                    if H_identifier == L_identifier:
                        spaced_sequences = get_spaced_sequence(Heavy_seq, Light_seq, scheme)
                        if spaced_sequences is not None:
                            output.write(str(">"+H_identifier+"_H|"+H_identifier+"\n"+spaced_sequences[0]+"\n>"+H_identifier+"_L|"+H_identifier+"\n"+spaced_sequences[1]+"\n"))
                elif "_L|" in line or "_VH|" in line:
                    Light_seq_identifier = line.strip()
                    L_identifier = Light_seq_identifier.split("|")[1]
                    Light_seq = f.readline().strip()
                    Heavy_seq_identifier = line.strip()
                    H_identifier = f.readline().split("|")[1]
                    Heavy_seq = f.readline().strip()
                    if H_identifier == L_identifier:
                        get_spaced_sequence(Heavy_seq, Light_seq, scheme)
                        if spaced_sequences is not None:
                            output.write(str(">"+H_identifier+"_H|"+H_identifier+"\n"+spaced_sequences[0]+"\n>"+H_identifier+"_L|"+H_identifier+"\n"+spaced_sequences[1]+"\n"))
