"""
Using Floor function.
"""

from copyreg import pickle
from math import ceil
import pickle
from definition_parser import *
import json
from fig3KEGGdefinitions import *
import argparse, sys



class NAMeD:
    def __init__(self,  metabDef: dict, MAGlist, masterdata, misclist=True) -> None:
        self.misc = misclist # generate misc table (TRUE) or output detail (FALSE)
        self.MAGNames = self.loadMagNames(MAGlist)

        self.mData = self.loadMasterdata(masterdata)

        self.graphs = metabDef
        
        self.output = []
        self.functionDetect()

        self.printOutput()


    def loadMagNames(self, MAGlist) -> dict:
        print("\n~~~~~  LOADING MAG Names  ~~~~~")
        # Takes in string of path of MAG names
        mag_list = {}
        with open(f"{MAGlist}",'r') as m:
            m = m.readlines()
            for mag in m:
                mag_list[mag.strip()] = []
        return mag_list
    
    def loadMasterdata(self, masterdata):
        print("\n~~~~~  LOADING MASTERDATA  ~~~~~")
        # Takes in string of path of masterdata
        # Only takes in KOFam calls.

        mag_gene_dict = {}
        with open(f"{masterdata}", 'r') as g:
            g = g.readlines()
            for line in g:
                line = line.split()

                if float(line[-2]) < 1e-20: # KOFam hits under 1e-50
                    if line[3] == "KOfam":
                        if line[2] in mag_gene_dict:
                            if "!!!" in line[4]:
                                for gene in line[4].split("!!!"):
                                    mag_gene_dict[line[2]].append(gene)
                            else:
                                mag_gene_dict[line[2]].append(line[4])
                        else: # First time mag name appears
                            if "!!!" in line[4]:
                                mag_gene_dict[line[2]] = line[4].split("!!!")
                            else:
                                mag_gene_dict[line[2]] = [line[4]]

        return mag_gene_dict
    
    def functionDetect(self):
        print("\n~~~~~  DETECTING METABOLISM  ~~~~~")

        for mag in self.MAGNames:
            if self.misc == True:
                magOutput = [mag]

            for gph in self.graphs:
                for pathway in self.graphs[gph].listOfPathways:
                    ignore = False
                    pathOutput = []
                    pwlen = len(pathway)
                    
                    threshold = math.floor(pwlen*THRESHOLD)
                    if pwlen < math.floor(1/THRESHOLD)+1:
                        threshold+=1 # Gives smaller pathways a fighting chance

                    missing = 0
                    count = 0
                    for node in pathway:
                        count += 1
                        if self.graphs[gph].node2gene[node] in self.mData[mag]:
                            pass
                        else:
                            # print(node_to_gene[node])
                            missing += 1
                            if missing == threshold+1:
                                ignore = True
                                
                    
                    if (ignore == False):
                        pathOutput = [mag, str(threshold), str(pwlen), gph, str(missing), str(count), "-".join([self.graphs[gph].node2gene[node] for node in pathway])]
                        if self.misc == False:
                            self.output.append("\t".join(pathOutput))
                        break

                if self.misc==True:
                    if ignore == False:
                        magOutput.append("Y")
                    else:
                        magOutput.append("N")

            if self.misc==True:
                self.output.append("\t".join(magOutput))    
                
    def printOutput(self):
        print("\n~~~~~  PRINTING OUTPUT  ~~~~~")

        if self.misc == False:
            header = "MAG\tTHRESHOLD\tPWLEN\tMETABOLISM\tMISSING\tGENECOUNT\tPATHWAY"
        else:
            header = "\t".join(['MAG_names']+[df for df in self.graphs])
        self.output.insert(0,header)

        if self.misc == True:
            filename = f"misc_table_{THRESHOLD}_mod.txt"
        else:
            filename = f"outputLog_{THRESHOLD}.txt"

        with open(filename, 'w') as t:
            t.write("\n".join(self.output))


if __name__=="__main__":
    helpText = "Network Algorithm for Metabolic Detection (NAMeD)\n\nThis algorithm takes a masterdata sheet of genes from your genomes and detects from a list of metabolisms whether your genomes have any of them present. It takes the KEGG definition of each metabolism and represents it as a network of nodes - with starting and target nodes. If it's possible to get from the starting to the target node of your genome, then the metabolism is present. Make sure to make use of the -t flag to control the percentage of the pathway you'll allow missing."

    # DEFAULTS
    MAGNAME_PATH = "./assets/tatoosh_MAG_names.txt"

    # Parsing args
    parser = argparse.ArgumentParser(description=helpText)
    parser.add_argument("-i","--input",required=True)
    parser.add_argument("-t","--threshold",default=0.33,type=float)
    parser.add_argument("-m", "--mode", default='single', choices=['single','multi','all'])
    parser.add_argument("-d", "--definition")
    args = parser.parse_args()
    doIt = True

    # SELECTING INPUT FILE
    if args.input:
        MASTERDATA_PATH = args.input
    else:
        print("INPUT FILE ERROR: Check input file path")

    # SETTING DICTIONARY OF METABOLISMS
    if args.definition: # not default
        if args.mode=='single':
            try: 
                META_DEF = {"Def1": Gph(args.definition)}
            except:
                print("Input Error: Verify input string is a KEGG definition")
                doIt = False
        elif args.mode=='all':
            with open(f"assets/KEGG_def_pickle", "rb") as u:
                META_DEF = pickle.load(u)      

            print("\n\nThis might take a while\n\n") 

    else:               # default
        if args.mode:
            if args.mode=='single':
                META_DEF = {"Denitrification": Gph(Denitrification)}
            elif args.mode=='multi':
                META_DEF = defaultMultiDef
            else: # args.mode=='all'
                # Loads a dictionary of "Module name": Gph(KEGG Def) that has been prepickled
                # and runs the Meta Detect Algorithm on it.
                with open(f"assets/KEGG_def_pickle", "rb") as u:
                    META_DEF = pickle.load(u)      

                print("\n\nThis might take a while\n\n") 
        else:
            META_DEF = {"Denitrification": Gph(Denitrification)}

    # SETTING THRESHOLD

    if args.threshold and args.threshold<=1 and args.threshold >=0:
        THRESHOLD=args.threshold
    else:
        THRESHOLD = 0.33 # DEFAULT

    print(f"\nUSING {THRESHOLD} AS THRESHOLD")
        
    # RUNNING CODE
    if doIt == True:
        NAMeD(META_DEF, MAGNAME_PATH, MASTERDATA_PATH)
    else:
        pass