"""
Using Floor function.
"""
from lib2to3.pytree import BasePattern
import os
from copyreg import pickle
from math import ceil
import pickle
from definition_parser import *
import json

BASEPATH = os.getcwd()
MAGNAME_PATH = "assets/tatoosh_MAG_names.txt"
MASTERDATA_PATH = "assets/tatoosh_masterdata.txt"
META_DEF = {
    "Assimilatory_Sulfur_Reduction": Gph(Asr),
    "Dissimilatory_Sulfur_Reduction": Gph(Dsr),
    "Thiosulfate_Oxidation": Gph(TsO),
    "Nitrogen_Fixation": Gph(Nif),
    "Assimilatory_Nitrate_Reduction": Gph(Anr),
    "Dissimilatory_Nitrate_Reduction": Gph(Dnr),
    "Denitrification": Gph(Denitrification),
    "Nitrification": Gph(Nitrification),
    "Comammox": Gph(coammox),
    "Anammox": Gph(anammox),
    "Vit_B1": Gph(thiamin),
    "Vit_B2":Gph(ribo),
    "Vit_B7":Gph(biotin),
    "Vit_B12_Aerobic":Gph(cob_aerobic),
    "Vit_B12_Anaerobic": Gph(cob_anaerobic)
    }

THRESHOLD = 0.25 # fraction allowed to go missing

class MetaDetect:
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
        with open(f"{BASEPATH}/{MAGlist}",'r') as m:
            m = m.readlines()
            for mag in m:
                mag_list[mag.strip()] = []
        return mag_list
    
    def loadMasterdata(self, masterdata):
        print("\n~~~~~  LOADING MASTERDATA  ~~~~~")
        # Takes in string of path of masterdata
        # Only takes in KOFam calls.

        mag_gene_dict = {}
        with open(f"{BASEPATH}/{masterdata}", 'r') as g:
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
                    threshold = math.floor(len(pathway)*THRESHOLD) #ALT
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
                        pathOutput = [mag, str(threshold), str(len(pathway)), gph, str(missing), str(count), "-".join([self.graphs[gph].node2gene[node] for node in pathway])]
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
            filename = f"{BASEPATH}/trail_misc_table_{THRESHOLD}.txt"
        else:
            filename = f"{BASEPATH}/output_{THRESHOLD}.txt"

        with open(filename, 'w') as t:
            t.write("\n".join(self.output))


if __name__=="__main__":
    # # MetaDetect(META_DEF, MAGNAME_PATH, MASTERDATA_PATH)

    with open(f'{BASEPATH}/assets/all_KEGG_defs.json','r') as t:
        data = json.load(t)

    print("\n\n\n\n\n\n\n")

    # pickling Graph db of all KEGG defs
    defDict = dict()
    for category in data.keys():
        catName = category.replace(" ","_")
        for path in data[category]:
            pathDesc = path+"_"+data[category][path]['name'].replace(" ","_")

            for module in data[category][path]['modules']:
                try:
                    modDesc = module + "_" + data[category][path]['modules'][module]['name'].replace(" ","_")

                    defDict[f"{catName}__{pathDesc}__{modDesc}"] = Gph(data[category][path]['modules'][module]['definition'])
                except:
                    print("YOOOO\n",module, data[category][path]['modules'][module])
                    pass # for Modules with empty definitions
    
    print("\t".join(defDict.keys()))

    # with open(f"{BASEPATH}/assets/KEGG_def_pickle", "wb") as def_pickle:
    #     pickle.dump(defDict, def_pickle)

    ###########
    # # Loads a dictionary of "Module name": Gph(KEGG Def) that has been prepickled
    # # and runs the Meta Detect Algorithm on it.
    with open(f"{BASEPATH}/assets/KEGG_def_pickle", "rb") as u:
        new_dict = pickle.load(u)

    MetaDetect(new_dict, MAGNAME_PATH, MASTERDATA_PATH)
