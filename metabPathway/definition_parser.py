"""
This parser takes in a KEGG definition and returns a metabolic network.
The object class Gph both generates the metabolic note and outputs a 
list of all possible paths from the source vertex to the terminal vertex.

As this has been tested on the KEGG definitions only in this package, 
I would recommend verifying the network by plotting (unhash the plot(g) 
line in genGraph() in the Gph Object). But when running, hash out the plot(g) 
function because it really increases runtime

Python 3.9.9
Author: Khashiff Miranda (kkmiranda.github.io)
"""


import re
from platform import node
from igraph import *
from KEGG_definitions import *


class Block:
    def __init__(self, input) -> None:
        self.input = input

        self.genes = set()
        self.edges = []
        self.srcNodes = []
        self.tgtNodes = []
        self.findEdges(input)

    def __str__(self) -> str:
        return "".join(self.input)

    def findEdges(self, input):
        b = 0
        i = 0
        tgt = True #flag that indicates whether a target has been found
        while i < len(input):
            if len(input[i]) == 6: #gene
                call = input[i]
                if call not in self.genes:
                    self.srcNodes.append(call)
                    
                self.genes.add(call)
                j = i + 1
                if j >= len(input):
                    self.tgtNodes.append(call)
                    break
                #find target (' '/'+') or nah (anything else)
                if input[j] == ' ':
                    tgt = False # edge found, is not a target node
                    j += 1 # move target pointer
                    if input[j] == "(":
                        j += 1
                        while input[j]!=')':
                            if len(input[j])==6:
                                self.edges.append((call,input[j]))
                                self.genes.add(input[j])
                            j += 1
                        
                    elif len(input[j]) == 6:
                        self.edges.append((call, input[j]))
                        self.genes.add(input[j])

                elif input[j] == '+':
                    tgt = False # edge found, is not a target node
                    j += 1
                    
                    self.edges.append((call,input[j]))
                    self.genes.add(input[j])
                    
                elif (input[j] == "," or input[j] == ")"):
                    try:
                        if b > 0 and j<len(input)-3:
                            while input[j] != " ":
                                j += 1
                            else:
                                j += 1
                                while input[j]!=')' and input[j]!=",":
                                    if len(input[j])==6:
                                        tgt = False # edge found, is not a target node
                                        self.edges.append((call,input[j]))
                                        self.genes.add(input[j])
                                    j += 1
                                    
                        else:
                            self.tgtNodes.append(call)
                    except:
                        # no clue why this works for thiamin and cobalamin
                        pass
                if tgt == True:
                    self.tgtNodes.append(call)

            else:
                if input[i] == "(":
                    b += 1
                if input[i] == ")":
                    b -= 1
            
            tgt = True
            i+=1 #move call pointer

class Pathway:
    def __init__(self, KEGGdef) -> None:
        self.defString = re.sub(r"\-K\d{5}","",KEGGdef)
        self.tokens = []
        self.tokenize()
        self.blox = []
        self.blockify()

        self.edges = set()
        self.genes = set()
        self.sourceNodes = []
        self.targetNodes = []
        self.unifyBlocks()

    # PARSING INFO
    def tokenize(self):
        i = 0
        while i < len(self.defString):
            if self.defString[i] == "(" or self.defString[i] == "+" or self.defString[i] == ")" or self.defString[i] == " " or self.defString[i] == ",":
                self.tokens.append(self.defString[i])
                i+=1
            else:
                self.tokens.append(self.defString[i:i+6])
                i+=6
        return

    def blockify(self):
        #goes inside pathway
        block = []
        b = 0
        for token in self.tokens:
            if token == "(":
                if b > 0:
                    block.append(token)
                b += 1
            elif token == ")":
                b -= 1
                if b > 0:
                    block.append(token)
            elif token == " " and b == 0:
                self.blox.append(Block(block))
                block = []
            else:
                block.append(token)

        self.blox.append(Block(block)) #last block.

    # CONSOLIDATION FX
    def unifyBlocks(self):
        # Bring together Block properties into Pathway properties
        if len(self.blox) == 1: # pathway just one gene long

            self.sourceNodes = self.blox[0].srcNodes
            self.targetNodes = self.blox[0].tgtNodes
            self.edges.add((self.blox[0].srcNodes[0],self.blox[0].tgtNodes[0]))
            self.genes.add(self.blox[0].srcNodes[0])
            


        for i in range(len(self.blox)):
            if i == 0:
                self.sourceNodes = self.blox[i].srcNodes
                target = self.blox[i].tgtNodes
                if len(self.blox)==1:
                    self.targetNodes = target
                
            else: # after first blox
                for tgtnode in target:
                    for srcnode in self.blox[i].srcNodes:
                        self.edges.add((tgtnode, srcnode))
                target = self.blox[i].tgtNodes

 
                if i == len(self.blox)-1:
                    self.targetNodes = target
            
            #add genes
            for gene in self.blox[i].genes:
                self.genes.add(gene)
            #add edges
            for edge in self.blox[i].edges:
                self.edges.add(edge)
        return

class Gph:
    def __init__(self, inputStr):
        self.pws = []
        self.defStr = self.parse(inputStr)
        self.edges = set()
        self.geneCount = 0
        self.gene2node = {}
        self.sourceNodes = set()
        self.targetNodes = set()
        self.singletons = []
        self.listOfPathways = []
        self.consolidatePathways()
        self.genGraph()
    
    # PARSING INFO
    def parse(self, inputStr):
        pathways = inputStr.split("\n")
        # for loop to load pathways into classes
        for pathway in pathways:
            self.pws.append(Pathway(pathway))
        return

    def consolidatePathways(self):
        for i in range(len(self.pws)):
            # Within a pathway now

            # Adding genes to gene dict
            for gene in self.pws[i].genes:
                if gene not in self.gene2node:
                    self.gene2node[gene] = self.geneCount
                    self.geneCount += 1
            
            for edge in self.pws[i].edges:
                self.edges.add((self.gene2node[edge[0]], self.gene2node[edge[1]]))
                if self.gene2node[edge[0]] == self.gene2node[edge[1]]:
                    self.singletons.append((self.gene2node[edge[0]], self.gene2node[edge[1]]))
                # self.edges.add(edge)
            
            for src in self.pws[i].sourceNodes:
                self.sourceNodes.add(src)
            for tgt in self.pws[i].targetNodes:
                self.targetNodes.add(tgt)
            
            #in the case of a solitary node
            if self.pws[i].targetNodes == [] and len(self.pws[i].sourceNodes) == 1:
                self.targetNodes.add(self.pws[i].sourceNodes[0])

        self.edges = list(self.edges)
        self.sourceNodes = list(self.sourceNodes)
        self.targetNodes = list(self.targetNodes)
        self.node2gene = dict((v,k) for k,v in self.gene2node.items())

    
    def genGraph(self):
        g = Graph(self.edges, directed=True)
        for i in list(self.node2gene.keys()):
            try:
                g.vs[i]["label"] = self.node2gene[i]
            except: # accounting for graphs with vertices that don't have edges
                g.add_vertex(i)
                g.vs[i]["label"] = self.node2gene[i]
            if g.degree(i) == 0:
                g.add_edge(i,i)
                self.singletons.append((i,i))
            
        # plot(g) #comment in if you want this graph plotted
        

        for s in self.sourceNodes:
            for t in self.targetNodes:
                a = g.get_all_simple_paths(self.gene2node[s],self.gene2node[t])
                
                for pathway in a:
                    self.listOfPathways.append(pathway)
  

        for edge in self.singletons:
            self.listOfPathways.append(edge)

if __name__=="__main__":
    Gph("(K01497,K14652) (K01498 K00082,K11752) (K22912,K20860,K20861,K20862,K21063,K21064)\n(K02858,K14652)\nK00794 K00793 (K20884 K22949,K11753)\nK01497 K14654 K14655\nK02858\nK00794 K00793 K00861 K00953")    #now just work on the rust function
    pass