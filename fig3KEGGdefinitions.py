"""
KEGG definitions import graphlib
from https://www.genome.jp/module/

If there are multiple pathways for a single metabolism, I concatenated
their definitions with a "\n" delimiter. These definitions get parsed into a graph
in defintion_parser.py

Author: Khashiff Miranda (kkmiranda.github.io)
"""

from definition_parser import *

# SULFUR METABOLISMS (https://www.genome.jp/pathway/map00920)
Asr = "(K13811,K00958+K00860,K00955+K00957,K00956+K00957+K00860) K00390 (K00380+K00381,K00392)" # M00176
Dsr = "K00958 (K00394+K00395) (K11180+K11181)" # M00596
TsO = "K17222+K17223+K17224-K17225-K22622+K17226+K17227" # M00595

# NITROGEN METABOLISMS https://www.genome.jp/pathway/map00910
Nif = "K02588+K02586+K02591-K00531,K22896+K22897+K22898+K22899" # M00175 
Anr = "(K00367,K10534,K00372-K00360) (K00366,K17877)" # https://www.genome.jp/module/M00531
Dnr = "(K00370+K00371+K00374,K02567+K02568) (K00362+K00363,K03385+K15876)" # Https://www.genome.jp/module/M00530
Denitrification = "(K00370+K00371+K00374,K02567+K02568) (K00368,K15864) (K04561+K02305) K00376" # https://www.genome.jp/module/M00529
Nitrification = "K10944+K10945+K10946 K10535" # https://www.genome.jp/module/M00528
coammox = "K10944+K10945+K10946 K10535 K00370+K00371" # https://www.genome.jp/module/M00804
anammox = "(K00368,K15864) (K20932,K20933,K20934) K20935" # derived from https://www.genome.jp/pathway/map00910

# B-VITAMIN PRODUCTION

## THIAMIN (Vit. B1) 
## M00898, M00897, M00896, M00895, M00127
thiamin = "K03146\nK18278 K00877 K14154\nK03147 K14153 K22911 K00949\nK22699\nK03147 (K00941 (K00788,K21220),K21219) K00946\nK03148+K03154 K03151\nK03153 K03149 K10810\nK03147 K00941 K00788 K00946\nK03148+K03154 K03151\nK03150 K03149\nK03147 ((K00941 K00788),K14153,K21219) K00946"

## RIBOFLAVIN (Vit. B2): 
## M00911,M00125
ribo = "(K01497,K14652) (K01498 K00082,K11752) (K22912,K20860,K20861,K20862,K21063,K21064)\n(K02858,K14652)\nK00794 K00793 (K20884 K22949,K11753)\nK01497 K14654 K14655\nK02858\nK00794 K00793 K00861 K00953"

## BIOTIN (Vit. B7)
## M00123
biotin = "K00652 ((K00833,K19563) K01935,K19562) K01012"

## COBALAMIN (Vit. B12)
## Anaerobic:   M00924
## Aerobic:     M00925
cob_aerobic = "(K02303,K13542) (K03394,K13540) K02229 (K05934,K13540,K13541) K05936 K02228 K05895 K00595 K06042 K02224 K02230+K09882+K09883"
cob_anaerobic = "(K02302,(K02303,K13542) (K02304,K24866)) (K02190,K03795,K22011) K03394 (K05934,K13541,K21479) K05936 (K02189,K13541) K02188 K05895 (K02191 K03399,K00595) K06042 K02224"

defaultMultiDef = {
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