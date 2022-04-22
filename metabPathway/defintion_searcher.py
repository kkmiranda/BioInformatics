"""
Developing this as a one time use to download all KEGG definitions

Downloads the JSON file from the genome.jp website. Meant to be used only once.
"""

from importlib.resources import path
import requests as req
import re
import json
import os


BASEPATH = "/Users/khashiffm/Documents/Research/Cathylab/metagenomes/nifH/masterdata_17May2021/heatmapGen/pathwayAlgorithm/"


with open(f'{BASEPATH}/KEGGPathways.json','r') as t:
    data = json.load(t)

# print(data["children"][0]["name"])

print("!!! Finding Pathway IDs")

# 0th item in list is metabolism pathways
pathways = dict()
for i in data["children"][0]["children"]:
    category = re.findall(r".+([A-Z]\w.+)", i["name"])[0]
    pathways[category] = dict()
    print(category)
    for j in i["children"]:
        try:
            match = re.match(r'(\d{5}).(.+)\s\[PATH',j["name"]).groups()
            pathways[category][match[0]] = dict()
            pathways[category][match[0]]["name"] = match[1]
            pathways[category][match[0]]["modules"] = dict()
        except:
            pass #harmless

print(pathways.keys())

def findDef(key: str) -> tuple:
    """
    Using KEGG Module number, scrapes out name and definition
    """

    resp = req.get(f"https://www.genome.jp/module/{key}")

    defName = re.findall(r"<td>Name<\/td>\s+?<td>(.+)<\/td>",resp.text)[0]
    KEGGdef = re.findall(r"<td>Definition<\/td>\s+?<td>\s+?(.+)<br", resp.text)

    #strip HTML tags
    strippedDef = re.sub('<[^<]+?>', '', KEGGdef[0].lstrip())

    return((defName,strippedDef))

print("!!! Finding Module IDs and Module Info")
# go to pathway to get all modules
failed = []
for metab in pathways.keys():
    for pw in pathways[metab].keys():
        resp = req.get(f"https://www.genome.jp/entry/pathway+ko{pw}")
        for mod in set(re.findall(r"(M\d{5})",resp.text)):
            pathways[metab][pw]["modules"][mod] = dict()
            try:
                modInfo = findDef(mod)
                pathways[metab][pw]["modules"][mod]["name"] = modInfo[0]
                pathways[metab][pw]["modules"][mod]["definition"] = modInfo[1]
            except:
                failed.append((metab,pw,mod))
                print("FAILED with ",metab,pw, mod)
print("weeeee\n\n\n")

json_obj = json.dumps(pathways, indent=4)
with open(f"{BASEPATH}/all_KEGG_defs.json",'w') as o:
    o.write(json_obj)

for i in failed:
    print(i)
