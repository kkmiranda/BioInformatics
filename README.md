# Network Algorithm for Metabolism Detection

This program uses KEGG definitions from the genome.jp database to identify the presence/absence of metabolisms across genomes. Prerequisites to using this program include using Prodigal to identify open-reading frames (ORFs) in metagenomic data and annotating these open reading frames using the KEGG database. To make life easier, we simply used the anvio environment to do the above as it also gave us a really functional interface that we used to bin the DNA into metagenome-assembled genomes (MAGs) (optional for this code). Once this was done, we extracted data from the anvio PROFILE and CONTIGS databases in which we stored our bins. For an example of the whole workflow, check out this <a href="https://github.com/kkmiranda/PNWMetagenomes/">project</a> where we applied this algorithm! 

## An overview

This metabolism pathway detection goes beyond the simple absence/presence of annotation softwares. Each metabolism is comprised of many genes that code for a protein which is the worker that accomplishes a metabolic task. If every metabolism was a simple linear assembly line, this code would be fairly useless, too complex and time-consuming as a simple excel-tally macro would get the job done. However, evolution has made the job so much more exciting for computer scientists like myself! Evolution has allowed metabolisms like Vitamin B12 synthesis diverge to having alternate pathways (<a href="https://www.genome.jp/module/M00925">aerobic</a> vs <a href = "https://www.genome.jp/module/M00924">anaerobic</a>). 

These alternate pathways allow for the same metabolism to be conducted in different environmental conditions, in this case across a range of oxygen environments. Similarly, evolution has allowed for degeneracy across metabolic intermediates where different genes could encode for the same protein OR different protein-intermediates could achieve the same metabolism. This makes the question of whether a metabolism exists in a genome a lot harder. 

## Graphical approach to metabolisms
What I employ here is a network based simple paths approach. Each metabolism has a KEGG definition. I parse this definition to generate a network-graph connecting all the genes in a metabolic pathway in <a href="./definition_parser.py">`definition_parser.py`</a>. Using g.get_all_simple_paths(), I can then get every possible path from the first genes to the final genes of the pathway. If my annotated Metagenome has all the genes necessary to traverse the entire graph (metabolism), the metabolism exists!

### A denitrifying example
Denitrification is a pathway (that my mum did her PhD. on in the 90's!) that microbes use to respire ammonium when oxygen isn't available. It converts NH4+ to N2 gas and is of interest in sewage treatment and remediating high nitrogen loading in groundwater. On a bioinformatics standpiont, its KEGG module number is M00529.

<img src="./assets/denitrification.png"/>

Two representations of the dentrification pathway: The boxy figure on the left is found on the genome.jp website while the network on the right is what NAMeD generates from the KEGG Defintion:

(K00370+K00371+K00374,K02567+K02568) (K00368,K15864) (K04561+K02305) K00376*

I set it up as the red-noded graph in igraph and from there I generate all possible paths from the source nodes (K00370, K02567) to the target nodes (K00376) using `get_all_simple_paths()` (L252 from <a href="./definition_parser.py">`definition_parser.py`</a>). In this specific metabolism, there are four possible pathways of genes which if present indicates that a genome may carry out denitrification:
<ul>
<li>['K02567', 'K02568', 'K00368', 'K04561', 'K02305', 'K00376']</li>

<li>['K02567', 'K02568', 'K15864', 'K04561', 'K02305', 'K00376']</li>

<li>['K00370', 'K00371', 'K00374', 'K00368', 'K04561', 'K02305', 'K00376']</li>

<li>['K00370', 'K00371', 'K00374', 'K15864', 'K04561', 'K02305', 'K00376']</li>
</ul>

## Missing threshholds
From here, we just need to see if all the genes in any of these pathways are present. Simple! However, a lot of our research is conducted in novel environments where the products of evolution are yet to be explored. Perhaps the gene homologous to K00370 (narG, narZ, nxrA) in the denitrification pathway occuring in Pacific Northwest seagrass rhizomes has evolved to such a state that it won't be picked up by Prodigal's gene annotation software. Or maybe KEGG hasn't accounted for this novel genetic diversity in their database. In this case, the genome <i>can</i> denitrify nitrate but our software would fail as it wouldn't pick up the gene doing K00370's job. I account for this by allowing NAMeD to say a metabolism is present even if a user-defined maximum __% of genes in a pathway is missing. 

For example: I set a threshhold allowing for at most 33% of genes missing (DEFAULT). In the last pathway above, the length is 7 genes (n = 7). So this threshold would allow at most 2 gene (m = 2) to be absent (2.31 gets rounded down to 2 to stay conservative!) This floor function tends to make this program sensitive to pathways of smaller lengths (if n = 2, m = 0) and more allowing towards pathways of longer lengths (n = 15, m = 3). Keep this in mind when interpreting your data!

<ul>
<li>`threshold = 0` to look for the presence of EVERY gene</li>
<li>`threshold = 0.33` DEFAULT</li>
<li>`threshold = 1` for chaos (every metabolism will show up as present)</li>
</ul>

## Output

I'm giving you a couple of choices with your output. Either you run the output by fiddling with the `miscList` parameter. If `miscList=True`, the output is a tab delim file that indicates 'Y' or 'N' for the presence of a metabolism. If `miscList=False`, the output is a tab delim file that contains a detailed list of all every pathway detected as present, along with how many genes were present/absent (useful for understanding your output).

This is super useful in understanding why the metabolism of your dreams isn't showing up. It's a good practice to try out different thresholds to the point you're comfortable**. 

The final miscTable can be worked up as a heatmap in R or Python with the aesthetics of your choice! 


## Modes

We've got a few default modes to run against your masterdata. If you've not got your masterdata file, head on over to our figshare 10.6084/m9.figshare.20152949 and pick up `tatoosh_masterdata.txt` and `tatoosh_MAG_names.txt` and place them in the `assets/` directory.

Run it using:
```
python NAMeD.py -i assets/tatoosh_masterdata.txt
```

You can toggle with thresholds using the `-t` flag. For example `-t 0.25` allows 25% of pathways to be missing for a metabolism to be declared present. 

### <i>DEFAULT: Single Metabolism mode</i>

```
python NAMeD.py -i assets/tatoosh_masterdata.txt -m single
```

### <i>Multiple Metabolisms mode</i>

```
python NAMeD.py -i assets/tatoosh_masterdata.txt -m multi
```


### <i>EVERY METABOLISM MODE (supercharge your analysis)</i>
I've taken the liberty to give a user the possibility of testing out their masterdata against every module that exists in the KEGG database. To save you the time and effort of having to download all definitions and convert them into networks within a usable structure, I did all that and then <a href="https://www.geeksforgeeks.org/understanding-python-pickling-example/">pickled</a> the output file. This pickled file massively saves on time loading this otherwise huge database.

```
python NAMeD.py -i assets/tatoosh_masterdata.txt -m all
```

### Customize your input
If you've got a single definition of interest you want to test out for example nitrogen fixation (M00175): 

K02588+K02586+K02591-K00531,K22896+K22897+K22898+K22899

You can input this definition as a string using the `-d` flag. Don't forget the quotation marks!

```
python NAMeD.py -i assets/tatoosh_masterdata.txt -m single -d "K02588+K02586+K02591-K00531,K22896+K22897+K22898+K22899"
```

OR

If you have a set of definitions you want to test for, refer to <a href="./fig3KEGGdefinitions.py">fig3KEGGdefinitions.py</a> and input your definitions in that style in a new python doc called `testDef.py`

Then using the <i>multi</i> mode and the `-d` flag:
```
python NAMeD.py -i assets/tatoosh_masterdata.txt -m multi -d testDef.py
```


## How is this useful/different?

I developed this code to account for a couple factors:

<ul>
<li><b>OPEN SOURCE: </b>Other metabolism detection codes weren't open source. I didn't know/could control what was going on behind the scenes which becomes dangerous in bridging the gap between the data and the biology</li>
<li><b>Meandering pathways: </b>From looking across various definitions in the KEGG module database, I noticed that complex pathways with multiple branches can't be solved with a simple tally count. <i>e.g.</i>: A star-like metabolic map with 5 terminal nodes connected by an edge each to 1 central node has a total of 6 genes. But the presence of just 2 would qualify the presence of the metabolism. When this happens in long pathways with many branches, this approach allows us to grapple the complexity by testing for every possibility</li>
<li><b>Connecting separated trees:</b> Another feature I noticed in KEGG definitions is that it would have multiple distinct unconnected trees even if genes were repeated across different trees. This algorithm connects all these trees into a single graph (as far as possible), allowing for even more possibilities. <i>eg:</i> Riboflavin synthesis has two Modules (M00911, M00125) that together comprise of 6 trees. I used this software to combine both modules into one (refer to `KEGG_definitions.py`) and it compiles these 6 trees into 2 (with one singleton). This allows for a larger range of possibilities to detect the presence of a metabolism! <img src="../assets/riboflavin.png"></li>
<li><b>Same metabolism, different path lengths:</b> From the above graph, you'd note that a genome may complete the same synthesis of Riboflavin using pathways of length 1, 3, 4 or 5. By generating each possible pathway in the now connected trees, I can apply my static thresholds (25%) across a range of pathway lengths. This allows for a more dynamic ability to detect a pathway. This is super useful when long pathways have shortcuts going through them.</li>
</ul>

<i>* The rules to interpret these definitions get hairy with '+' and '-' and positions of commas, spaces and parentheses. Refer to the <a href="https://www.genome.jp/kegg/module.html">KEGG Module documentation</a></i>

<i>** Most thresholds in bioinformatics are arbitrary! Let reason dictate your way forward</i>