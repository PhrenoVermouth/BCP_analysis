![graphviz](https://github.com/user-attachments/assets/b18166e1-f707-41a7-b86c-576b95cca189)# BCP Analysis Pipeline

A Nextflow-based pipeline for single-cell RNA sequencing data preprocessing, specifically designed for BCP (Billion Cell Program) analysis. The pipeline performs comprehensive single-cell data processing including alignment, ambient RNA removal, doublet detection, and basic quality control metrics.
![Upload<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="417pt" height="787pt" viewBox="0.00 0.00 417.16 787.40">
<g id="graph0" class="graph" transform="scale(1 1) rotate(0) translate(4 783.4)">
<title>dag</title>
<polygon fill="white" stroke="none" points="-4,4 -4,-783.4 413.16,-783.4 413.16,4 -4,4"/>
<!-- p0 -->
<g id="node1" class="node">
<title>p0</title>
<ellipse fill="black" stroke="black" cx="213.95" cy="-760.8" rx="1.8" ry="1.8"/>
<text text-anchor="middle" x="161.02" y="-766.8" font-family="Times,serif" font-size="14.00">Channel.fromPath</text>
</g>
<!-- p1 -->
<g id="node2" class="node">
<title>p1</title>
<ellipse fill="none" stroke="black" cx="213.95" cy="-718.4" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="187.4" y="-726.2" font-family="Times,serif" font-size="14.00">splitCsv</text>
</g>
<!-- p0&#45;&gt;p1 -->
<g id="edge1" class="edge">
<title>p0-&gt;p1</title>
<path fill="none" stroke="black" d="M213.95,-758.8C213.95,-755.3 213.95,-743.26 213.95,-733.4"/>
<polygon fill="black" stroke="black" points="217.45,-733.55 213.95,-723.55 210.45,-733.55 217.45,-733.55"/>
</g>
<!-- p2 -->
<g id="node3" class="node">
<title>p2</title>
<ellipse fill="none" stroke="black" cx="213.95" cy="-674.2" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="198.3" y="-682" font-family="Times,serif" font-size="14.00">map</text>
</g>
<!-- p1&#45;&gt;p2 -->
<g id="edge2" class="edge">
<title>p1-&gt;p2</title>
<path fill="none" stroke="black" d="M213.95,-714.37C213.95,-709.21 213.95,-698.41 213.95,-689.46"/>
<polygon fill="black" stroke="black" points="217.45,-689.5 213.95,-679.5 210.45,-689.5 217.45,-689.5"/>
</g>
<!-- p3 -->
<g id="node4" class="node">
<title>p3</title>
<ellipse fill="none" stroke="black" cx="213.95" cy="-599.8" rx="65.77" ry="18"/>
<text text-anchor="middle" x="213.95" y="-595.6" font-family="Times,serif" font-size="14.00">STAR_SOLO</text>
</g>
<!-- p2&#45;&gt;p3 -->
<g id="edge3" class="edge">
<title>p2-&gt;p3</title>
<path fill="none" stroke="black" d="M213.95,-670.45C213.95,-663.63 213.95,-645.34 213.95,-629.21"/>
<polygon fill="black" stroke="black" points="217.45,-629.55 213.95,-619.55 210.45,-629.55 217.45,-629.55"/>
<text text-anchor="middle" x="256.71" y="-640" font-family="Times,serif" font-size="14.00">ch_input_reads</text>
</g>
<!-- p5 -->
<g id="node5" class="node">
<title>p5</title>
<ellipse fill="none" stroke="black" cx="213.95" cy="-526.8" rx="105.46" ry="18"/>
<text text-anchor="middle" x="213.95" y="-522.6" font-family="Times,serif" font-size="14.00">GZIP_SOLO_OUTPUT</text>
</g>
<!-- p3&#45;&gt;p5 -->
<g id="edge4" class="edge">
<title>p3-&gt;p5</title>
<path fill="none" stroke="black" d="M213.95,-581.61C213.95,-574.03 213.95,-564.9 213.95,-556.34"/>
<polygon fill="black" stroke="black" points="217.45,-556.34 213.95,-546.34 210.45,-556.34 217.45,-556.34"/>
</g>
<!-- p13 -->
<g id="node6" class="node">
<title>p13</title>
<ellipse fill="none" stroke="black" cx="86.95" cy="-526.8" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="72.46" y="-534.6" font-family="Times,serif" font-size="14.00">mix</text>
</g>
<!-- p3&#45;&gt;p13 -->
<g id="edge5" class="edge">
<title>p3-&gt;p13</title>
<path fill="none" stroke="black" d="M173.83,-585.23C150.69,-576.14 121.91,-562.57 99.95,-544.8 98.2,-543.39 96.55,-541.69 95.06,-539.93"/>
<polygon fill="black" stroke="black" points="97.95,-537.96 89.34,-531.78 92.22,-541.97 97.95,-537.96"/>
</g>
<!-- p4 -->
<g id="node7" class="node">
<title>p4</title>
<ellipse fill="black" stroke="black" cx="338.95" cy="-526.8" rx="1.8" ry="1.8"/>
</g>
<!-- p3&#45;&gt;p4 -->
<g id="edge6" class="edge">
<title>p3-&gt;p4</title>
<path fill="none" stroke="black" d="M255.68,-585.44C278.81,-576.56 307.03,-563.13 327.95,-544.8 329.92,-543.07 331.69,-540.9 333.19,-538.69"/>
<polygon fill="black" stroke="black" points="336.11,-540.65 337.79,-530.19 329.95,-537.32 336.11,-540.65"/>
</g>
<!-- p6 -->
<g id="node8" class="node">
<title>p6</title>
<ellipse fill="none" stroke="black" cx="174.95" cy="-453.8" rx="60.4" ry="18"/>
<text text-anchor="middle" x="174.95" y="-449.6" font-family="Times,serif" font-size="14.00">SCRUBLET</text>
</g>
<!-- p5&#45;&gt;p6 -->
<g id="edge7" class="edge">
<title>p5-&gt;p6</title>
<path fill="none" stroke="black" d="M204.51,-508.61C200.1,-500.59 194.75,-490.85 189.82,-481.87"/>
<polygon fill="black" stroke="black" points="192.91,-480.24 185.03,-473.16 186.78,-483.61 192.91,-480.24"/>
</g>
<!-- p7 -->
<g id="node14" class="node">
<title>p7</title>
<ellipse fill="none" stroke="black" cx="295.95" cy="-453.8" rx="42.7" ry="18"/>
<text text-anchor="middle" x="295.95" y="-449.6" font-family="Times,serif" font-size="14.00">SOUPX</text>
</g>
<!-- p5&#45;&gt;p7 -->
<g id="edge13" class="edge">
<title>p5-&gt;p7</title>
<path fill="none" stroke="black" d="M233.38,-508.97C244.1,-499.69 257.57,-488.03 269.26,-477.91"/>
<polygon fill="black" stroke="black" points="271.37,-480.71 276.63,-471.52 266.78,-475.42 271.37,-480.71"/>
</g>
<!-- p14 -->
<g id="node11" class="node">
<title>p14</title>
<ellipse fill="none" stroke="black" cx="143.95" cy="-395.2" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="129.46" y="-403" font-family="Times,serif" font-size="14.00">mix</text>
</g>
<!-- p13&#45;&gt;p14 -->
<g id="edge20" class="edge">
<title>p13-&gt;p14</title>
<path fill="none" stroke="black" d="M87,-522.8C87.3,-511 89.53,-466.86 105.95,-435.8 112.46,-423.48 123.94,-412.31 132.58,-405"/>
<polygon fill="black" stroke="black" points="134.56,-407.89 140.23,-398.95 130.22,-402.4 134.56,-407.89"/>
</g>
<!-- p9 -->
<g id="node9" class="node">
<title>p9</title>
<ellipse fill="none" stroke="black" cx="266.95" cy="-395.2" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="252.46" y="-403" font-family="Times,serif" font-size="14.00">join</text>
</g>
<!-- p6&#45;&gt;p9 -->
<g id="edge8" class="edge">
<title>p6-&gt;p9</title>
<path fill="none" stroke="black" d="M200.58,-437.03C218.18,-426.21 240.76,-412.31 254.53,-403.84"/>
<polygon fill="black" stroke="black" points="256.04,-407.02 262.72,-398.8 252.37,-401.06 256.04,-407.02"/>
</g>
<!-- p17 -->
<g id="node10" class="node">
<title>p17</title>
<ellipse fill="none" stroke="black" cx="204.95" cy="-370" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="189.3" y="-377.8" font-family="Times,serif" font-size="14.00">map</text>
</g>
<!-- p6&#45;&gt;p17 -->
<g id="edge9" class="edge">
<title>p6-&gt;p17</title>
<path fill="none" stroke="black" d="M182.75,-435.48C185.17,-429.89 187.77,-423.62 189.95,-417.8 194.04,-406.85 198.05,-394.21 200.89,-384.83"/>
<polygon fill="black" stroke="black" points="204.19,-386.02 203.67,-375.44 197.48,-384.04 204.19,-386.02"/>
</g>
<!-- p6&#45;&gt;p14 -->
<g id="edge10" class="edge">
<title>p6-&gt;p14</title>
<path fill="none" stroke="black" d="M165.65,-435.82C160.87,-427.1 155.16,-416.67 150.79,-408.69"/>
<polygon fill="black" stroke="black" points="153.97,-407.21 146.1,-400.12 147.83,-410.57 153.97,-407.21"/>
</g>
<!-- p15 -->
<g id="node12" class="node">
<title>p15</title>
<ellipse fill="none" stroke="black" cx="152.95" cy="-330.4" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="138.46" y="-338.2" font-family="Times,serif" font-size="14.00">mix</text>
</g>
<!-- p6&#45;&gt;p15 -->
<g id="edge11" class="edge">
<title>p6-&gt;p15</title>
<path fill="none" stroke="black" d="M172.87,-435.5C170.65,-417.99 166.78,-390.18 161.95,-366.4 160.51,-359.34 158.55,-351.57 156.83,-345.15"/>
<polygon fill="black" stroke="black" points="160.26,-344.43 154.2,-335.74 153.52,-346.31 160.26,-344.43"/>
</g>
<!-- p16 -->
<g id="node13" class="node">
<title>p16</title>
<ellipse fill="none" stroke="black" cx="152.95" cy="-271.8" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="138.46" y="-279.6" font-family="Times,serif" font-size="14.00">mix</text>
</g>
<!-- p6&#45;&gt;p16 -->
<g id="edge12" class="edge">
<title>p6-&gt;p16</title>
<path fill="none" stroke="black" d="M156.38,-436.27C146.8,-426.46 136.06,-413.16 130.95,-398.8 116.76,-358.93 135.89,-308.58 146.6,-285.46"/>
<polygon fill="black" stroke="black" points="149.66,-287.16 150.94,-276.64 143.38,-284.06 149.66,-287.16"/>
</g>
<!-- p10 -->
<g id="node16" class="node">
<title>p10</title>
<ellipse fill="none" stroke="black" cx="266.95" cy="-330.4" rx="50.21" ry="18"/>
<text text-anchor="middle" x="266.95" y="-326.2" font-family="Times,serif" font-size="14.00">SAM_QC</text>
</g>
<!-- p9&#45;&gt;p10 -->
<g id="edge16" class="edge">
<title>p9-&gt;p10</title>
<path fill="none" stroke="black" d="M266.95,-391.17C266.95,-385.45 266.95,-372.53 266.95,-360.15"/>
<polygon fill="black" stroke="black" points="270.45,-360.22 266.95,-350.22 263.45,-360.22 270.45,-360.22"/>
</g>
<!-- p18 -->
<g id="node20" class="node">
<title>p18</title>
<ellipse fill="none" stroke="black" cx="188.95" cy="-221.4" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="174.46" y="-229.2" font-family="Times,serif" font-size="14.00">mix</text>
</g>
<!-- p17&#45;&gt;p18 -->
<g id="edge24" class="edge">
<title>p17-&gt;p18</title>
<path fill="none" stroke="black" d="M204.61,-365.94C202.8,-349.32 194.09,-269.49 190.5,-236.6"/>
<polygon fill="black" stroke="black" points="194.01,-236.51 189.45,-226.95 187.05,-237.27 194.01,-236.51"/>
</g>
<!-- p14&#45;&gt;p15 -->
<g id="edge21" class="edge">
<title>p14-&gt;p15</title>
<path fill="none" stroke="black" d="M144.38,-391.17C145.54,-383.08 148.77,-360.54 150.92,-345.57"/>
<polygon fill="black" stroke="black" points="154.37,-346.15 152.32,-335.75 147.44,-345.15 154.37,-346.15"/>
</g>
<!-- p15&#45;&gt;p16 -->
<g id="edge22" class="edge">
<title>p15-&gt;p16</title>
<path fill="none" stroke="black" d="M152.95,-326.36C152.95,-319.14 152.95,-300.6 152.95,-287.36"/>
<polygon fill="black" stroke="black" points="156.45,-287.37 152.95,-277.37 149.45,-287.37 156.45,-287.37"/>
</g>
<!-- p16&#45;&gt;p18 -->
<g id="edge23" class="edge">
<title>p16-&gt;p18</title>
<path fill="none" stroke="black" d="M154.88,-268.2C159.62,-261.84 172,-245.19 180.46,-233.82"/>
<polygon fill="black" stroke="black" points="183.11,-236.11 186.27,-226 177.5,-231.93 183.11,-236.11"/>
</g>
<!-- p7&#45;&gt;p9 -->
<g id="edge14" class="edge">
<title>p7-&gt;p9</title>
<path fill="none" stroke="black" d="M287.25,-435.82C282.83,-427.2 277.56,-416.91 273.5,-408.98"/>
<polygon fill="black" stroke="black" points="276.64,-407.45 268.97,-400.14 270.41,-410.64 276.64,-407.45"/>
</g>
<!-- p8 -->
<g id="node15" class="node">
<title>p8</title>
<ellipse fill="black" stroke="black" cx="295.95" cy="-395.2" rx="1.8" ry="1.8"/>
</g>
<!-- p7&#45;&gt;p8 -->
<g id="edge15" class="edge">
<title>p7-&gt;p8</title>
<path fill="none" stroke="black" d="M295.95,-435.51C295.95,-426.91 295.95,-416.69 295.95,-408.83"/>
<polygon fill="black" stroke="black" points="299.45,-408.98 295.95,-398.98 292.45,-408.98 299.45,-408.98"/>
</g>
<!-- p19 -->
<g id="node17" class="node">
<title>p19</title>
<ellipse fill="none" stroke="black" cx="232.95" cy="-246.6" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="217.3" y="-254.4" font-family="Times,serif" font-size="14.00">map</text>
</g>
<!-- p10&#45;&gt;p19 -->
<g id="edge17" class="edge">
<title>p10-&gt;p19</title>
<path fill="none" stroke="black" d="M259.9,-312.45C253.54,-297.13 244.29,-274.89 238.45,-260.84"/>
<polygon fill="black" stroke="black" points="241.75,-259.66 234.68,-251.77 235.29,-262.35 241.75,-259.66"/>
</g>
<!-- p11 -->
<g id="node18" class="node">
<title>p11</title>
<ellipse fill="black" stroke="black" cx="266.95" cy="-271.8" rx="1.8" ry="1.8"/>
</g>
<!-- p10&#45;&gt;p11 -->
<g id="edge18" class="edge">
<title>p10-&gt;p11</title>
<path fill="none" stroke="black" d="M266.95,-312.11C266.95,-303.51 266.95,-293.29 266.95,-285.43"/>
<polygon fill="black" stroke="black" points="270.45,-285.58 266.95,-275.58 263.45,-285.58 270.45,-285.58"/>
</g>
<!-- p20 -->
<g id="node21" class="node">
<title>p20</title>
<ellipse fill="none" stroke="black" cx="205.95" cy="-177.2" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="191.46" y="-185" font-family="Times,serif" font-size="14.00">mix</text>
</g>
<!-- p19&#45;&gt;p20 -->
<g id="edge26" class="edge">
<title>p19-&gt;p20</title>
<path fill="none" stroke="black" d="M231.79,-242.7C228.28,-233.96 217.65,-207.41 211.17,-191.23"/>
<polygon fill="black" stroke="black" points="214.6,-190.4 207.64,-182.42 208.11,-193 214.6,-190.4"/>
</g>
<!-- p12 -->
<g id="node19" class="node">
<title>p12</title>
<ellipse fill="black" stroke="black" cx="86.95" cy="-599.8" rx="1.8" ry="1.8"/>
<text text-anchor="middle" x="42.57" y="-605.8" font-family="Times,serif" font-size="14.00">Channel.empty</text>
</g>
<!-- p12&#45;&gt;p13 -->
<g id="edge19" class="edge">
<title>p12-&gt;p13</title>
<path fill="none" stroke="black" d="M86.95,-597.57C86.95,-591.08 86.95,-560.53 86.95,-542.07"/>
<polygon fill="black" stroke="black" points="90.45,-542.25 86.95,-532.25 83.45,-542.25 90.45,-542.25"/>
</g>
<!-- p18&#45;&gt;p20 -->
<g id="edge25" class="edge">
<title>p18-&gt;p20</title>
<path fill="none" stroke="black" d="M190.06,-217.63C192.18,-212.39 196.89,-200.7 200.63,-191.4"/>
<polygon fill="black" stroke="black" points="203.75,-193.02 204.24,-182.44 197.26,-190.41 203.75,-193.02"/>
</g>
<!-- p21 -->
<g id="node22" class="node">
<title>p21</title>
<ellipse fill="none" stroke="black" cx="205.95" cy="-133" rx="3.6" ry="3.6"/>
<text text-anchor="middle" x="183.69" y="-140.8" font-family="Times,serif" font-size="14.00">collect</text>
</g>
<!-- p20&#45;&gt;p21 -->
<g id="edge27" class="edge">
<title>p20-&gt;p21</title>
<path fill="none" stroke="black" d="M205.95,-173.17C205.95,-168.01 205.95,-157.21 205.95,-148.26"/>
<polygon fill="black" stroke="black" points="209.45,-148.3 205.95,-138.3 202.45,-148.3 209.45,-148.3"/>
</g>
<!-- p23 -->
<g id="node23" class="node">
<title>p23</title>
<ellipse fill="none" stroke="black" cx="259.95" cy="-58.6" rx="55.03" ry="18"/>
<text text-anchor="middle" x="259.95" y="-54.4" font-family="Times,serif" font-size="14.00">MULTIQC</text>
</g>
<!-- p21&#45;&gt;p23 -->
<g id="edge28" class="edge">
<title>p21-&gt;p23</title>
<path fill="none" stroke="black" d="M205.13,-129.05C203.47,-122.4 200.36,-105.97 206.64,-94.6 209.65,-89.15 213.91,-84.4 218.72,-80.31"/>
<polygon fill="black" stroke="black" points="220.61,-83.27 226.62,-74.55 216.48,-77.62 220.61,-83.27"/>
<text text-anchor="middle" x="249.79" y="-98.8" font-family="Times,serif" font-size="14.00">ch_for_multiqc</text>
</g>
<!-- p25 -->
<g id="node25" class="node">
<title>p25</title>
<ellipse fill="black" stroke="black" cx="247.95" cy="-1.8" rx="1.8" ry="1.8"/>
</g>
<!-- p23&#45;&gt;p25 -->
<g id="edge30" class="edge">
<title>p23-&gt;p25</title>
<path fill="none" stroke="black" d="M256.15,-40.26C254.36,-32.07 252.26,-22.48 250.63,-15.06"/>
<polygon fill="black" stroke="black" points="254.09,-14.49 248.53,-5.47 247.25,-15.99 254.09,-14.49"/>
</g>
<!-- p24 -->
<g id="node26" class="node">
<title>p24</title>
<ellipse fill="black" stroke="black" cx="272.95" cy="-1.8" rx="1.8" ry="1.8"/>
</g>
<!-- p23&#45;&gt;p24 -->
<g id="edge31" class="edge">
<title>p23-&gt;p24</title>
<path fill="none" stroke="black" d="M264.06,-40.26C266.01,-32.07 268.28,-22.48 270.04,-15.06"/>
<polygon fill="black" stroke="black" points="273.41,-16.01 272.32,-5.47 266.6,-14.39 273.41,-16.01"/>
</g>
<!-- p22 -->
<g id="node24" class="node">
<title>p22</title>
<ellipse fill="black" stroke="black" cx="314.95" cy="-133" rx="1.8" ry="1.8"/>
<text text-anchor="middle" x="262.02" y="-139" font-family="Times,serif" font-size="14.00">Channel.fromPath</text>
</g>
<!-- p22&#45;&gt;p23 -->
<g id="edge29" class="edge">
<title>p22-&gt;p23</title>
<path fill="none" stroke="black" d="M314.48,-131.03C312.32,-126.61 303.08,-108.1 292.95,-94.6 290.39,-91.19 287.54,-87.76 284.62,-84.43"/>
<polygon fill="black" stroke="black" points="287.24,-82.11 277.89,-77.12 282.09,-86.85 287.24,-82.11"/>
<text text-anchor="middle" x="356.28" y="-98.8" font-family="Times,serif" font-size="14.00">ch_multiqc_config</text>
</g>
</g>
</svg>ing graphviz.svg…]()



## Features

- **Single-cell alignment** using STAR
- **Ambient RNA removal** for cleaner expression profiles
- **Doublet detection** to filter multiplets
- **Comprehensive QC metrics** with MultiQC reporting
- **Embedded QC plots** included in MultiQC report
- **Automated preprocessing** with minimal manual intervention

## Future Development

Parameters will be optimized and made configurable for different species and organ types to enhance pipeline flexibility and accuracy.

## Environment Configuration

### Current Setup
The pipeline currently requires manual conda environment configuration for testing and development purposes. 

### Planned Updates
Singularity container support will be implemented upon completion of the development phase to ensure better reproducibility and easier deployment across different computing environments.

## Installation

```bash
# Clone the repository
git clone https://github.com/PhrenoVermouth/BCP_analysis.git
cd BCP_analysis

# Create conda environment
mamba env create -f bin/environment.yml
```

## Usage

```bash
# Run the pipeline
nextflow run ~/BCP_analysis/main.nf -profile standard
```

## Output

The pipeline generates:
- Processed single-cell count matrices
- Quality control reports via MultiQC
- Doublet detection results
- Ambient RNA removal metrics
