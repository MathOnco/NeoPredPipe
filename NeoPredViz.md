## Web Based Visualizations using NeoPredViz

#### Additional Dependencies For Visualizations

1. [bokeh](http://bokeh.pydata.org/en/latest/docs/installation.html)

2. [holoviews](http://holoviews.org) (pip install holoviews)

### Executing NeoPredViz

You can use NeoPredViz to visualize some or all of the information in the output files consistent with the modularity of execution.

```bash
python NeoPredViz.py --help
```

```bash
usage: NeoPredViz.py [-h] [-s SUMMARYFILE] [-o OUTPUT] [-n NEOSFILE]
                     [-r NEORECOFILE]

optional arguments:
  -h, --help      show this help message and exit

Required arguments:
  -s SUMMARYFILE  Neoantigens Summary File
  -o OUTPUT       Output Directory Path

Optional Arguments:
  -n NEOSFILE     Neoantigen Predictions Output File. Default: None.
  -r NEORECOFILE  Neoantigen Recognition Potentials Output File. Default:
                  None.

```

The NeoPredViz output will produce a simple web based visualization for you to explore outputs quickly.
The use all of the information that NeoPredPipe can output simply run:

```bash
python NeoPredViz.py -o NeoantigenVizDashboard.html -s SummaryTable.txt -n NeoantigensOutput.txt -r NeoRecoPoTable.txt
```

Finally, in order to visualize, open the output file (NeoantigenVizDashboard.html) with your preferred web browser. This file can also be shared as an html with others.

It's possible to export these plots, but for this, please refer to [bokeh documentation](http://bokeh.pydata.org/). 
