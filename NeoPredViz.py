import sys
from itertools import combinations
import argparse

import pandas as pd
import numpy as np

# bokeh
from bokeh.io import show, output_file
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure
from bokeh.core.properties import value
from bokeh.models.widgets import DataTable, TableColumn, Panel, Div, Select
from bokeh.layouts import widgetbox
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar

# holoviews
import holoviews as hv
from holoviews import opts, dim

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument("-s", dest="SummaryFile", default=None, type=str, help="Neoantigens Summary File")
    requiredNamed.add_argument("-o", dest="Output", default="./NeoPredViz_Output.html", type=str, help="Output Directory Path")
    postProcess = parser.add_argument_group('Optional Arguments')
    postProcess.add_argument("-n", dest="NeosFile", default=None, type=str, help="Neoantigen Predictions Output File. Default: None.")
    postProcess.add_argument("-r", dest="NeoRecoFile", default=None, type=str, help="Neoantigen Recognition Potentials Output File. Default: None.")
    Options = parser.parse_args()  # main user args

    if not Options.SummaryFile:
        parser.error("Some of the required arguments were not provided. Please check required arguments.")

    if Options.Output == "./NeoPredViz_Output.html":
        print("INFO: Vizualization webpage placed in execution directory as NeoPredViz_Output.html")
    if not Options.NeosFile or not Options.NeoRecoFile:
        print("INFO: Only visualization partial information as at least one input file was omitted by user.")

    return(Options)

class Data:

    def __init__(self, summaryNeosFile=None, allNeosFile=None, recopo=None):
        self.summaryData = pd.read_csv(summaryNeosFile, sep="\t")
        if allNeosFile is not None:
            self.neosData = pd.read_csv(allNeosFile, sep="\t")
        if recopo is not None:
            self.recopoData = pd.read_csv(recopo, sep="\t")
        else:
            self.recopoData = None
        self.sharedNeos = None
        self.matchedNeos = None
        if allNeosFile is not None:
            with open(allNeosFile, 'r') as inFile:
                self.allLines = [line.replace('\n','') for line in inFile.readlines()]
            self.deconstructedRegions, self.sharedCount = self.deconstruct()

    def deconstruct(self):
        out = [ ]

        regions = {}
        for line in self.allLines:
            if line.startswith('Sample') == False:
                line = line.split('\t')
                n = len(line)

                regionCount = 0
                done=False
                for i, item in enumerate(line):
                    if item == '1' or item == '0' and done==False and item != '-1':
                        regionCount+=1
                    elif item.startswith('line') == True and done==False:
                        done = True
                        end = i
                        regions.update({line[0]:regionCount})

                for i in range(1,regionCount+1):
                    if line[i] == '1':
                        out.append([line[0]+'_R'+ str(i)] + line[end:len(line)])


        samples = [item[0] for item in out]
        epitopes = [item[10] for item in out]
        samplePairs = [",".join(map(str, comb)) for comb in combinations(list(set(samples)), 2)]
        shared = dict.fromkeys(samplePairs)
        for key in shared:
            shared[key]=0

        # For each epitope get all samples with that epitope
        epiFound = dict.fromkeys(list(set(epitopes)))
        for key in epiFound:
            epiFound[key]=[]
        for i, epi in enumerate(epitopes):
            epiFound[epi].append( samples[i] )

        # get only only those that are shared across two regions minimum
        keepEpis = {}
        for samples in epiFound:
            if len(epiFound[samples]) > 1:
                keepEpis.update({samples : epiFound[samples]})
        for epi in keepEpis:
            # Get all pairs of the samples present
            thepairs = [",".join(map(str, comb)) for comb in combinations(list(set(keepEpis[epi])), 2)]
            for pair in thepairs:
                try:
                    shared[pair] += 1
                except KeyError:
                    shared[pair.split(',')[1]+','+pair.split(',')[0]]

        return(out, shared)

    def EpitopeTable(self):
        Columns = [TableColumn(field=Ci, title=Ci) for Ci in self.neosData.columns]  # bokeh columns
        data_table = DataTable(columns=Columns, source=ColumnDataSource(self.neosData) ,width=1200, height=200)  # bokeh table

        return(data_table)

    def SummaryTable(self):

        Columns = [TableColumn(field=Ci, title=Ci) for Ci in self.summaryData.columns]  # bokeh columns
        data_table = DataTable(columns=Columns, source=ColumnDataSource(self.summaryData) ,width=1200, height=200)  # bokeh table

        return(data_table)

    def SummaryBarChart(self):
        self.summaryData.sort_values(by=['Total'], inplace=True)
        self.summaryData.reset_index(drop=True, inplace=True)

        # Get factors for each Sample
        factors = []
        shared = []
        clonal = []
        subclonal = []
        self.summaryData.rename(index=self.summaryData.Sample, inplace=True)
        for sam in self.summaryData.Sample:
            factors.append( (sam,"Total") )
            clonal.append(self.summaryData.loc[sam,'Clonal'])
            subclonal.append(self.summaryData.loc[sam,'Subclonal'])
            shared.append(self.summaryData.loc[sam,'Shared'])

            factors.append( (sam,"WB") )
            clonal.append(self.summaryData.loc[sam,'Clonal_WB'])
            subclonal.append(self.summaryData.loc[sam,'Subclonal_WB'])
            shared.append(self.summaryData.loc[sam,'Shared_WB'])

            factors.append( (sam,"SB") )
            clonal.append(self.summaryData.loc[sam,'Clonal_SB'])
            subclonal.append(self.summaryData.loc[sam,'Subclonal_SB'])
            shared.append(self.summaryData.loc[sam,'Shared_SB'])

        clonality = ['clonal','subclonal','shared']

        source = ColumnDataSource(data=dict(
            x=factors,
            clonal=clonal,
            subclonal=subclonal,
            shared=shared
        ))

        TOOLTIPS = [('Clonal', '@clonal'), ('Subclonal', '@subclonal'), ('Shared', '@shared')]

        p = figure(x_range=FactorRange(*factors), tooltips=TOOLTIPS, height=400)
        clonalityColors = ["#00a4ed","#ef4f25","#ede614"]
        p.vbar_stack(clonality, x='x',width=0.9,color=clonalityColors, source=source, legend=[value(x) for x in clonality])
        p.xaxis.major_label_orientation = np.pi / 2
        p.xaxis.group_label_orientation = np.pi / 2
        p.legend.orientation = "horizontal"
        p.legend.location = "top_center"
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        p.xaxis.major_tick_line_color = None
        p.yaxis.axis_label="Neoantigens"
        p.x_range.range_padding = 0.1
        return(p)

    def ChordDiagram(self):
        # Step 1 Get Data
        self.sharedNeos = self.GetShared()
        self.matchedNeos = self.GetMatchedNeos()

        hv.output(size=200)

        source = []
        target = []
        value = []
        for i,sam in enumerate(self.sharedNeos):
            for pair in self.sharedNeos[sam]:
                source.append( sam+"_"+pair.split(',')[0] )
                target.append( sam+"_"+pair.split(',')[1] )
                value.append(self.sharedNeos[sam][pair])

        for matched in self.matchedNeos:
            source.append(matched.split(',')[0])
            target.append(matched.split(',')[1])
            value.append(self.matchedNeos[matched])

        links = pd.DataFrame({ 'source': source,'target': target,'value': value })

        chord = hv.Chord(links)

        # chord = hv.Chord((links, nodes)).select(value=(5, None))
        chord.opts(
            opts.Chord(cmap='Category20', edge_cmap='Category20',
                       labels='index', node_color=dim('index').str()))

        p = hv.render(chord)

        select = Select(title="Option:", value="foo", options=["foo", "bar", "baz", "quux"])

        return(p, select)

    def GetShared(self):
        # Step 1 get the idx where the regions are present
        for i, v in enumerate(list(self.neosData.columns)):
            if v == "Sample":
                start = i + 1
            elif v == "LineID":
                end = i
        # Step 2 Pull columns for the regions (if those regions exist across all samples
        regions = self.neosData.iloc[:, start-1:end].copy()
        regions = regions.replace(-1, np.NaN)
        regionsBySample = [x for _, x in regions.groupby('Sample')]
        regionsBySample = [x.dropna(axis=1) for x in regionsBySample]

        allInfo = {}
        for region in regionsBySample:
            samNam = region.iloc[0,0]
            # Step 3 Set idx of start and end
            samCols = [list(region.columns.values)[i] for i in range(1, region.shape[1])]
            # Step 4 Get every pair of these column names
            regionPairs = [",".join(map(str, comb)) for comb in combinations(samCols, 2)]
            regionSharedDict = dict.fromkeys(regionPairs, 0)

            for pair in regionPairs:
                rowSum = region[[pair.split(',')[0], pair.split(',')[1]]].sum(axis=1)
                regionSharedDict[pair] += rowSum.loc[rowSum==2,].count()

            allInfo.update({samNam:regionSharedDict})

        return(allInfo)

    def GetMatchedNeos(self):
        # Step 1 get the idx where the regions are present
        for i, v in enumerate(list(self.neosData.columns)):
            if v == "Sample":
                start = i + 1
            elif v == "LineID":
                end = i
        # Step 2 Pull columns for the regions (if those regions exist across all samples
        regions = self.neosData.iloc[:, start - 1:end].copy()
        regions = regions.replace(-1, np.NaN)

        regionsBySample = [x for _, x in regions.groupby('Sample')]
        regionsBySample = [x.dropna(axis=1) for x in regionsBySample]

        # Get the number of neoantigen matches
        neoMatches = pd.concat([regions.reset_index(drop=True), self.neosData[['core']].reset_index(drop=True)], axis=1)
        neoMatches = [x for _, x in neoMatches.groupby('core') if x.shape[0] > 1]
        neoMatches = [x.dropna(axis=1) for x in neoMatches]

        edgeListAndValues = {}
        for it in neoMatches:
            allRegions = it.iloc[:, 1 : it.shape[1]-1].copy()
            # Create list of all possible regions that HAVE the core epitope!
            epitopeInTheseSamples = []
            for rowIdx in range(0, allRegions.shape[0]):
                for colIdx in range(0, allRegions.shape[1]):
                    if allRegions.iloc[ rowIdx , colIdx ] == 1:
                        neoNode = it.iloc[ rowIdx , 0 ] + "_" + list(allRegions.columns.values)[colIdx]
                        epitopeInTheseSamples.append(neoNode)

            edges = [",".join(map(str, comb)) for comb in combinations(epitopeInTheseSamples, 2)]
            for edge in edges:
                try:
                    # If Pair Exists Already
                    edgeListAndValues[edge] += 1
                except KeyError:
                    # Add the new edges
                    edgeListAndValues.update({edge : 1})

        return(edgeListAndValues)

    def Histogram(self):
        """ Bokeh alternative histogram """
        # hist, edges = np.histogram(self.neosData[['Affinity']], density=True, bins=50)
        # p = figure(height=400)
        # p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
        #        fill_color="black", line_color="white", alpha=0.5)
        # p.xaxis.axis_label = "Binding Affinity (nM)"
        # p.yaxis.axis_label = "Frequency"

        data = self.neosData.loc[self.neosData['Affinity']<=500].copy()
        df = hv.Dataset(data)
        o = df.hist(dimension='Affinity', groupby='Sample', bins=50, adjoin=False)
        o.opts(opts.Histogram(alpha=0.9, height=200))

        p = hv.render(o)

        return(p)

    def NeoRecoScatter(self):
        # # A vs R colored by RecoPo
        # p = figure(height=400)
        # p.circle(x=x,y=y, size=3, color="navy", alpha=0.5)

        out = self.recopoData.loc[self.recopoData.NeoantigenRecognitionPotential != 0.0].copy()
        # out.NeoantigenRecognitionPotential = np.log10(out.NeoantigenRecognitionPotential)

        p = hv.BoxWhisker(out, 'Sample', 'NeoantigenRecognitionPotential')
        p.opts(logy=True, xrotation=90, height=200, box_fill_color='Sample', show_legend=False, cmap='Category20')

        # p.opts(style=dict(box_color=hv.Cycle('Set1')))
        o = hv.render(p)
        return(o)

    def HeatTable(self):
        ready = {}
        for item in self.sharedCount:
            ready.update({ (item.split(',')[0],item.split(',')[1]): self.sharedCount[item] })

        k = np.array([item for item in ready])
        v = np.array([ready[item] for item in ready])

        unq_keys, key_idx = np.unique(k, return_inverse=True)
        key_idx = key_idx.reshape(-1, 2)
        n = len(unq_keys)
        adj = np.zeros((n, n), dtype=v.dtype)
        adj[key_idx[:, 0], key_idx[:, 1]] = v
        adj += adj.T

        adj = adj.astype(float)
        for i in range(0,adj.shape[0]):
            for k in range(0,adj.shape[1]):
                if k<=i:
                    adj[i,k]=np.nan

        dfout = pd.DataFrame(data=np.log10(adj+0.01),index = unq_keys, columns = unq_keys)
        dfout.index.name = 'Sam1'
        dfout.columns.name = 'Sam2'
        df = pd.DataFrame(dfout.stack(), columns=['Neoantigens']).reset_index()

        source = ColumnDataSource(df)

        import bokeh.palettes as p
        colors = p.Plasma256
        mapper = LinearColorMapper(palette=colors, low=df.Neoantigens.min(), high=df.Neoantigens.max())

        p = figure(title = "log10( Shared Neoantigens )", plot_height=400, plot_width=400, x_range=list(dfout.index), y_range=list(reversed(dfout.columns)),
           toolbar_location=None, tools="", x_axis_location="below")
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "5pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi / 3
        p.rect(x='Sam1',y='Sam2', source=source, width=1, height=1, fill_color={'field':'Neoantigens','transform':mapper}, line_color=None)

        color_bar = ColorBar(color_mapper=mapper, location=(0, 0),
                             ticker=BasicTicker(desired_num_ticks=int(len(colors)/10)))

        p.add_layout(color_bar, 'right')
        return(p)

if __name__=="__main__":
    Options = Parser()
    hv.extension('bokeh')

    # Used for all options
    headerDiv = Div(text="""<h1>NeoPredViz</h1>""", height=15)
    footerDiv = Div(
        text="""<p>NeoPredViz was developed by <a href="https://ryanoschenck.com">Ryan O Schenck</a> for <a href="https://github.com/MathOnco/NeoPredPipe">NeoPredPipe</a></p>""",
        height=15)
    summaryHead = Div(text="""<h5>Summary Table</h5>""", height=5)
    neoshead = Div(text="""<h5>Neoantigens</h5>""", height=5)
    ###

    if not Options.NeoRecoFile and not Options.NeosFile:
        usrData = Data(summaryNeosFile=Options.SummaryFile, allNeosFile=None, recopo=None)
        bar = usrData.SummaryBarChart()
        summaryTable = widgetbox(usrData.SummaryTable())
        usersChildren = [[headerDiv],
                          [bar],
                          [None],
                          [summaryHead],
                          [summaryTable],
                          [footerDiv]]

    elif not Options.NeoRecoFile:
        usrData = Data(summaryNeosFile=Options.SummaryFile, allNeosFile=Options.NeosFile, recopo=None)
        bar = usrData.SummaryBarChart()
        hist = usrData.Histogram()
        heat = usrData.HeatTable() # If NeosFile is provided
        summaryTable = widgetbox(usrData.SummaryTable())
        neoTable = widgetbox(usrData.EpitopeTable())
        usersChildren = [[headerDiv],
                          [bar, hist],
                          [heat, None],
                          [summaryHead],
                          [summaryTable],
                          [neoshead],
                          [neoTable],
                          [footerDiv]]

    elif not Options.NeosFile:
        usrData = Data(summaryNeosFile=Options.SummaryFile, allNeosFile=None, recopo=Options.NeoRecoFile)
        bar = usrData.SummaryBarChart()
        recoBox = usrData.NeoRecoScatter()
        summaryTable = widgetbox(usrData.SummaryTable())
        usersChildren = [[headerDiv],
                          [bar, recoBox],
                          [summaryHead],
                          [summaryTable],
                          [footerDiv]]

    else: # All files provided
        usrData = Data(summaryNeosFile=Options.SummaryFile, allNeosFile=Options.NeosFile, recopo=Options.NeoRecoFile)
        bar = usrData.SummaryBarChart()
        hist = usrData.Histogram()
        heat = usrData.HeatTable() # If NeosFile is provided
        recoBox = usrData.NeoRecoScatter()
        summaryTable = widgetbox(usrData.SummaryTable())
        neoTable = widgetbox(usrData.EpitopeTable())
        usersChildren = [[headerDiv],
                          [bar, hist],
                          [heat, recoBox],
                          [summaryHead],
                          [summaryTable],
                          [neoshead],
                          [neoTable],
                          [footerDiv]]


    mainPageLayout = gridplot(children = usersChildren,
                              sizing_mode="scale_width")

    # Finally yeild outputs
    output_file(Options.Output, title="NeoPredViz")
    show(mainPageLayout)
