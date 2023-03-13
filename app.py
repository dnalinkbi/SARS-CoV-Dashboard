# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import os
import glob
import json
import warnings
import numpy as np
import pandas as pd
import plotly.express as px
from collections import OrderedDict
from dash import Dash, dcc, html, Input, Output, dash_table, State, callback_context
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash_bootstrap_templates import ThemeChangerAIO, template_from_url, load_figure_template

warnings.filterwarnings(action='ignore')


#22_용역23차(PAC)_87 -> 23_87
#23_023_PAC_03 -> 23_023_03

def input_file_list(folder_path):
	#path = "./input_report_files/2022/"
	file_list = glob.glob(folder_path + "/./*xlsx")
	return file_list

def covid_table(file_list):
	# data merge
	merged_data = pd.DataFrame()
	for i in file_list:
		merged_data = pd.concat([merged_data, pd.read_excel(i, header=None).iloc[2:,:]], ignore_index=True)

	# convert column name
	merged_data.columns = 'batch', 'date', 'sample', 'QC', 'result', 'day', 'totalreads', 'mappedreads', 'coverage', 'badbases', 'depth', 'length', 'clade', 'pango'

	# [batch] column rename
	#print(merged_data['batch'].astype('str').str.replace("\(PAC3\)","").str.replace("\(PAC\)","").str.replace("용역","").str.replace("차","").str.split("_", expand = True))
	merged_data['batch'] = merged_data['batch'].astype('str').str.replace("\(PAC3\)","").str.replace("\(PAC\)","").str.replace("용역","").str.replace("차","").str.split("_", expand = True)[1]
	
	# [clade] column rename
	merged_data['clade'] = merged_data['clade'].astype('str').str.replace(" \(Omicron\)", "")

	# dataframe pivot
	merged_data = merged_data[['batch', 'sample', 'QC', 'totalreads', 'mappedreads', 'coverage', 'badbases', 'depth', 'length', 'clade', 'pango']]

	# fill nan value
	merged_data = merged_data.fillna("NA")

	# add pass/fail column
	#merged_data['pf'] = merged_data['coverage'].apply(lambda x: "PASS" if x >= 0.9 else "FAIL")
	merged_data.insert(3, 'P/F', merged_data['coverage'].apply(lambda x: "Pass" if x >= 0.9 else "Fail"))

	# dataframe sort by batch
	merged_data = merged_data.astype({'batch': int}).sort_values(by='batch')
	return merged_data

def QC_table(full_table, qcstate):
	qc_table = full_table[full_table['QC'] == str(qcstate)]
	return qc_table

def groupby_clade(df):
	#df = df.replace({'clade' : {'20C':'etc', '20I (Alpha, V1)':'etc', '21C (Epsilon)':'etc', '21I (Delta)':'etc',
	#							'21K (Omicron)':'etc', '22A (Omicron)':'etc', '22C (Omicron)':'etc', 'recombinant':'etc'}})
	df = df.replace({'clade' : {'20C':'etc', '20I (Alpha, V1)':'etc', '21C (Epsilon)':'etc', '21I (Delta)':'etc',
								'21K':'etc', '22A':'etc', '22C':'etc', 'recombinant':'etc', '20A':'etc', '19A':'etc', '23A':'etc'}})
	df = df.groupby(['batch','clade'], sort=False)['clade'].count().unstack(fill_value=0).stack().to_frame()
	df.reset_index(inplace=True)
	df.columns = 'batch', 'clade', 'count'
	df['percent'] = (df['count'] / df.groupby('batch')['count'].transform('sum')) * 100
	return df

def total_summary(df):
	total_summary = df.agg(sample_count=("sample", "count"), mean_reads=("totalreads", "mean"), median_reads=("totalreads", "median"), max_reads=("totalreads", "max"), min_reads=("totalreads", "min"), 
						   mean_mapped_reads=("mappedreads", "mean"), median_mapped_reads=("mappedreads", "median"), max_mapped_reads=("mappedreads", "max"), min_mapped_reads=("mappedreads", "min"),
						   mean_cov=("coverage", "mean"), max_cov=("coverage", "max"), min_cov=("coverage", "min"), 
						   mean_badbases=("badbases", "mean"), median_badbases=("badbases", "median"), max_badbases=("badbases", "max"), min_badbases=("badbases", "min"),
						   mean_depth=("depth", "mean"), median_depth=("depth", "median"), max_depth=("depth", "max"), min_depth=("depth", "min"))
	return total_summary

def clade_count(df):
	df = groupby_clade(df).groupby('clade')['count'].sum().to_frame()
	df.reset_index(inplace=True)
	df.columns = 'clade', 'count'
	return df





# COVID DATAFRAME
covid_2022 = covid_table(input_file_list("/denovo/workspace.bsy/work/SEQUEL/COVID/final_report/SARS-CoV-Dashboard/input_report_files/2022/"))
covid_2023 = covid_table(input_file_list("/denovo/workspace.bsy/work/SEQUEL/COVID/final_report/SARS-CoV-Dashboard/input_report_files/2023/"))
covid_2022.insert(0, 'year', '2022')
covid_2023.insert(0, 'year', '2023')
covid = pd.concat([covid_2022,covid_2023])
covid_groupby_clade = groupby_clade(covid)
covid_qcgood = QC_table(covid, "good")
covid_qcgood_groupby_clade = groupby_clade(covid_qcgood)
etc_list = ['20C', '20I (Alpha, V1)', '21C (Epsilon)', '21I (Delta)', '21K (Omicron)', '22A (Omicron)', '22C (Omicron)', 'recombinant']


# ===== YEAR FILTER ===== #
year_filter_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Year Select'),
		dbc.RadioItems(
			np.append(['All'], sorted(covid['year'].unique())),
			'All',
			id='filter-year',
		)
	]), className="m-2 shadow"
)

# ===== QC FILTER ===== #
qc_filter_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='QC filter'),
		dbc.Checklist(
			np.append(['All'], sorted(covid['QC'].unique())),
			['All'],
			id='filter-qc',
		)
	]), className="m-2 shadow",
)

# ===== COVERAGE FILTER ===== #
coverage_filter_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Coverage filter'),
		dbc.RadioItems(
			['All', '>=90%', '<90%'],
			'All',
			id='filter-cov',
		)
	]), className="m-2 shadow",
)

# ===== CLADE FILTER ===== #
clade_filter_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Clade Check'),
		dbc.Checklist(
			np.append(['All'], sorted(covid_groupby_clade['clade'].unique())),
			['All'],
			id='check-clade'
		)
	]), className="m-2 shadow",
)

# ===== BATCH RANGE ===== #
batch_range_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Batch Range'),
		dcc.RangeSlider(
			covid_groupby_clade['batch'].min(),
			covid_groupby_clade['batch'].max(),
			step=1,
			value=[covid_groupby_clade['batch'].min(),covid_groupby_clade['batch'].max()],
			allowCross=False,
			marks={str(batch): str(batch) for batch in range(covid_groupby_clade['batch'].min(), covid_groupby_clade['batch'].max()+1, 30)},
			tooltip={"placement": "bottom", "always_visible": True},
			id='batch_slider_value',
		)
	]), className="m-2 shadow",
)

# ===== CLADE DATATABLE ===== #
total_clade_count_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Total Clade Count'),
		dash_table.DataTable(
			data=clade_count(covid).to_dict('records'),
			columns=[{'id':c, 'name':c} for c in clade_count(covid).columns],
			fixed_rows={'headers': True},
			sort_action="native",
			sort_mode='multi',
			editable=True,
		)
	]), className="m-2 shadow"
)

# ===== SUNBRUST PLOT ===== #
sunbrust_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Clade & Pango Sunbrust Plot'),
		dcc.Graph(id='sunburst-graph', style={'height': '45vh'})
	]), className="m-2 shadow"
)

# ===== DATATABLE ===== #
datatable_card = dbc.Card(
	dbc.CardBody([
		html.H4(children='Data Table'),

		dcc.Download(id="download-full"),
		dcc.Download(id="download-filtered"),
		
		dbc.Row([
			dbc.Col(
				html.Div(
					dcc.Dropdown(
						options=[
							{"label": "Excel file", "value": "excel"},
	                    	{"label": "CSV file", "value": "csv"},
	                    ],
						id="download-filtered-dropdown",
						placeholder="Choose download file type."
					),
				), md=4
			),
			dbc.Col(
				html.Div(dbc.Button("Download Total Data", id="btn_full", color="primary"),), md=4
			),
			dbc.Col(
				html.Div(dbc.Button("Download Filtered Data", id="btn_filtered", color="success"),), md=4
			), 
		], className="m-1", align="center"),

		# ===== FILTERED TABLE ===== #
		dash_table.DataTable(
			id='filtered-table',
			columns=[{'id':c, 'name':c} for c in covid.columns],
			fixed_rows={'headers': True},
			sort_action="native",
			sort_mode='multi',
			page_size=20,
			style_table={'overflowX': 'auto'}
		),
	]), className="m-2 shadow"
)

# ===== BAR PLOT ===== #
barplot = dbc.Card(
	dbc.CardBody([
		html.H4(children="Clade BarPlot"),
		dbc.Row([
			# ===== BATCH X COUNT BAR PLOT ===== #
			html.Div([
				dcc.Graph(id='groupby_clade_count', style={'height': '36vh'}),
			]),
			# ===== BATCH X PERCENT BAR PLOT ===== #
			html.Div([
				dcc.Graph(id='groupby_clade_percent', style={'height': '36vh'}),
			])
		])
	]), className="m-2 shadow"
)

# ===== BOX PLOT ===== #
boxplot = dbc.Card(
	dbc.CardBody([
		html.H4(children="Depth & Reads BoxPlot"),
		dbc.Row([
			# ===== DEPTH BOX PLOT ===== #
			html.Div([
				dcc.Graph(id='depth-boxplot', style={'height': '36vh'}),
			]),
			# ===== TOTAL READS BOX PLOT ===== #
			html.Div([
				dcc.Graph(id='reads-boxplot', style={'height': '36vh'}),
			])
		])
	]), className="m-2 shadow"
)

# ===== PARALLEL COORDINATES PLOT ===== #
parallel_coordinates_plot = dbc.Card(
	dbc.CardBody([
		html.H4(children="Parallel Coordinates Plot"),
		dbc.Row([
			html.Div([
				dcc.Graph(id='parallel_coordinates-plot'),
			])
		])
	]), className="m-2 shadow"
)

#===============================#
# ===== DASH APP INSTANCE ===== #
#===============================#
dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY, dbc_css])

# ===== DASH LAYOUT ===== #
app.layout = dbc.Container([
	# ===== Title ===== #
	html.Div([
		html.H1(children='DNALINK SARS-CoV-Analysis DASH BOARD', className="bg-primary bg-gradient text-white text-center fw-bold p-4 m-2"),
	]),
	
	
	html.Div([
		
		# ===== LEFT LAYOUT ===== #
		html.Div(children=[
			# ===== BATCH RANGE ===== #
			batch_range_card,

			
			html.Div([
				dbc.Row([
					dbc.Col([
						# ===== YEAR FILTER ===== #
						year_filter_card,
						# ===== COVERAGE FILTER ===== #
						coverage_filter_card,
						# ===== QC FILTER ===== #
						qc_filter_card,
						# ===== CLADE FILTER ===== #
						clade_filter_card,
						# ===== TOTAL CLADE COUNT TABLE ===== #
						#total_clade_count_card,
					], width=5),
					dbc.Col([
						sunbrust_card,
					], width=7),
					#dcc.Graph(id='indicator')
				], className="g-1")
			]),

			# ===== DATA TABLE ===== #
			#html.Div([
			#	datatable_card
			#]),
		], style={'width': '40%', 'padding': 5, 'flex': 1}),


		# ===== RIGHT LAYOUT =====
		html.Div([
			dbc.Tabs([
				dbc.Tab(label="Parallel Coordinates Plot", children=[parallel_coordinates_plot]),
				dbc.Tab(label="Bar Plot", children=[barplot]),
				dbc.Tab(label="Box Plot", children=[boxplot]),
				dbc.Tab(label="Data Table", children=[datatable_card]),
			]), 
		], style={'width': '50%', 'padding': 13, 'flex': 1})

	], style={'display': 'flex', 'flex-direction': 'row'})
], fluid=True, className="dbc")







# ===== STACKED BAR PLOT CALLBACK ===== #
@app.callback(
	Output('groupby_clade_count', 'figure'),
	Output('groupby_clade_percent', 'figure'),
	#Output('indicator', 'figure'),
	Input('filter-year', 'value'),
	[Input('check-clade', 'value')],
	Input('filter-cov', 'value'),
	Input('filter-qc', 'value'),
	[Input('batch_slider_value', 'value')])
def clade_graph(year, clade, cov, qc, batch):
	# year filter
	if year == 'All':
		df_year = covid
	else:
		df_year = covid[covid['year'] == year]

	# coverage filter
	if cov == '>=90%':
		df_cov = df_year[df_year['coverage'] >= 0.9]
	elif cov == '<90%':
		df_cov = df_year[df_year['coverage'] < 0.9]
	elif cov == 'All':
		df_cov = df_year

	# qc filter
	#if qc == 'All':
	if 'All' in qc:
		df_qc = groupby_clade(df_cov)
	else:
		df_qc = groupby_clade(df_cov[df_cov['QC'].isin(qc)])

	# batch & clade filter
	dfp = df_qc[(df_qc['batch'] >= batch[0]) & (df_qc['batch'] <= batch[1])]
	if 'All' in clade:
		dfc = dfp
	else:
		dfc = dfp[dfp['clade'].isin(clade)]

	#fig_indi = go.Figure(go.Indicator(mode="delta", value=dfc[(dfc['batch'] == batch[1]) & (dfc['clade'] == "22B (Omicron)")]['count'].values[0], delta={"reference": dfc[(dfc['batch'] == batch[0]) & (dfc['clade'] == "22B (Omicron)")]['count'].values[0], "relative": True}))

	# make plot	
	fig_count = px.bar(dfc, x='batch', y='count', color='clade', title='Clade Count BarPlot', color_discrete_map={'22B': 'gray', 'etc': 'yellow'})
	fig_percent = px.bar(dfp, x='batch', y='percent', color='clade', title='Clade Percent BarPlot', color_discrete_map={'22B': 'gray', 'etc': 'yellow'})
	return fig_count, fig_percent



# ===== BOX PLOT & SUNBRUST PLOT CALLBACK ===== #
@app.callback(
	Output('depth-boxplot', 'figure'),
	Output('reads-boxplot', 'figure'),
	Output('filtered-table', 'data'),
	Output("sunburst-graph", "figure"),
	Output("parallel_coordinates-plot", "figure"),
	Input('filter-year', 'value'),
	[Input('check-clade', 'value')],
	Input('filter-cov', 'value'),
	Input('filter-qc', 'value'),
	[Input('batch_slider_value', 'value')])
def boxplot(year, clade, cov, qc, batch):
	# year filter
	if year == 'All':
		df_year = covid
	else:
		df_year = covid[covid['year'] == year]

	# coverage filter
	if cov == '>=90%':
		df_cov = df_year[df_year['coverage'] >= 0.9]
	elif cov == '<90%':
		df_cov = df_year[df_year['coverage'] < 0.9]
	elif cov == 'All':
		df_cov = df_year

	# qc filter
	#if qc == 'All':
	if 'All' in qc: 
		df_qc = df_cov
	else:
		df_qc = df_cov[df_cov['QC'].isin(qc)]

	# clade
	if 'All' in clade:
		dff = df_qc[(df_qc['batch'] >= batch[0]) & (df_qc['batch'] <= batch[1])]
	elif 'etc' in clade:
		dff = df_qc[(df_qc['batch'] >= batch[0]) & (df_qc['batch'] <= batch[1]) & (df_qc['clade'].isin(clade + etc_list))]
	else:
		dff = df_qc[(df_qc['batch'] >= batch[0]) & (df_qc['batch'] <= batch[1]) & (df_qc['clade'].isin(clade))]

	# sunbrust-data-parsing
	dff_sun = dff.groupby(['clade', 'pango'])['sample'].count().to_frame()
	dff_sun.reset_index(inplace=True)
	dff_sun.columns = 'clade', 'pango', 'count'
	

	# parallel-coordinates plot
	dff = dff.replace({'clade' : {'20C':'etc', '20I (Alpha, V1)':'etc', '21C (Epsilon)':'etc', '21I (Delta)':'etc',
								'21K':'etc', '22A':'etc', '22C':'etc', 'recombinant':'etc', '20A':'etc', '19A':'etc', '23A':'etc'}})
	dff_parallel = dff[['batch', 'totalreads', 'coverage', 'depth']]
	dff_parallel.insert(0, 'QC', dff['QC'].map({j: i for i,j in enumerate(sorted(dff['QC'].unique()))}))
	dff_parallel.insert(4, 'clade', dff['clade'].map({j: i for i,j in enumerate(sorted(dff['clade'].unique()))}))
	dff_parallel.insert(0, 'P/F', dff['P/F'].map({j: i for i,j in enumerate(sorted(dff['P/F'].unique()))}))
	


	# make figure
	fig_depth = px.box(dff, x='batch', y='depth', title='Depth BoxPlot')
	fig_reads = px.box(dff, x='batch', y='totalreads', title='Total Reads BoxPlot')
	fig_sun = px.sunburst(dff_sun, path=['clade', 'pango'], values='count', color='clade')
	#fig_parallel = px.parallel_coordinates(dff_parallel, color='batch', color_continuous_scale=px.colors.sequential.Viridis,
	#									   dimensions=['P/F', 'QC', 'totalreads', 'coverage', 'depth', 'clade'])

	fig_parallel = go.Figure(data=
		go.Parcoords(
			line = dict(color=dff_parallel['batch'], colorscale='matter', showscale=True),
			dimensions = list([
				dict(label='Total Reads', range=[dff_parallel['totalreads'].min(), dff_parallel['totalreads'].max()],
					 values=dff_parallel['totalreads']),
				dict(label='Pass/Fail', range=[-1,2], tickvals=[i for i,j in enumerate(sorted(dff['P/F'].unique()))],
					 ticktext=[j for i,j in enumerate(sorted(dff['P/F'].unique()))], values=dff_parallel['P/F']),
				dict(label='Coverage', range=[0, 1],
					 values=dff_parallel['coverage']),
				dict(label='QC', range=[-1,4], 
					 tickvals=[i for i,j in enumerate(sorted(dff['QC'].unique()))],
					 ticktext=[j for i,j in enumerate(sorted(dff['QC'].unique()))], values=dff_parallel['QC']),
				dict(label='Depth', range=[dff_parallel['depth'].min(), dff_parallel['depth'].max()],
					 values=dff_parallel['depth']),
				dict(label='Clade', range=[dff_parallel['clade'].min(), dff_parallel['clade'].max()], 
					 tickvals=[i for i,j in enumerate(sorted(dff['clade'].unique()))],
					 ticktext=[j for i,j in enumerate(sorted(dff['clade'].unique()))], values=dff_parallel['clade']),
			])
		)
	)


	return fig_depth, fig_reads, dff.to_dict('records'), fig_sun, fig_parallel



# ===== FULL DATATABLE DOWNLOAD BTN ===== #
@app.callback(
	Output("download-full", "data"),
	Input("btn_full", "n_clicks"),
	State("download-filtered-dropdown", "value"),
	prevent_initial_call=True,
)
def func(n_clicks_btn, download_type):
	if download_type == "csv":
		return dcc.send_data_frame(covid.to_csv, "SARS-COVID.Full.csv")
	else:
		return dcc.send_data_frame(covid.to_excel, "SARS-COVID.Full.xlsx")


# ===== FILTERED DATATABLE DOWNLOAD BTN ===== #
@app.callback(
	Output("download-filtered", "data"),
	Input("btn_filtered", "n_clicks"),
	State("filtered-table", "derived_virtual_data"),
	State("download-filtered-dropdown", "value"),
	prevent_initial_call=True,
)
def func(n_clicks_btn, df, download_type):
	filtered_df = pd.DataFrame(df)
	if download_type == "csv":
		return dcc.send_data_frame(filtered_df.to_csv, "SARS-COVID.Filtered.csv")
	else:
		return dcc.send_data_frame(filtered_df.to_excel, "SARS-COVID.Filtered.xlsx")





if __name__ == '__main__':
    app.run_server(host='211.174.205.41', port='8060', debug=True, use_reloader=True)
