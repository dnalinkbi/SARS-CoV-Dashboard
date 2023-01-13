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
from dash_bootstrap_templates import ThemeChangerAIO, template_from_url, load_figure_template

warnings.filterwarnings(action='ignore')


def input_file_list():
	path = "./input_report_files/"
	file_list = glob.glob(path + "./*트.xlsx")
	return file_list

def covid_table(file_list):
	# data merge
	merged_data = pd.DataFrame()
	for i in file_list:
		merged_data = pd.concat([merged_data, pd.read_excel(i, header=None).iloc[2:,:]], ignore_index=True)

	# convert column name 
	merged_data.columns = 'round', 'date', 'sample', 'QC', 'result', 'day', 'totalreads', 'mappedreads', 'coverage', 'badbases', 'depth', 'length', 'clade', 'pango'

	# [round] column rename
	merged_data['round'] = merged_data['round'].astype('str').str.replace("(","").str.replace(")","").str.replace("PAC3","").str.replace("PAC","").str.replace("22_용역","").str.replace("차","").str.split("_", expand = True)[0]

	# dataframe pivot
	merged_data = merged_data[['round', 'sample', 'QC', 'totalreads', 'mappedreads', 'coverage', 'badbases', 'depth', 'length', 'clade', 'pango']]

	# fill nan value
	merged_data = merged_data.fillna("NA")

	# dataframe sort by round
	merged_data = merged_data.astype({'round': int}).sort_values(by='round')
	return merged_data

def QC_table(full_table, qcstate):
	qc_table = full_table[full_table['QC'] == str(qcstate)]
	return qc_table

def groupby_clade(df):
	df = df.replace({'clade' : {'20C':'etc', '20I (Alpha, V1)':'etc', '21C (Epsilon)':'etc', '21I (Delta)':'etc',
								'21K (Omicron)':'etc', '22A (Omicron)':'etc', '22C (Omicron)':'etc', 'recombinant':'etc'}})
	df = df.groupby(['round','clade'], sort=False)['clade'].count().unstack(fill_value=0).stack().to_frame()
	df.reset_index(inplace=True)
	df.columns = 'round', 'clade', 'count'
	df['percent'] = (df['count'] / df.groupby('round')['count'].transform('sum')) * 100
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
covid = covid_table(input_file_list())
covid_groupby_clade = groupby_clade(covid)
covid_qcgood = QC_table(covid, "good")
covid_qcgood_groupby_clade = groupby_clade(covid_qcgood)
etc_list = ['20C', '20I (Alpha, V1)', '21C (Epsilon)', '21I (Delta)', '21K (Omicron)', '22A (Omicron)', '22C (Omicron)', 'recombinant']


# ===== QC FILTER ===== #
qc_filter_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='QC filtering'),
		dcc.RadioItems(
			np.append(['All'], covid['QC'].unique()),
			'All',
			id='filter-qc',
		)
	]), className="mb-3 shadow",
)

# ===== COVERAGE FILTER ===== #
coverage_filter_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='Coverage filtereing'),
		dcc.RadioItems(
			['All', '>=90%', '<90%'],
			'All',
			id='filter-cov',
		)
	]), className="mb-3 shadow",
)

# ===== CLADE FILTER ===== #
clade_filter_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='Clade CheckList'),
		dbc.Checklist(
			np.append(['All'], sorted(covid_groupby_clade['clade'].unique())),
			['All'],
			id='check-clade'
		)
	]), className="mb-3 shadow",
)

# ===== ROUND RANGE ===== #
round_range_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='Round Range'),
		dcc.RangeSlider(
			covid_groupby_clade['round'].min(),
			covid_groupby_clade['round'].max(),
			step=1,
			value=[covid_groupby_clade['round'].min(),covid_groupby_clade['round'].max()],
			allowCross=False,
			marks={str(round): str(round) for round in range(covid_groupby_clade['round'].min(), covid_groupby_clade['round'].max()+1, 30)},
			tooltip={"placement": "bottom", "always_visible": True},
			id='round_slider_value',
		)
	]), className="mb-3 shadow",
)

# ===== CLADE DATATABLE ===== #
total_clade_count_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='Total Clade Count'),
		dash_table.DataTable(
			data=clade_count(covid).to_dict('records'),
			columns=[{'id':c, 'name':c} for c in clade_count(covid).columns],
			fixed_rows={'headers': True},
			sort_action="native",
			sort_mode='multi',
			editable=True,
		)
	]), className="mb-3 shadow"
)

# ===== SUNBRUST PLOT ===== #
sunbrust_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='Clade & Pango Sunbrust Plot'),
		dcc.Graph(id='sunburst-graph')
	]), className="mb-3 shadow"
)

# ===== DATATABLE ===== #
datatable_card = dbc.Card(
	dbc.CardBody([
		html.H2(children='Data Table'),

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
				), width=6
			),
			dbc.Col(
				html.Div(dbc.Button("Download Total Data", id="btn_full", color="primary"),)
			),
			dbc.Col(
				html.Div(dbc.Button("Download Filtered Data", id="btn_filtered", color="success"),)
			)
		], className="mb-2", justify='end'),

		# ===== FILTERED TABLE ===== #
		dash_table.DataTable(
			id='filtered-table',
			columns=[{'id':c, 'name':c} for c in covid.columns],
			fixed_rows={'headers': True},
			sort_action="native",
			sort_mode='multi',
			page_size=10,
			style_table={'overflowX': 'auto'}
		),
	]), className="mb-3 shadow"
)




# ===== DASH APP INSTANCE ===== #
dbc_css = "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
app = Dash(__name__, external_stylesheets=[dbc.themes.SANDSTONE, dbc_css])

# ===== DASH LAYOUT ===== #
app.layout = dbc.Container([
	# ===== Title ===== #
	html.Div([
		html.H1(children='DNA LINK SARS-CoV-Analysis DASH BOARD', className="bg-primary text-white text-center fw-bold"),
	]),
	
	html.Div([
		# ===== LEFT LAYOUT ===== #
		html.Div(children=[
			
			# ===== QC FILTER ===== #
			qc_filter_card,
			
			# ===== COVERAGE FILTER ===== #
			coverage_filter_card,
			
			# ===== CLADE FILTER ===== #
			clade_filter_card,

			# ===== ROUND RANGE ===== #
			round_range_card,
			
			# ===== CLADE DATATABLE & SUNBRUST PLOT ===== #
			html.Div([
				dbc.Row([
					dbc.Col(
						total_clade_count_card, width=4
					),
					dbc.Col(
						sunbrust_card,					
					)
				])
			]),

			# ===== DATA TABLE ===== #
			html.Div([
				datatable_card
			])
		], style={'width': '35%', 'padding': 5, 'margin': 20, 'flex': 1}),


		# ===== RIGHT LAYOUT =====
		html.Div(children=[

			# ===== ROUND X COUNT BAR PLOT ===== #
			html.Div([
				dcc.Graph(id='groupby_clade_count', style={'height': '40vh'}),
			]),

			# ===== ROUND X PERCENT BAR PLOT ===== #
			html.Div([
				dcc.Graph(id='groupby_clade_percent', style={'height': '40vh'}),
			]),

			# ===== DEPTH BOX PLOT ===== #
			html.Div([
				dcc.Graph(id='depth-boxplot', style={'height': '40vh'}),
			]),

			# ===== TOTAL READS BOX PLOT ===== #
			html.Div([
				dcc.Graph(id='reads-boxplot', style={'height': '40vh'}),
			]),

		], style={'width': '50%', 'padding': 5, 'flex': 1})
	], style={'display': 'flex', 'flex-direction': 'row'})
], fluid=True, className="dbc")







# ===== STACKED BAR PLOT CALLBACK ===== #
@app.callback(
	Output('groupby_clade_count', 'figure'),
	Output('groupby_clade_percent', 'figure'),
	[Input('check-clade', 'value')],
	Input('filter-cov', 'value'),
	Input('filter-qc', 'value'),
	[Input('round_slider_value', 'value')])
def clade_graph(clade, cov, qc, round):
	# coverage filter
	if cov == '>=90%':
		df_cov = covid[covid['coverage'] >= 0.9]
	elif cov == '<90%':
		df_cov = covid[covid['coverage'] < 0.9]
	elif cov == 'All':
		df_cov = covid

	# qc filter
	if qc == 'All':
		df_qc = groupby_clade(df_cov)
	else:
		df_qc = groupby_clade(QC_table(df_cov, qc))

	# round & clade filter
	dfp = df_qc[(df_qc['round'] >= round[0]) & (df_qc['round'] <= round[1])]
	if 'All' in clade:
		dfc = dfp
	else:
		dfc = dfp[dfp['clade'].isin(clade)]
	
	# make plot	
	fig_count = px.bar(dfc, x='round', y='count', color='clade', color_discrete_map={'22B (Omicron)': 'gray', 'etc': 'yellow'})
	#fig_count = px.line(dfc, x='round', y='count', color='clade')
	fig_percent = px.bar(dfp, x='round', y='percent', color='clade', color_discrete_map={'22B (Omicron)': 'gray', 'etc': 'yellow'})
	return fig_count, fig_percent



# ===== BOX PLOT & SUNBRUST PLOT CALLBACK ===== #
@app.callback(
	Output('depth-boxplot', 'figure'),
	Output('reads-boxplot', 'figure'),
	Output('filtered-table', 'data'),
	Output("sunburst-graph", "figure"),
	[Input('check-clade', 'value')],
	Input('filter-cov', 'value'),
	Input('filter-qc', 'value'),
	[Input('round_slider_value', 'value')])
def boxplot(clade, cov, qc, round):
	# coverage filter
	if cov == '>=90%':
		df_cov = covid[covid['coverage'] >= 0.9]
	elif cov == '<90%':
		df_cov = covid[covid['coverage'] < 0.9]
	elif cov == 'All':
		df_cov = covid

	# qc filter
	if qc =='All':
		df_qc = df_cov
	else:
		df_qc = QC_table(df_cov, qc)

	# clade
	if 'All' in clade:
		dff = df_qc[(df_qc['round'] >= round[0]) & (df_qc['round'] <= round[1])]
	elif 'etc' in clade:
		dff = df_qc[(df_qc['round'] >= round[0]) & (df_qc['round'] <= round[1]) & (df_qc['clade'].isin(clade + etc_list))]
	else:
		dff = df_qc[(df_qc['round'] >= round[0]) & (df_qc['round'] <= round[1]) & (df_qc['clade'].isin(clade))]

	# sunbrust-data-parsing
	dff_sun = dff.groupby(['clade', 'pango'])['sample'].count().to_frame()
	dff_sun.reset_index(inplace=True)
	dff_sun.columns = 'clade', 'pango', 'count'
	
	# make figure
	fig_depth = px.box(dff, x='round', y='depth')
	fig_reads = px.box(dff, x='round', y='totalreads')
	fig_sun = px.sunburst(dff_sun, path=['clade', 'pango'], values='count', color='clade')
	
	return fig_depth, fig_reads, dff.to_dict('records'), fig_sun



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
    app.run_server(host='211.174.205.41', port='8050', debug=True)
