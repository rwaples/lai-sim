import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# input
flare_report = str(snakemake.input.flare_report)
mosaic_report = str(snakemake.input.mosaic_report)
rfmix_report = str(snakemake.input.rfmix_report)
# output
plot_path = str(snakemake.output.plot_path)


flare_qq = pd.read_csv(flare_report, sep='\t')
rfmix_qq = pd.read_csv(rfmix_report, sep='\t')
mosaic_qq = pd.read_csv(mosaic_report, sep='\t')


def plot_qq_reports(flare, rfmix, mosaic, plot_path, reflines=True):
	with sns.plotting_context('notebook', font_scale=1.3):
		fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 6), sharex=True, sharey='row')
		names = ['flare', 'rfmix', 'mosaic']
		for i, report in enumerate([flare, rfmix, mosaic]):
			name = names[i]
			if reflines:
				point1 = [0, 0]
				point2 = [.1, 0]
				x_values = [point1[0], point2[0]]
				y_values = [point1[1], point2[1]]
				ax[0, i].plot(x_values, y_values, c='k', ls='--')
				point1 = [.9, 1]
				point2 = [1.1, 1]
				x_values = [point1[0], point2[0]]
				y_values = [point1[1], point2[1]]
				ax[0, i].plot(x_values, y_values, c='k', ls='--')
				point1 = [1.9, 2]
				point2 = [2.0, 2]
				x_values = [point1[0], point2[0]]
				y_values = [point1[1], point2[1]]
				ax[0, i].plot(x_values, y_values, c='k', ls='--')
			ax[0, i].scatter(
				(report['bin'] + .5) / len(report) * 2,
				report['mean'],
				c=['r'] + ['b'] * int(np.floor((len(report) - 3) / 2)) + ['r'] + ['b'] * int(np.ceil((len(report) - 3) / 2)) + ['r']
			)
			ax[0, i].set_title(name)
			# ax[0,i].set_xlabel('inferred ancestry dosage bin')
			# ax[0,i].set_ylabel('mean of true ancestry dosage')
			ax[1, i].scatter(
				(report['bin'] + .5) / len(report) * 2,
				report['n'],
				c=['r'] + ['b'] * int(np.floor((len(report) - 3) / 2)) + ['r'] + ['b'] * int(np.ceil((len(report) - 3) / 2)) + ['r']
			)

			ax[1, i].set_xlabel('inferred ancestry dosage bin')
			# ax[1,i].set_ylabel('count')
			ax[1, i].semilogy()
	ax[0, 0].set_ylabel('mean true \nancestry dosage')
	ax[1, 0].set_ylabel('count')
	sns.despine()
	plt.tight_layout()
	plt.savefig(plot_path, dpi=300)


plot_qq_reports(flare_qq, rfmix_qq, mosaic_qq, plot_path)
