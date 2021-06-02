rule download_mosaic:
	output: 'programs/MOSAIC/MOSAIC_1.3.7.tar.gz'
	shell:
		"""
		wget https://maths.ucd.ie/~mst/MOSAIC/MOSAIC_1.3.7.tar.gz --directory-prefix programs/MOSAIC
		"""


rule install_mosaic:
	input: "programs/MOSAIC/MOSAIC_1.3.7.tar.gz"
	shell:
		"R CMD INSTALL programs/MOSAIC/MOSAIC_1.3.7.tar.gz"
