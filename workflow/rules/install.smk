rule download_mosaic:
	output: 'programs/MOSAIC/MOSAIC_1.3.7.tar.gz'
	shell:
		"""
		wget https://maths.ucd.ie/~mst/MOSAIC/MOSAIC_1.3.7.tar.gz --directory-prefix programs/MOSAIC
		"""


rule install_mosaic:
	input:
		"programs/MOSAIC/MOSAIC_1.3.7.tar.gz"
	output:
		'programs/MOSAIC/MOSAIC/mosaic.R'
	shell:
		"""
		R CMD INSTALL programs/MOSAIC/MOSAIC_1.3.7.tar.gz

		tar -xvf programs/MOSAIC/MOSAIC_1.3.7.tar.gz --directory programs/MOSAIC

		# patch mosaic to allow setting seed on the cmd line
		patch -b programs/MOSAIC/MOSAIC/mosaic.R ./mosaic.patchfile

		"""
